using FFTW
include("data_reader.jl")
using HDF5
using LinearAlgebra
using SparseArrays
"""
    This function is where the entire Program spends most of it's time.
    Computes the Tikhonov-regularized inverse operator (regularized pseudoinverse):

        (V'V + λL'L)^(-1) V'

    This operator can be applied to a data vector `b` to obtain the Tikhonov-regularized solution `x` to the least squares problem:

        x = argmin_x ||Vx - b||^2 + λ||Lx||^2

    where:
      - `V` is the system matrix,
      - `L` is the regularization matrix (default: identity),
      - `λ` is the regularization parameter,
      - `b` is the data vector,
      - `x` is the solution vector that balances data fit and regularization.

    Returns the operator matrix, not the solution itself.
    
    This function should be multithreaded by default due to internal use of LAPLACK and BLAS methods. 
    You may want to double check on your system, I am not very sure about the details of this, and have not managed to make it run faster.

    You can call the amount of cores with LinearAlgebra.BLAS.get_num_threads() and set with LinearAlgebra.BLAS.set_num_threads(n)

    Note: Running it with 16 Bit Floats makes it run 10-20 times slower due to CPU specialisation.
    The F32 implementation make it it roughly 2x slower compared to F64, but avoids the Memory bottleneck. 
"""
function _tikhonov_inverse(V::Array{ComplexF32}; L::Union{Array{ComplexF32}, UniformScaling{Bool}, SparseMatrixCSC{Float32}} = I, lambda::Float32, debug_mode::Bool = false)
    if debug_mode
        println("Size of V in memory: ", Base.summarysize(V) / 1024^2, " MB")
    end

    start_time = time()
    X = (V)' * V

    if L isa UniformScaling
        @inbounds for i in 1:size(X, 1)
            X[i, i] += lambda
        end
    elseif L isa Diagonal
        for i in 1:size(X, 1)
            X[i, i] += L[i,i] * lambda
        end
    else
        X .+= lambda .* (L' * L)
    end
    
    L = nothing

    if debug_mode
        println("Time taken for Matmul computation: (X = V' * V + lambda * L' * L) ", time() - start_time, " seconds")
    end

    start_time = time()
    Res = X \ (V)'
    if debug_mode
        println("Time taken for Tikhonov Inversion: ", time() - start_time, " seconds")
    end
    return Res
end

function Tikhonov_invert_matrices(object_name::String, matrix::Array{ComplexF32,2}; lambdas = [0.1, 0.01, 0.001], file_name::String = "Inverted_matrices.h5", debug_mode = false)

    # Calculate the Tikhonov-regularized inverse operator for each of the lambda values. 

    target_dir = "data/$object_name/Inverted_Matrices"

    lambdas = convert(Array{Float32}, lambdas)
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    target_file = target_dir * "/" * file_name
    if !isfile(target_file)
        h5open(target_file, "w") do h5f
            println("Created new HDF5 file: $target_file")
        end
    end
    h5file = h5open(target_file, "r+")
    for lambda in lambdas
        dataset_name = "Lambda_$lambda"
        if isfile(target_file)
            if haskey(h5file, dataset_name)
                println("Dataset $dataset_name already exists in $target_file. Skipping...")
                continue
            end
        end
    
        println("Calculating Tikhonov inverse with lambda = $lambda")
        # Create a sparse diagonal matrix with 1s except at sensor positions

        n_space = size(get_xyz(object_name), 1)
        sensors = get_sensors(object_name)
                
        # Get matrix size from the actual matrix
        _, matrix_width = size(matrix)
        defect_dimensionality = matrix_width ÷ n_space
                
        # Create expanded diagonal elements
        diagonal_elements = ones(Float32, matrix_width)
                
        # For each sensor position, set all corresponding dimensional elements to zero
        for sensor in sensors
            for d in 1:defect_dimensionality
                idx = (sensor-1) * defect_dimensionality + d
                diagonal_elements[idx] = 0.0
            end
        end
                
        L = sparse(Diagonal(diagonal_elements))

        if debug_mode
            println("Created regularization matrix L with $(nnz(L)) non-zeros")
            println("Sensor positions set to zero: $(length(sensors))")
        end

        timing = @timed mat = _tikhonov_inverse(matrix, L = L, lambda = lambda, debug_mode = debug_mode )
        println("Time taken for Tikhonov inverse: ", timing.time, " seconds")
        h5write(target_file, dataset_name, mat)
    end
end

function get_system_matrix_1pulse(object_name::String, emitter_idx::Int, sensor_idx_array::Array{Int}, defect_projector; Greens_block::String = "nothing")::Array{ComplexF32}
    # Calculates the System matrix as defined in the script. The Diagonal property of certain submatrix blocks 
    # in frequency domain have been used, to optimize expression: submatrix_block[:, i] = G1[n,:] .* G2[n,:] .* defect_projector_f[:, i]
    # This seems to be decently faster than using SparseArrays for some reason. (5s vs 3.5hrs)


    if ndims(defect_projector) == 1
        defect_projector = reshape(Float32.(defect_projector), :, 1)
    else
        defect_projector = Float32.(defect_projector)
    end

    defect_projector_f = rfft(defect_projector, 1)
    defect_projector_f = defect_projector_f[1:(get_freq_range(object_name)),:]
    defect_projector_f[1:(get_low_freq_idx(object_name)-1)] .= 0
    n_f , defect_dimensionality = size(defect_projector_f)
    @assert n_f == get_freq_range(object_name)

    xyz = get_xyz(object_name)
    n_space = size(xyz, 1)
    sensors = get_sensors(object_name)
    sensors = sensors[sensor_idx_array]
    if is_debugmode(object_name) == true
        println("creating_matrix with sensors ", sensors, " indexes ", sensor_idx_array)
    end
    n_sensors = length(sensors)
    pulse_f = get_source(object_name, emitter_idx)
    
    full_matrix = zeros(ComplexF32, n_sensors * n_f, n_space * defect_dimensionality)
    submatrix_block = zeros(ComplexF32, n_f, defect_dimensionality)
    start_time = time()


    greens_block_pulse = get_greens_block(object_name, emitter_idx, Greens_block = Greens_block)
    
    G2 = _get_field_times_source(greens_block_pulse, pulse_f)

    for (m, m_idx) in enumerate(sensor_idx_array)
       G1 = get_greens_block(object_name, m_idx, Greens_block = Greens_block)
       Threads.@threads for n in 1:n_space
            for i in 1:defect_dimensionality
                submatrix_block[:, i] = G1[n,:] .* G2[n,:] .* defect_projector_f[:, i]
            end
            # now fill the full_matrix with the submatrix_block. Doing this directly would likely increase efficiency significantly
            full_matrix[(m - 1) * n_f + 1:m * n_f, (n - 1) * defect_dimensionality + 1:n * defect_dimensionality] = submatrix_block
        end
    end
    if is_debugmode(object_name)
        println("Matrix generation time: ", time() - start_time, " seconds")
        println("Matrix size: ", size(full_matrix))
    end
    return full_matrix
end