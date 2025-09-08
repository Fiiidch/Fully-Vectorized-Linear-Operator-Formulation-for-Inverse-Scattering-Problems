include("data_reader.jl")
using HDF5

"""
Home of function pairs create_***_measuring_vector() and evaluate_***_matrix(). 
When the evaluate... funcitons are called, they load their corresponding inverted matrices from object data.
They call the their corresponding create... functions to form all the defect measurements and refernce measurements into a vector.

The Matrix is evaluated through matrix multiplication, and the results saved as h5 files. 

These functions are only called within the evaluate_matrix() funciton from command_center.jl

To understand the meaning of the naming "single" "full" and "reduced" by which the functions are named, see functions create_***_system_matrix() command_center.jl
for a description.
"""

function create_single_measuring_vector(object_name::String, defect_set::String, matrix_name::String)

    # Open the matrix file to get emitter_idx and sensor_idx_array
    matrix_file = "data/$object_name/Inverted_Matrices/$matrix_name.h5"
    h5f = h5open(matrix_file, "r")

    emitter_idx = read(h5f, "emitter_idx")
    sensor_idx_array = read(h5f, "sensor_idx_array")
    if is_debugmode(object_name)
        println("called with sensor array", sensor_idx_array, "defect_set ", defect_set, "emitter_idx ", emitter_idx)
        println("creating measuring vector for $matrix_name")
    end
    close(h5f)
    # For each sensor_idx, get the source-defect difference vector
    diffs = Vector{Any}()
    for sensor_idx in sensor_idx_array
        push!(diffs, get_source_defect_diff(object_name, emitter_idx, sensor_idx, defect_set))
    end
    measuring_vector = vcat(diffs...)

    return measuring_vector

end

function evaluate_single_matrix(object_name::String, inverted_matrix::String, defect_set::String)


    results_dir = "data/$object_name/results/"
    if !isdir(results_dir)
        mkdir(results_dir)
    end

    results_file = "$results_dir/results_$(inverted_matrix)_$defect_set.h5"
    if isfile(results_file)
        rm(results_file)
    end

    matrix_file = "data/$object_name/Inverted_Matrices/$inverted_matrix.h5"
    h5file = h5open(matrix_file, "r")
    sensor_idx_array = read(h5file, "sensor_idx_array")
    lambda_keys = filter(key -> startswith(key, "Lambda_"), keys(h5file))
    measuring_vec = create_single_measuring_vector(object_name, defect_set, inverted_matrix)
    for (k,key) in enumerate(lambda_keys)
        h5file = h5open(matrix_file, "r")
        
        matrix = read(h5file, key)
    
        close(h5file)
        if k == 1
            println("Size of matrix to be evaluated: ", size(matrix))
        end
        
        
        t_start = time()
        sol = matrix * measuring_vec
        t_end = time()
        if is_debugmode(object_name)
            println("Matrix multiplication for $key took $(t_end - t_start) seconds.")
        end
           
        h5write(results_file, key, sol)
    end

    h5file = h5open(matrix_file, "r")
    defect_projector = read(h5file, "defect_projector")
    h5write(results_file, "defect_projector", defect_projector)
    h5write(results_file, "defect_set", defect_set)
    h5write(results_file, "defect_name", defect_set)



    close(h5file)

end



function create_full_measuring_vector(object_name::String, defect_set::String, matrix_name::String)
      # Open the matrix file to get emitter_idx and sensor_idx_array
    matrix_file = "data/$object_name/Inverted_Matrices/$matrix_name.h5"

    sensor_idx_array = get_sensors(object_name)


    diffs = [get_source_defect_diff(object_name, e_idx, s_idx, defect_set) 
             for (s_idx,_) in enumerate(sensor_idx_array), (e_idx,_) in enumerate(sensor_idx_array)]
    measuring_vector = vcat(diffs...)
    return measuring_vector
end

function evaluate_full_matrix(object_name::String, inverted_matrix::String, defect_set::String)

    results_dir = "data/$object_name/results/"
    if !isdir(results_dir)
        mkdir(results_dir)
    end

    results_file = "$results_dir/results_$(inverted_matrix)_$defect_set.h5"
    if isfile(results_file)
        rm(results_file)
    end

    matrix_file = "data/$object_name/Inverted_Matrices/$inverted_matrix.h5"
    h5file = h5open(matrix_file, "r")
    lambda_keys = filter(key -> startswith(key, "Lambda_"), keys(h5file))
    measuring_vec = create_full_measuring_vector(object_name, defect_set, inverted_matrix)

      for (k,key) in enumerate(lambda_keys)
        h5file = h5open(matrix_file, "r")
        
        matrix = read(h5file, key)
    
        close(h5file)
        if k == 1
            println("Size of matrix to be evaluated: ", size(matrix))
        end
        
        t_start = time()
        sol = matrix * measuring_vec
        t_end = time()
        if is_debugmode(object_name)
            println("Matrix multiplication for $key took $(t_end - t_start) seconds.")
        end
           
        h5write(results_file, key, sol)
    end

    h5file = h5open(matrix_file, "r")
    if haskey(h5file, "defect_projector")
        defect_projector = read(h5file, "defect_projector")
        h5write(results_file, "defect_projector", defect_projector)
    end

    
    h5write(results_file, "defect_name", defect_set)
    h5write(results_file, "defect_set", defect_set)

    close(h5file)
end



function create_reduced_measuring_vector(object_name::String, defect_set::String, matrix_name::String)
    matrix_file = "data/$object_name/Inverted_Matrices/$matrix_name.h5"
    h5f = h5open(matrix_file, "r")
    emitter_idx_array = read(h5f, "emitter_idx_array")
    sensor_idx_array_per_emitter = read(h5f, "sensor_idx_array_per_emitter")
    close(h5f)

    # sensor_idx_array_per_emitter is a 2D array (padded with -1)
    diffs = Vector{Any}()
    for e in emitter_idx_array
        sensors = sensor_idx_array_per_emitter[:, e]
        # Remove padding (-1)
        sensors = sensors[sensors .!= -1]
        for (s) in sensors
            push!(diffs, get_source_defect_diff(object_name, e, s, defect_set))
        end
    end
    measuring_vector = vcat(diffs...)
    return measuring_vector
end

function evaluate_reduced_matrix(object_name::String, inverted_matrix::String, defect_set::String)
    
    results_dir = "data/$object_name/results/"
    if !isdir(results_dir)
        mkdir(results_dir)
    end

    results_file = "$results_dir/results_$(inverted_matrix)_$defect_set.h5"
    if isfile(results_file)
        rm(results_file)
    end

    matrix_file = "data/$object_name/Inverted_Matrices/$inverted_matrix.h5"
    h5file = h5open(matrix_file, "r")
    lambda_keys = filter(key -> startswith(key, "Lambda_"), keys(h5file))
    close(h5file)

    measuring_vec = create_reduced_measuring_vector(object_name, defect_set, inverted_matrix)

    for (k,key) in enumerate(lambda_keys)
        h5file = h5open(matrix_file, "r")
        matrix = read(h5file, key)
        close(h5file)
        if k == 1
            println("Size of matrix to be evaluated: ", size(matrix))
        end

        t_start = time()
        sol = matrix * measuring_vec
        t_end = time()
        if is_debugmode(object_name)
            println("Matrix multiplication for $key took $(t_end - t_start) seconds.")
        end
           
        h5write(results_file, key, sol)
    end

    h5file = h5open(matrix_file, "r")
    if haskey(h5file, "defect_projector")
        defect_projector = read(h5file, "defect_projector")
        h5write(results_file, "defect_projector", defect_projector)
    end
    h5write(results_file, "defect_name", defect_set)
    h5write(results_file, "defect_set", defect_set)
    close(h5file)

end