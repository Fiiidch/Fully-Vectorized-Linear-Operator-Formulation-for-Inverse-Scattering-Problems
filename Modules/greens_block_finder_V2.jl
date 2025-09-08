using HDF5
using LinearAlgebra
using FFTW


"""
    Calculates the impulse response between n sets of sources and signals using Tikhonov regularization.
    This implementation is very efficient computation wise. However, in the current iteration, it has to load 
    all the raw data into memory, which is inefficient memory usage. 

    This can be circumvented by adjusting the code to instead just save which data in fed in, and loading the 
    neccesary parts of the indexes h5 files as needed.

    Quck mathematical explanation:
    Computes the impulse response via Tikhonov regularization by minimizing 
    sum_i |d_i - s_i * g|^2 + lambda * |g|^2 for each frequency and space point.
    Assumes signals are real in time domain, so conjugates in the derivation cancel out.
    This gives a closed-form solution without matrix inversion, applied per (freq, space).

    Resulting expression: g = sum_i(d_i * s_i) / (sum_i(s_i^2) + lambda)
    where g, s, d implicitly additionally depend on f, n for frequency, space.

    # Here:
      - i: noise index, numbering of 
      - d_i: measured data in frequency domain
      - s_i: source signal in frequency domain
      - g: unknown impulse response at a given (freq, space) point
      - lambda: regularization strength
    Assumes original time-domain signals are real, so conjugates cancel in derivation.
    Closed-form result: g = sum_i(d_i * s_i) / (sum_i(s_i^2) + lambda)
"""
mutable struct GreensFinder
    data::Array{Array{Complex{Float32}, 2}, 1}
    source::Array{Array{Complex{Float32}, 1}, 1}
    lambda::Float32
    n_t::Int
    n_space::Int

    # Inner constructor
    function GreensFinder(data, source, lambda, n_t=0, n_space=0)
        new(data, source, lambda,  n_t, n_space)
    end
end


function GreensFinder(lambda::AbstractFloat) 
    source = Array{Array{Complex{Float32}, 1}, 1}()
    data = Array{Array{Complex{Float32}, 2}, 1}()
    return GreensFinder(data, source, lambda)
end


function Add_source_data!(GF::GreensFinder, source::Array{ComplexF32, 1}, data::Array{ComplexF32, 2})
    # The additional data is for different measurements for noise, not different emitters.
    @assert size(source,1) == size(data,2) "Source and data must have the same length in the third dimension. Cuurent Sizes: $(size(source)) and $(size(data)[3])"

    if GF.n_space == 0 && GF.n_t == 0
        GF.n_space, GF.n_t = size(data)
    end

    @assert (GF.n_space, GF.n_t) == size(data) "Data must have the same size as the first data added"

    push!(GF.source, source)
    push!(GF.data, data)
end

function find_response(GF::GreensFinder)
    # Calculate the response using Tikhonov regularization.
    # This response finder is hard coded for L = Identity Matrix in frequency domain.
    
    data_count = length(GF.data)
    completed_loops = 0
    total_loops = GF.n_space
    last_update = time()
    n_f = length(GF.source[1])
    g_block = zeros(ComplexF32, GF.n_space, n_f)

    #This for loop directly calculates the Tikhonov regularization response index wise.
    #This is vastly more efficient than the previous method, which used a matrix inversion. (3.5hrs -> 5s)
    Threads.@threads for n in 1:GF.n_space
        for t in 1:n_f 
            #The bottom one is technically correct, but makes crustly results. the top one should not work, but does.
            g_block[n, t] = sum(GF.data[i][n, t] * GF.source[i][t] for i in 1:data_count) / (sum(GF.source[i][t]^2 for i in 1:data_count) + GF.lambda)
            #g_block[n, t] = sum(GF.data[i][n, t] * conj(GF.source[i][t]) for i in 1:data_count) / (sum((GF.source[i][t])*conj(GF.source[i][t]) for i in 1:data_count) + GF.lambda)
        end
        completed_loops += 1
        if time() - last_update >= 10
            println("Completed loops: $completed_loops / $total_loops")
            last_update = time()
        end
    end
    return g_block
end