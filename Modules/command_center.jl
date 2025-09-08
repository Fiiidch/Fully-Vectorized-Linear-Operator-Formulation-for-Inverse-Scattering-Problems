using HDF5
using FFTW
include("custom_plotters.jl")
include("data_reader.jl")
include("greens_block_finder_V2.jl")
include("MatrixGenerator.jl")
include("evaluate_inverted_matrix.jl")
"""
This module is for creating and managing the object data files for the user.
Interaction from user should be only done through this module.
"""

"""
Initializes a new object with geometry, sensor info, and parameters. Creates all required folders and saves metadata.
"""
function initialize_new_object(object_name::String, sensors::Vector{Int}, XYZ::AbstractArray{<:AbstractFloat, 2},
    low_freq::Real, max_freq::Real, n_t::Int, data_freq::Real; debug_mode::Bool = false)
    if size(XYZ, 2) != 3
        error("The second dimension of XYZ must be 3. Current size: $(size(XYZ, 2))")
    end

    if max_freq * 2 > data_freq
        @warn "max_freq ($max_freq) is larger than allowed: data_freq ($data_freq / 2). Setting max_freq to data_freq / 2 ($(data_freq / 2))."
        max_freq = data_freq / 2
    end

    if max_freq < 0 || low_freq < 0 || low_freq >= max_freq
        error("Bandpass frequencies max_freq ($max_freq) and low_freq ($low_freq) are silly.")
    end

    if !isdir("data")
        mkdir("data")
    end

    if isdir("data/$(object_name)")
        error("Directory 'data/$(object_name)' already exists. Delete object using delete_object(object_name) or choose a different object name.")
    end

    mkdir("data/$(object_name)")
    mkdir("data/$(object_name)/raw_data")
    mkdir("data/$(object_name)/defect_measurements")

    if !isdir("Plots")
        mkdir("Plots")
    end

    h5write("data/$(object_name)/object_data.h5", "XYZ", convert(Array{Float32, 2}, XYZ))
    h5write("data/$(object_name)/object_data.h5", "sensors", convert(Array{Int64}, sensors))
    h5write("data/$(object_name)/object_data.h5", "low_freq", convert(Float32, low_freq)) # Low frequency of bandpass filter. Yes, the name is misleading
    h5write("data/$(object_name)/object_data.h5", "max_freq", convert(Float32, max_freq)) # High bar of bandpass filter
    h5write("data/$(object_name)/object_data.h5", "n_t", convert(Int64, n_t)) # Total time samples of input data
    h5write("data/$(object_name)/object_data.h5", "data_freq", convert(Float32, data_freq)) # Sampling frequency of input data
    h5write("data/$(object_name)/object_data.h5", "debug_mode", debug_mode) # Sampling frequency of input data

    println("Object $object_name created.")
    println()
end

"""
Adds a raw measurement (time domain) for a given sensor to the object. Converts to frequency domain and saves.
"""
function add_raw_data(object_name::String, raw_data::AbstractArray{<:AbstractFloat, 2}, sensor_location::Int)
    
    file_path = "data/$(object_name)/object_data.h5"
    if !isfile(file_path)
        error("File $(file_path) does not exist. Object may not be initialized.")
    end

    sensors = h5read(file_path, "sensors")
    if sensor_location ∉ sensors
        error("Sensor $(sensor) does not exist in the defined sensors array. Cannot add raw data.")
    end

    n_t = get_n_t(object_name)
    if size(raw_data, 2) != n_t
        error("Raw data time dimension ($(size(raw_data, 2))) does not match expected n_t ($n_t).")
    end

    sensor_index = findfirst(==(sensor_location), sensors)

    raw_data_file = "data/$(object_name)/raw_data/raw_data_from_sensor_$(sensor_index).h5"

    if isfile(raw_data_file)
        raw_h5file = h5open(raw_data_file)
        set_index = length(keys(raw_h5file)) + 1
        close(raw_h5file)
    else
        set_index = 1
    end
    
    raw_data = convert(Array{Float32}, raw_data)
    max_freq = get_max_freq(object_name)
    data_freq = get_data_freq(object_name)
    freqs = FFTW.rfftfreq(n_t, data_freq)
    raw_data_f = rfft(raw_data, 2)

    max_freq_idx = findfirst(freqs .>= max_freq)
    low_freq_idx = get_low_freq_idx(object_name)

    raw_data_f[:, 1:low_freq_idx-1] .= 0
    raw_data_f = raw_data_f[:, 1:(max_freq_idx-1)]

    n, f = size(raw_data_f)

    h5file = h5open("data/$(object_name)/object_data.h5", "r")
    if "freq_range" in keys(h5file)
        existing_f = read(h5file, "freq_range")
        close(h5file)
        if existing_f != f
            error("New frequency range ($f) does not match existing range ($existing_f). The data likely has differing time samples.")
        end
    else
        close(h5file)
        h5write("data/$(object_name)/object_data.h5", "freq_range", f)
    end

    h5write(raw_data_file, "$(set_index)", raw_data_f)
    
end

"""
Adds a defect measurement (time domain) for a given emitter/sensor pair and defect group. Converts to frequency domain and saves.
"""
function add_measurements(object_name::String, defect_group::String, sensor_location::Int, measurement_location::Int, 
                            measurements::Vector)
    file_path = "data/$(object_name)/object_data.h5"
    sensors = h5read(file_path, "sensors")

    if sensor_location ∉ sensors
        error("Sensor $(sensor_location) does not exist in the sensors array.")
    end

    if measurement_location ∉ sensors
        error("Sensor $(measurement_location) does not exist in the sensors array.")
    end

    n_t = get_n_t(object_name)
    if length(measurements) != n_t
        error("Measurement time dimension ($(length(measurements))) does not match expected n_t ($n_t)")
    end

    target_dir = "data/$(object_name)/defect_measurements/$defect_group"
    if !isdir(target_dir)
        mkdir(target_dir)
    end

    measurements = convert(Array{Float32}, measurements)
    max_freq = get_max_freq(object_name)
    data_freq = get_data_freq(object_name)

    freqs = FFTW.rfftfreq(n_t, data_freq)
    measurements_f = rfft(measurements)
    max_freq_idx = findfirst(freqs .>= max_freq)
    low_freq_idx = get_low_freq_idx(object_name)
    measurements_f[1:low_freq_idx-1] .= 0
    measurements_f = measurements_f[1:max_freq_idx-1]

    sensor_origin_index = findfirst(==(sensor_location), sensors)
    measurement_location_index = findfirst(==(measurement_location), sensors)

    measurement_data_file_location = "data/$(object_name)/defect_measurements/$defect_group/measurements_from_emitter_$(sensor_origin_index).h5"

    if !isfile(measurement_data_file_location)
        # This is a workaround to initiate an empty h5 file.
        h5write(measurement_data_file_location, "init", 0)
    end

    measurement_data_file = h5open(measurement_data_file_location, "r+")
    
    local group_name = "sensor_$(measurement_location_index)"
    if !(group_name in keys(measurement_data_file))
        create_group(measurement_data_file, group_name)
    end

    set_index = length(keys(measurement_data_file[group_name])) + 1
    h5write(measurement_data_file_location, "$group_name/$set_index", measurements_f)

    if is_debugmode(object_name)
        println("Measurement data added to $measurement_data_file_location under group $group_name/$set_index '$defect_group'. With parameters: emitter $sensor_location, sensor $measurement_location.")
    end
    close(measurement_data_file)
end

"""
Deletes all data and folders for the specified object.
"""
function delete_object(object_name::String)
    if isdir("data/$(object_name)")
        rm("data/$(object_name)", recursive=true)
        println("Deleted object: $(object_name)")
    else
        println("Object '$(object_name)' does not exist.")
    end
    println()
end

"""
Calculates and saves Green's functions (impulse responses) for all sensors using Tikhonov regularization.
"""
function calculate_greens_function(object_name::String; block_name::String = "Greens_block_$(object_name)_L_0.h5", lambda = 0.0000001)
    println("Starting impulse response calculation for object '$object_name'...")

    t_start = time()
    _find_source_pulses(object_name)

    raw_data_path = "data/$object_name/raw_data"
    if !isdir(raw_data_path) || isempty(readdir(raw_data_path))
        throw(ErrorException("Raw data folder is empty or doesn't exist. Please add raw data before calculating Green's function."))
    end
    greens_dir = "data/$object_name/Greens_Blocks"
    
    if !isdir(greens_dir)
        mkdir(greens_dir)
    end
    greens_block_path = greens_dir * "/" * block_name
    if isfile(greens_block_path)
        println("Warning: File $greens_block_path already exists. Overwriting...")
        rm(greens_block_path)
    end

    sensors = get_sensors(object_name)

    for (i, sensor) in enumerate(sensors)
        raw_data_file = "data/$object_name/raw_data/raw_data_from_sensor_$i.h5"
        if !isfile(raw_data_file)
            continue
        end
        h5file = _open_h5file(raw_data_file)
        sets = keys(h5file)
        finder = GreensFinder(lambda)

        for (s, set_number) in enumerate(sets)
            Add_source_data!(finder, get_single_location_from_raw(object_name, i, sensor, s), get_raw_data(object_name, i, s))
        end
        
        block_response = find_response(finder)
        block_response[:, 1:get_low_freq_idx(object_name)-1] .= 0 
        if is_debugmode(object_name)
            println("Calculated Green's functions for sensor $i (sensor id: $sensor) and saved to $greens_dir/$block_name")
        end

        h5write(greens_dir*"/$block_name", "sensor_$i", block_response)

        finder = nothing
    end
    t_end = time()
    println("Finished impulse response calculation for '$block_name' in $(round(t_end - t_start, digits=2)) seconds.")
    println()
end

"""
Alternate Green's function calculation with amplitude bounding (testing/experimental).
"""
function calculate_greens_function_bounded(object_name::String; block_name::String = "Greens_block_$(object_name)_L_0_bounded.h5", lambda = 0.0)
    raw_data_path = "data/$object_name/raw_data"
    if !isdir(raw_data_path) || isempty(readdir(raw_data_path))
        throw(ErrorException("Raw data folder is empty or doesn't exist. Please add raw data before calculating Green's function."))
    end
    greens_dir = "data/$object_name/Greens_Blocks"

    if !isdir(greens_dir)
        mkdir(greens_dir)
    end
    greens_block_path = greens_dir * "/" * block_name
    if isfile(greens_block_path)
        println("Warning: File $greens_block_path already exists. Overwriting...")
        rm(greens_block_path)
    end


    sensors = get_sensors(object_name)

    for (i, sensor) in enumerate(sensors)
        raw_data_file = "data/$object_name/raw_data/raw_data_from_sensor_$i.h5"
        if !isfile(raw_data_file)
            continue
        end
        h5file = _open_h5file(raw_data_file)
        sets = keys(h5file)
        finder = GreensFinder(lambda)

        for (s, set_number) in enumerate(sets)
            set_number_int = parse(Int, set_number)
            Add_source_data!(finder, get_single_location_from_raw(object_name, i, sensor, set_number_int), get_raw_data(object_name, i, set_number_int))
        end

        block_response = find_response(finder)
        block_response[:, 1:get_low_freq_idx(object_name)-1] .= 0

        # Downscale all frequencies with abs > 1 to 1
        inds = abs.(block_response) .> 1
        block_response[inds] .= block_response[inds] ./ abs.(block_response[inds])

        h5write(greens_dir*"/$block_name", "sensor_$i", block_response)
    end
    return true
end

"""
Following are a few matrix generators. Currently, they are hard coded in a way to each have their own generator. 
They save the generated matrices under /Inverted_Matrices with a name defined by the user. All of these generators 
really only stack the results of get_system_matrix_1pulse() from MatrixGenerator.jl in different ways. 

General inputs: 
matrix_name: Name of the matrix. This is needed to call it with the evaluater.

Greens_block: If multiple Grens functions have been defined, it may be defined with this. By default, it's named "nothing", which makes the function use the first
Greens function it can find.

lambdas: A list of Ints which defines with which regularization parameters the inversion is performed. 
each inversion may be computationally expensive. Hence it's recomended a simple defect matrix is used initially, to find for which lambda values the inversion is stable.
As a rough guideline, a system with 10'000 points requres roughly a lambda parameter of 0.1-10, a 50x50 grid may only need 0.1-0.001.
The parameter being too low returns nonsense noise data. At a treashold, the data makes sense. Going above that changes it only slightly, smoothing the function out. 


defect_projector: A (n_t, d) dimensional matrix, where n_t is the time axis, and d the defect dimensionality. 
This defines how the defect parameters act on each location. Example:
def_mat = zeros(Float32, n_t, 1)
def_mat[1, 1] = 1
This is a basic 1-dimensional defect matrix which allows each defect to act by refect the incoming wave as it with a given amplitude.
A multidimensional matrix may be:
def_mat = zeros(Float32, n_t, 3)
def_mat[1, 1] = 1
def_mat[2, 2] = 1
def_mat[3, end] = 1
This allows each defect to additionally delay or increment the incoming wave by one time step before re-emitting it. 

Additional parameters are explained in the respective functions
"""

# These functions really should be within MatrixGenerator.jl, for structural clarity.

"""
Creates a system matrix for a single emitter and a set of sensors, then inverts it for multiple lambdas.
"""
function create_single_system_matrix(object_name::String, matrix_name::String, emitter_idx::Int, sensor_idx_array::Array{Int},  
                                    defect_projector; lambdas = [10.0,1.0,0.1], Greens_block::String = "nothing")
    println("Creating single matrix with sensor_array", sensor_idx_array, " and emitter ", emitter_idx)
    
    if is_debugmode(object_name)
        println("Creating single system matrix '$matrix_name' for object '$object_name'")
        println("  emitter_idx: $emitter_idx")
        println("  sensor_idx_array: $sensor_idx_array")
        println("  lambdas: $lambdas")
    end

    t_start = time()

    _find_source_pulses(object_name)

    raw_mat = get_system_matrix_1pulse(object_name, emitter_idx, sensor_idx_array, defect_projector, Greens_block = Greens_block)

    Tikhonov_invert_matrices(object_name, raw_mat,  lambdas = lambdas, file_name = "$matrix_name" * ".h5", debug_mode = is_debugmode(object_name))

    target_dir = "data/$object_name/Inverted_Matrices"
    target_file = "$target_dir/$matrix_name.h5"

    h5open(target_file, "r+") do h5f
        if !("sensor_idx_array" in keys(h5f))
            h5write(target_file, "sensor_idx_array", sensor_idx_array)
        end
        if !("emitter_idx" in keys(h5f))
            h5write(target_file, "emitter_idx", emitter_idx)
        end
        if !("defect_projector" in keys(h5f))
            h5write(target_file, "defect_projector", defect_projector)
        end
        if !("matrix_type" in keys(h5f))
            h5write(target_file, "matrix_type", "single")
        end
    end
    t_end = time()
    println("Finished creating single system matrix '$matrix_name' in $(round(t_end - t_start, digits=2)) seconds.")
    println()
end



"""
Creates a system matrix for all emitters and sensors (full system), then inverts for multiple lambdas.
"""
function create_full_system_matrix(object_name::String, matrix_name::String, defect_projector; 
                                    lambdas = [10, 1, 0.1], Greens_block::String = "nothing")

    println("Creating full matrix with sensor_array")

    if is_debugmode(object_name)
        println("Creating full system matrix '$matrix_name' for object '$object_name'")
        println("  lambdas: $lambdas")
    end
    t_start = time()

    _find_source_pulses(object_name)

    sensors = get_sensors(object_name)

    s = collect(1:length(sensors))

    if ndims(defect_projector) == 1
        defect_projector = reshape(Float32.(defect_projector), :, 1)
    else
        defect_projector = Float32.(defect_projector)
    end

    sensor_array = []
    full_matrix = nothing
    for (e, emitter) in enumerate(sensors)
        raw_mat = get_system_matrix_1pulse(object_name, e, s, defect_projector, Greens_block = Greens_block)
        if full_matrix === nothing
        full_matrix = raw_mat
        else
        full_matrix = vcat(full_matrix, raw_mat)
        end
        sensor_array = append!(sensor_array, sensors)
    end
    
    Tikhonov_invert_matrices(object_name, full_matrix,  lambdas = lambdas, file_name = "$matrix_name" * ".h5", debug_mode = is_debugmode(object_name))

    target_dir = "data/$object_name/Inverted_Matrices"
    target_file = "$target_dir/$matrix_name.h5"
    
    h5open(target_file, "r+") do h5f
        if !("defect_projector" in keys(h5f))
            h5write(target_file, "defect_projector", defect_projector)
        end
        if !("matrix_type" in keys(h5f))
        h5write(target_file, "matrix_type", "full")
         end
    end
    t_end = time()
    println("Finished creating full system matrix '$matrix_name' in $(round(t_end - t_start, digits=2)) seconds.")
    println()
end

"""
Creates a reduced system matrix (no self-measurement), then inverts for multiple lambdas.
"""
function create_reduced_system_matrix(object_name::String, matrix_name::String, defect_projector; 
                                    lambdas = [10, 1, 0.1], Greens_block::String = "nothing")
    # Creates a reduced system matrix where each transducer emits once and all others act as sensors (no self-measurement).
    # Saves the matrix and records emitter and sensor indices for reference.

    println("Creating reduced matrix with sensor_array")

    if is_debugmode(object_name)
        println("Creating reduced system matrix '$matrix_name' for object '$object_name'")
        println("  lambdas: $lambdas")
    end
    t_start = time()
    _find_source_pulses(object_name)

    sensors = get_sensors(object_name)
    s = collect(1:length(sensors))

    if ndims(defect_projector) == 1
        defect_projector = reshape(Float32.(defect_projector), :, 1)
    else
        defect_projector = Float32.(defect_projector)
    end

    emitter_idx_array = Int[]
    sensor_idx_array_per_emitter = Vector{Int}[]

    full_matrix = nothing
    for (e, emitter_idx) in enumerate(sensors)
        s = filter(x -> x != e, s)
        raw_mat = get_system_matrix_1pulse(object_name, e, s, defect_projector, Greens_block = Greens_block)
        if full_matrix === nothing
        full_matrix = raw_mat
        else
        full_matrix = vcat(full_matrix, raw_mat)
        end
        push!(emitter_idx_array, e)
        push!(sensor_idx_array_per_emitter, s)
    end

    
    Tikhonov_invert_matrices(object_name, full_matrix,  lambdas = lambdas, file_name = "$matrix_name" * ".h5", debug_mode = is_debugmode(object_name))

    target_dir = "data/$object_name/Inverted_Matrices"
    target_file = "$target_dir/$matrix_name.h5"
    
    h5open(target_file, "r+") do h5f
        if !("defect_projector" in keys(h5f))
            h5write(target_file, "defect_projector", defect_projector)
        end
        if !("matrix_type" in keys(h5f))
        h5write(target_file, "matrix_type", "reduced")
        end
        if !("emitter_idx_array" in keys(h5f))
            h5write(target_file, "emitter_idx_array", emitter_idx_array)
        end
        if !("sensor_idx_array_per_emitter" in keys(h5f))
            # Save as a 2D array (pad with -1 if needed for rectangular shape)
            maxlen = maximum(length.(sensor_idx_array_per_emitter))
            padded = [vcat(s, fill(-1, maxlen - length(s))) for s in sensor_idx_array_per_emitter]
            h5write(target_file, "sensor_idx_array_per_emitter", hcat(padded...))
        end
    end
    t_end = time()
    println("Finished creating reduced system matrix '$matrix_name' in $(round(t_end - t_start, digits=2)) seconds.")
    println()
end

"""
Evaluates an inverted matrix on a defect set and plots the results. Dispatches to the correct evaluation function based on matrix type.
"""
function evaluate_matrix(object_name::String, inverted_matrix::String, defect_set::String; color_bar = :ice)
    # Each of the matrix generators above also save their type in the h5file. This function reads it, and calls the corresponding evaluation function.
    # The result is saved as h5file, but also Plotted immediately. The Plots differ to the saved result by only being able to be deleted manually. 
    println("Evaluating reduced matrix $inverted_matrix with defect data $defect_set")

    t_start = time()

    matrix_file = "data/$object_name/Inverted_Matrices/$inverted_matrix.h5"
    h5file = h5open(matrix_file, "r")
    matrix_type = read(h5file, "matrix_type")
    result_name = "$(inverted_matrix)"

    if matrix_type == "single"
        evaluate_single_matrix(object_name, result_name, defect_set)
    elseif matrix_type == "full"
        evaluate_full_matrix(object_name, result_name, defect_set)
    elseif matrix_type == "reduced"
        evaluate_reduced_matrix(object_name, result_name, defect_set)
    else
        error("Unknown matrix_type: $matrix_type")
    end

    result_name = "$(inverted_matrix)_$defect_set"

    t_end = time()
    println("Evaluation Completed in $(round(t_end - t_start, digits=2)) seconds.")
    println()

    plot_results(object_name, result_name, color_bar = color_bar)
end
