using HDF5
using FFTW

"""
Functions for reading and converting data from data/"object_name" folders.
Underscore prefix: internal/helper/debug functions.
"""

"""
Extracts and saves source pulses and reference (no-defect) measurements from raw data.
Skips if already done.
"""
function _find_source_pulses(object_name::String)
    

    target_dir = "data/$(object_name)/defect_measurements/source_pulses"

    if isdir(target_dir)
        return 
    end

    if is_debugmode(object_name)
        println("_find_source_pulses() called")
    end

    mkdir(target_dir)

    ref_dir = "data/$(object_name)/defect_measurements/reference_defect_measurements"
    
    if !isdir(ref_dir)
        mkdir(ref_dir)
    end

    sensors = get_sensors(object_name)

    # Calculate average sources from raw data at each sensor. 
    for (s, sensor) in enumerate(sensors)  
        pulse = zeros(ComplexF32, get_freq_range(object_name))
        raw_data_file = "data/$object_name/raw_data/raw_data_from_sensor_$s.h5"

        if !isfile(raw_data_file)
            println("Warning: No raw data found for sensor $s, (location $sensor)")
            continue
        end

        for i in 1:_max_set_number_of_raw_data(object_name, s)
            pulse .+= get_single_location_from_raw(object_name, s, sensor, i)
        end

        pulse ./= _max_set_number_of_raw_data(object_name, s)
        # Transform to time domain
        #pulse_time = irfft(pulse, length(pulse) * 2 - 1)
        # Keep only the first 1/16th entries
        #n_keep = length(pulse_time) ÷ 16
        #pulse_time[n_keep+1:end] .= 0.0f0
        # Transform back to frequency domain
        #pulse = rfft(pulse_time)

        h5write("$target_dir/source_pulse_from_sensor_$s.h5", "pulse", pulse)
    end

    n_sensors = length(sensors)

    # save reference data, as defined above in the same format as the defect measurments have been saved.
    for emitter_idx in 1:n_sensors
        emitter_sensor = sensors[emitter_idx]
        n_sets = _max_set_number_of_raw_data(object_name, emitter_idx)
        ref_file = "$ref_dir/measurements_from_emitter_$(emitter_idx).h5"
        
        if !isfile(ref_file)
            h5write(ref_file, "init", 0)
        end

        h5f = h5open(ref_file, "r+")

        
        for receiver_idx in 1:n_sensors
            receiver_sensor = sensors[receiver_idx]
            group_name = "sensor_$(receiver_idx)"
            if !(group_name in keys(h5f))
                create_group(h5f, group_name)
            end
            group_keys = keys(h5f[group_name])
            # Add all sets for this emitter/receiver pair
            for set_number in 1:n_sets
                # Only add if not already present (by set index)
                set_index = length(keys(h5f[group_name])) + 1
                measurement = get_single_location_from_raw(object_name, emitter_idx, receiver_sensor, set_number)
                h5write(ref_file, "$group_name/$set_index", measurement)
            end
        end

        close(h5f)
        
    end
end

"""
Open HDF5 file with error message if missing.
"""
function _open_h5file(directory::String)
    if !isfile(directory)
        throw(ArgumentError("File $directory does not exist."))
    end
    return h5open(directory, "r")
end
"""
Read a key from an open HDF5 file, close after. Better error on missing key.
"""
function _read_h5file(h5file, key)
    # Reads a key from an open HDF5 file and closes it.
    all_keys = keys(h5file)
    if !("$key" in all_keys)
        throw(ArgumentError("Key $key not found in h5file. Available keys are: $all_keys"))

    end
    data = read(h5file, "$key")
    close(h5file)
    return data
end 
"""
# Returns true if debug mode is enabled for this object.
"""
function is_debugmode(object_name)
    
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    mode = _read_h5file(h5file, "debug_mode")
    return mode
end
"""
# Returns list of sensor indices for the object.
"""
function get_sensors(object_name::String)::Vector{Int}
    
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    sensors = _read_h5file(h5file, "sensors")
    return sensors
end

"""
Returns 2D array of [x, y, z] coordinates for each point/sensor.
Each row is a location.
"""
function get_xyz(object_name::String)::Array{Float32}
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    data = _read_h5file(h5file, "XYZ")
    return data
end

"""
Returns the low end of the bandpass filter (Float32).
"""
function get_low_freq(object_name::String)::Float32
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    low_freq = _read_h5file(h5file, "low_freq")
    return low_freq
end
"""
Returns the high end of the bandpass filter (Float32).
"""
function get_max_freq(object_name::String)::Float32
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    max_freq = _read_h5file(h5file, "max_freq")
    return max_freq
end
"""
Returns n_t, the length of the input data in time domain.
"""
function get_n_t(object_name::String)::Int
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    n_t = _read_h5file(h5file, "n_t")
    return n_t
end
"""
Returns data_freq, the sampling frequency of the input data.
"""
function get_data_freq(object_name::String)::Float32
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    data_freq = _read_h5file(h5file, "data_freq")
    return data_freq
end
"""
Returns the index where the bandpass filter starts (low freq cutoff).
"""
function get_low_freq_idx(object_name::String)::Int
    low_freq = get_low_freq(object_name)
    data_freq = get_data_freq(object_name)
    n_t = get_n_t(object_name)
    freqs = FFTW.rfftfreq(n_t, data_freq)
    low_freq_idx = findfirst(freqs .>= low_freq)
    return low_freq_idx
end
"""
Returns the length of the frequency domain data (bandpassed).
"""
function get_freq_range(object_name::String)::Int
    directory = "data/$object_name/object_data.h5"
    h5file = _open_h5file(directory)
    f = _read_h5file(h5file, "freq_range")
    return f
end

"""
Returns raw data as a 2D complex array (location, frequency) for a given emitter and set.
"""
function get_raw_data(object_name::String, sensor_Index::Int, set_number::Int)::Array{ComplexF32, 2}
    directory = "data/$object_name/raw_data/raw_data_from_sensor_$sensor_Index.h5"
    h5file = _open_h5file(directory)
    data = _read_h5file(h5file, set_number)
    return data
end
"""
Loads a specific location from raw data for all frequencies.
"""
function get_single_location_from_raw(object_name::String, origin_sensor_index::Int, location::Int, set_number::Int)::Array{ComplexF32}
    data_subset = get_raw_data(object_name, origin_sensor_index, set_number)[location,:]
    return data_subset
end
"""
Returns how many sets of raw data exist for a given sensor.
"""
function _max_set_number_of_raw_data(object_name::String, sensor_index::Int)::Int
    directory = "data/$object_name/raw_data/raw_data_from_sensor_$sensor_index.h5"
    h5file = _open_h5file(directory)
    sets = keys(h5file)
    close(h5file)
    return length(sets)
end
"""
Returns the Green's function block as a 2D complex array (location, frequency).
Greens_block can be specified or defaults to the first found.
"""
function get_greens_block(object_name::String, sensor::Int; Greens_block::String = "nothing")::Array{ComplexF32,2}
    greens_block_directory = "data/$object_name/Greens_Blocks"
    
    if Greens_block == "nothing"
        # Get first .h5 file in directory
        h5_files = filter(x -> endswith(x, ".h5"), readdir(greens_block_directory))
        if isempty(h5_files)
            throw(ArgumentError("No .h5 files found in $greens_block_directory"))
        end
        Greens_block = first(h5_files)
    else
        Greens_block = "$Greens_block"
        if !endswith(Greens_block, ".h5")
            Greens_block *= ".h5"
        end
    end

    directory = joinpath(greens_block_directory, Greens_block)
    h5file = _open_h5file(directory)
    data = _read_h5file(h5file, "sensor_$sensor")
    return data
end

function get_real_greens(object_name::String, sensor::Int, location::Int; Greens_block::String = "nothing")
    greens_block_directory = "data/$object_name/Greens_Blocks"
    
    if Greens_block == "nothing"
        # Get first .h5 file in directory
        h5_files = filter(x -> endswith(x, ".h5"), readdir(greens_block_directory))
        if isempty(h5_files)
            throw(ArgumentError("No .h5 files found in $greens_block_directory"))
        end
        Greens_block = first(h5_files)
    else
        Greens_block = "$Greens_block"
        if !endswith(Greens_block, ".h5")
            Greens_block *= ".h5"
        end
    end

    directory = joinpath(greens_block_directory, Greens_block)
    h5file = _open_h5file(directory)
    data = _read_h5file(h5file, "sensor_$sensor")[location,:]
    data = irfft(data, length(data) * 2 - 1)
    return data
end



"""
Returns the average source pulse for a sensor.
"""
function get_source(object_name::String, sensor_index::Int)::Array{ComplexF32}

    _find_source_pulses(object_name) #Calculates the needed data if it hasn't done so before
    
    directory = "data/$object_name/defect_measurements/source_pulses/source_pulse_from_sensor_$sensor_index.h5"
    h5file = _open_h5file(directory)
    data = _read_h5file(h5file, "pulse")
    return data
end

"""
Multiplies each row of a Green's block by a source (elementwise).
Mostly for debugging. (Plot to compare Greens's function vs raw data)
"""
function _get_field_times_source(greens_block::AbstractArray{<:ComplexF32, 2}, source::Array{ComplexF32})::Array{ComplexF32,2}
    #really only for debugging and quality control. Takes a greens block, and convolves every time row with the source.
    data = zeros(ComplexF32, size(greens_block, 1), size(greens_block, 2))
    for i in 1:size(greens_block, 1)
        data[i, :] = (greens_block[i, :]) .* (source)
    end
    return data
end

"""
Returns the average reference (no-defect) measurement for emitter/sensor pair.
"""
function get_avg_source(object_name::String, emitter_idx::Int, sensor_idx::Int)::Array{ComplexF32,1}
    # Fetch and average all reference (no-defect) measurements for emitter_idx and sensor_idx
    ref_file = "data/$object_name/defect_measurements/reference_defect_measurements/measurements_from_emitter_$(emitter_idx).h5"
    group_name = "sensor_$(sensor_idx)"
    if !isfile(ref_file)
        throw(ArgumentError("Reference file $ref_file does not exist."))
    end
    h5f = h5open(ref_file, "r")
    
    if !(group_name in keys(h5f))
        throw(ArgumentError("Group $group_name not found in $ref_file"))
    end
    group = h5f[group_name]
    keys_list = sort(collect(keys(group)), by=x->parse(Int, x))
    n = length(keys_list)
    if n == 0
        throw(ArgumentError("No measurements found in $group_name of $ref_file"))
    end
    sum_vec = nothing
    for k in keys_list
        vec = read(group[k])
        if sum_vec === nothing
            sum_vec = zeros(eltype(vec), size(vec))
        end
        sum_vec .+= vec
    end
    avg_vec = sum_vec ./ n
    return avg_vec

    close(h5f)
    
end

"""
Returns the average defect measurement for emitter/sensor pair and defect set.
"""
function get_avg_defect_measurement(object_name::String, emitter_idx::Int, sensor_idx::Int, defect_set::String)::Array{ComplexF32}
    # Fetch and average all defect measurements for emitter_idx and sensor_idx in the given defect_set
    file = "data/$object_name/defect_measurements/$defect_set/measurements_from_emitter_$(emitter_idx).h5"
    group_name = "sensor_$(sensor_idx)"
    if !isfile(file)
        throw(ArgumentError("Defect measurement file $file does not exist."))
    end
    h5f = h5open(file, "r")
    
    if !(group_name in keys(h5f))
        throw(ArgumentError("Group $group_name not found in $file"))
    end
    group = h5f[group_name]
    keys_list = sort(collect(keys(group)), by=x->parse(Int, x))
    n = length(keys_list)
    if n == 0
        throw(ArgumentError("No measurements found in $group_name of $file"))
    end
    sum_vec = nothing
    for k in keys_list
        vec = read(group[k])
        if sum_vec === nothing
            sum_vec = zeros(eltype(vec), size(vec))
        end
        sum_vec .+= vec
    end
    avg_vec = sum_vec ./ n
    return avg_vec

    close(h5f)
    
end

"""
Returns the difference between average defect and reference measurement.
"""
function get_source_defect_diff(object_name::String, emitter_idx::Int, sensor_idx::Int, defect_set::String)::Array{ComplexF32}
    avg_source = get_avg_source(object_name, emitter_idx, sensor_idx)
    avg_defect = get_avg_defect_measurement(object_name, emitter_idx, sensor_idx, defect_set)
    return avg_defect .- avg_source
end

"""
Returns a single reference measurement (no averaging). For debugging.
"""
function _get_single_source(object_name::String, emitter_idx::Int, sensor_idx::Int, set_index::Int)
    # Fetch a single reference (no-defect) measurement for emitter_idx and sensor_idx at set_index, for debugging
    ref_file = "data/$object_name/defect_measurements/reference_defect_measurements/measurements_from_emitter_$(emitter_idx).h5"
    group_name = "sensor_$(sensor_idx)"
    if !isfile(ref_file)
        throw(ArgumentError("Reference file $ref_file does not exist."))
    end
    h5f = h5open(ref_file, "r")
    try
        if !(group_name in keys(h5f))
            throw(ArgumentError("Group $group_name not found in $ref_file"))
        end
        group = h5f[group_name]
        keys_list = sort(collect(keys(group)), by=x->parse(Int, x))
        if set_index < 1 || set_index > length(keys_list)
            throw(ArgumentError("set_index $set_index out of bounds for $group_name in $ref_file"))
        end
        vec = read(group[keys_list[set_index]])
        return vec
    finally
        close(h5f)
    end
end

"""
Returns a single defect measurement (no averaging). For debugging.
"""
function _get_single_defect_measurement(object_name::String, emitter_idx::Int, sensor_idx::Int, defect_set::String, set_index::Int)
    # Fetch a single defect measurement for emitter_idx and sensor_idx in the given defect_set at set_index, for debugging
    file = "data/$object_name/defect_measurements/$defect_set/measurements_from_emitter_$(emitter_idx).h5"
    group_name = "sensor_$(sensor_idx)"
    if !isfile(file)
        throw(ArgumentError("Defect measurement file $file does not exist."))
    end
    h5f = h5open(file, "r")
    try
        if !(group_name in keys(h5f))
            throw(ArgumentError("Group $group_name not found in $file"))
        end
        group = h5f[group_name]
        keys_list = sort(collect(keys(group)), by=x->parse(Int, x))
        if set_index < 1 || set_index > length(keys_list)
            throw(ArgumentError("set_index $set_index out of bounds for $group_name in $file"))
        end
        vec = read(group[keys_list[set_index]])
        return vec
    finally
        close(h5f)
    end
end