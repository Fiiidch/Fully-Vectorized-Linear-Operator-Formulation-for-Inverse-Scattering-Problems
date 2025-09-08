using Plots
using HDF5
using FFTW
include("data_reader.jl")
using Measures

# Custom mean function to replace Statistics.mean
function calculate_mean(x)
    return sum(x) / length(x)
end

#map data in the format (n,t) or (n,f) to a 3d matrix of size (x,y,t). Plot along the t axis.
# It implicitly projects the data to the xy plane, so it's only usable for 2d data.

color = :ice

"""
Plots a time-series heatmap animation of data on the xy-plane
"""
function heatmap_3d_plot(object_name::String, data::AbstractArray{<:AbstractFloat, 2}, range; title::String="Heatmap", 
                    xlabel::String="X-axis [m]", ylabel::String="Y-axis [m]", file_name::String="heatmap.gif", fps::Int = 30)
    coords = h5read("data/$(object_name)/object_data.h5", "XYZ")
    _find_source_pulses(object_name)
    if size(data, 1) != size(coords, 1)
        error("Data and coordinates must have the same number of rows.")
    end

    n_t = (length(get_avg_source(object_name, 1, 1)) - 1) * 2
    max_freq = get_max_freq(object_name)
    Δt = 1 /(2 * max_freq)
    # Check if the provided range is within bounds
    if maximum(range) > n_t || minimum(range) < 1
        @warn "Provided range ($(first(range)) to $(last(range))) is out of bounds. Valid range is 1:$n_t. Range is set to 1_$n_t."
        range = max(1, first(range)):min(n_t, last(range))
    end

    time_range = (0:Δt:Δt*(length(range)-1)) .* 1000

    xs = coords[:, 1]
    ys = coords[:, 2]

    # Use tolerance to group coordinates (same as plot_results)
    tol = 0.004
    function group_coords(coords)
        sorted = sort(coords)
        groups = [[sorted[1]]]
        for x in sorted[2:end]
            if isapprox(x, groups[end][1]; atol=tol)
                push!(groups[end], x)
            else
                push!(groups, [x])
            end
        end
        [calculate_mean(g) for g in groups]
    end

    x_unique = group_coords(xs)
    y_unique = group_coords(ys)

    # Map coordinates to grid indices using tolerance
    function find_index(val, unique_vals)
        for (i, u) in enumerate(unique_vals)
            if isapprox(val, u; atol=tol)
                return i
            end
        end
        return nothing
    end

    num_measurements = length(range)
    grid = fill(NaN, length(y_unique), length(x_unique), num_measurements)  # rows = y, cols = x

    for n in 1:size(coords, 1)
        xi = find_index(xs[n], x_unique)
        yi = find_index(ys[n], y_unique)
        if isnothing(xi) || isnothing(yi)
            @warn "Point ($(xs[n]), $(ys[n])) not mapped to grid."
            continue
        end
        grid[yi, xi, :] .= data[n, range]
    end

    global_min = minimum(grid[div(length(y_unique),2), div(length(x_unique),2), :])*2
    global_max = maximum(grid[div(length(y_unique),2), div(length(x_unique),2), :])*2

    anim = @animate for (i, time) in enumerate(time_range)
        t_str = string(round(time; digits=2))
        heatmap(
            x_unique, y_unique, grid[:, :, i],
            aspect_ratio=1,
            c=:balance,
            title="$title at time $t_str ms",
            xlabel=xlabel,
            ylabel=ylabel,
            size=(550, 500),
            clim=(global_min, global_max),
            colorbar_title="\n z-velocity [arbritrary units]",
            right_margin=5mm,
        )
    end

    if !isdir("Plots/$object_name")
        mkpath("Plots/$object_name")
    end

    gif(anim, "Plots/$object_name/$file_name", fps=fps)

end

"""
Same as above, but for complex data (converts to real via irfft).
"""
function heatmap_3d_plot(object_name::String, data::AbstractArray{<:ComplexF32, 2}, time_range; title::String="Heatmap", 
    xlabel::String="X-axis [m]", ylabel::String="Y-axis [m]", file_name::String="heatmap.gif", fps::Int = 30)
    
    target_length = (size(data, 2) - 1) * 2
    data = irfft(data, target_length, 2)
    heatmap_3d_plot(object_name, data, time_range; title=title, xlabel=xlabel, ylabel=ylabel, file_name=file_name, fps=fps)
    return true
end

"""
Plots a heatmap animation of raw data for a given emitter and set.
"""
function heatmap_3d_plot_raw_data(object_name::String, emitter_Index::Int, time_range; set_number = 1)
    data = get_raw_data(object_name, emitter_Index, set_number)
    heatmap_3d_plot(object_name, data, time_range, title = "Raw_data with emitter $emitter_Index", file_name = "Raw_data with emitter $emitter_Index set $set_number.gif")
end
"""
Plots a heatmap animation of raw data with sensor annotations.
"""
function heatmap_3d_plot_raw_data_sensors(object_name::String, emitter_Index::Int, time_range; set_number = 1)
    data = get_raw_data(object_name, emitter_Index, set_number)
    xyz = get_xyz(object_name)
    sensors = get_sensors(object_name)
    
    n_t = (length(get_avg_source(object_name, 1, 1)) - 1) * 2
    max_freq = get_max_freq(object_name)
    Δt = 1 /(2 * max_freq)
    
    # Convert data from frequency domain to time domain
    data = irfft(data, n_t, 2)
    
    # Check if the provided range is within bounds
    if maximum(time_range) > n_t || minimum(time_range) < 1
        @warn "Provided range ($(first(time_range)) to $(last(time_range))) is out of bounds. Valid range is 1:$n_t. Range is set to 1_$n_t."
        time_range = max(1, first(time_range)):min(n_t, last(time_range))
    end

    time_values = (0:Δt:Δt*(length(time_range)-1)) .* 1000

    xs = xyz[:, 1]
    ys = xyz[:, 2]

    # Use tolerance to group coordinates
    tol = 0.004
    function group_coords(coords)
        sorted = sort(coords)
        groups = [[sorted[1]]]
        for x in sorted[2:end]
            if isapprox(x, groups[end][1]; atol=tol)
                push!(groups[end], x)
            else
                push!(groups, [x])
            end
        end
        [calculate_mean(g) for g in groups]
    end

    x_unique = group_coords(xs)
    y_unique = group_coords(ys)

    # Map coordinates to grid indices using tolerance
    function find_index(val, unique_vals)
        for (i, u) in enumerate(unique_vals)
            if isapprox(val, u; atol=tol)
                return i
            end
        end
        return nothing
    end

    num_measurements = length(time_range)
    grid = fill(NaN, length(y_unique), length(x_unique), num_measurements)

    for n in 1:size(xyz, 1)
        xi = find_index(xs[n], x_unique)
        yi = find_index(ys[n], y_unique)
        if isnothing(xi) || isnothing(yi)
            @warn "Point ($(xs[n]), $(ys[n])) not mapped to grid."
            continue
        end
        grid[yi, xi, :] .= data[n, time_range]
    end

    global_min = minimum(grid[div(length(y_unique),2), div(length(x_unique),2), :])*2
    global_max = maximum(grid[div(length(y_unique),2), div(length(x_unique),2), :])*2
    
    # Get sensor annotations
    sensor_annotations = []
    for sensor_idx in sensors
        sensor_x = xyz[sensor_idx, 1]
        sensor_y = xyz[sensor_idx, 2]
        xi = find_index(sensor_x, x_unique)
        yi = find_index(sensor_y, y_unique)
        if !isnothing(xi) && !isnothing(yi)
            push!(sensor_annotations, (x_unique[xi], y_unique[yi]))
        end
    end
    
    anim = @animate for (i, time) in enumerate(time_values)
        t_str = string(round(time; digits=2))
        p = heatmap(
            x_unique, y_unique, grid[:, :, i],
            aspect_ratio=1,
            c=:balance,
            title="Raw data with emitter $emitter_Index at time $t_str ms",
            xlabel="X-axis [m],  ◆ = sensor locations",
            ylabel="Y-axis [m]",
            size=(550, 500),
            clim=(global_min, global_max),
            colorbar_title="\n",
            right_margin=5mm,
        )
        
        # Add sensor annotations
        scatter!(p, [x for (x,y) in sensor_annotations], 
                   [y for (x,y) in sensor_annotations], 
                   color=:red, 
                   markershape=:diamond, 
                   markersize=5, 
                   label="")

        # Add sensor number annotations - positioned slightly higher above the sensor marker
        for (i, (x, y)) in enumerate(sensor_annotations)
            annotate!(p, [(x, y + 0.010, text("$i", :red, 8, :center, :bottom))])
        end
    end

    if !isdir("Plots/$object_name")
        mkpath("Plots/$object_name")
    end

    gif(anim, "Plots/$object_name/Raw_data_with_sensors_emitter_$(emitter_Index)_set_$(set_number).gif", fps=30)
end

"""
Plots the solution(s) from a result file as 2D heatmaps, one per lambda.
"""
function plot_results(object_name::String, result::String; color_bar = color)
    xyz = get_xyz(object_name)
    sensors = get_sensors(object_name)
    result_file = "data/$object_name/results/results_$(result).h5"

    if !isfile(result_file)
        error("Result file $result_file does not exist.")
    end

    h5file = h5open(result_file, "r")
    lambda_keys = filter(key -> startswith(key, "Lambda_"), keys(h5file))
    plot_dir = "Plots/$object_name/results/$result"
    defect_projector = read(h5file, "defect_projector")
    defect_dimensionality = size(defect_projector, 2)
    d = size(defect_projector, 2)
    def_name = read(h5file, "defect_name")
    #   heatmap_from_data(data::Vector{<:Number}, xyz::Vector{SVector{3, <:Number}}, d::Int)

    if !isdir(plot_dir)
        mkpath(plot_dir)
    end

    

    for key in lambda_keys
        data = read(h5file, key)
        n = size(data, 1) ÷ d

        @assert size(xyz, 1) == n "xyz must have n entries"

        vals = [norm(data[(i-1)*d+1:i*d]) for i in 1:n]

        xs = [xyz[i,1] for i in 1:n]
        ys = [xyz[i,2] for i in 1:n]

        # Use tolerance to group coordinates. This is hard coded and not good
        tol = 0.004
        function group_coords(coords)
            sorted = sort(coords)
            groups = [[sorted[1]]]
            for x in sorted[2:end]
                if isapprox(x, groups[end][1]; atol=tol)
                    push!(groups[end], x)
                else
                    push!(groups, [x])
                end
            end
            [calculate_mean(g) for g in groups]
        end

        x_unique = group_coords(xs)
        y_unique = group_coords(ys)

        # Map coordinates to grid indices using tolerance
        function find_index(val, unique_vals)
            for (i, u) in enumerate(unique_vals)
                if isapprox(val, u; atol=tol)
                    return i
                end
            end
            return nothing
        end

        grid = fill(NaN, length(y_unique), length(x_unique))
        mapped = Dict{Tuple{Int,Int},Int}()
        
        for i in 1:n
            xi = find_index(xs[i], x_unique)
            yi = find_index(ys[i], y_unique)
            if isnothing(xi) || isnothing(yi)
                @warn "Point ($(xs[i]), $(ys[i])) not mapped to grid."
                continue
            end
            idx = (yi, xi)
            #=
            if haskey(mapped, idx)
                @warn "Grid cell $idx double-mapped by points $mapped[idx] and $i."
            end
            =#
            grid[yi, xi] = vals[i]
            mapped[idx] = i
        end

        # Mark sensor locations on the grid by setting their values to 0
        sensor_annotations = []
        for (i,sensor_idx) in enumerate(sensors)
            sensor_x = xyz[sensor_idx, 1]
            sensor_y = xyz[sensor_idx, 2]
            
            xi = find_index(sensor_x, x_unique)
            yi = find_index(sensor_y, y_unique)
            
            if !isnothing(xi) && !isnothing(yi)
                grid[yi, xi] = 0.0  # Set sensor locations to 0
            end

            # Also track the position for annotation
            
            push!(sensor_annotations, (x_unique[xi], y_unique[yi]))
            
        end

        key_label = replace(key, "Lambda" => "L")
        p = heatmap(x_unique, y_unique, grid, xlabel="x-axis [m],  ◆ = sensor locations", ylabel="y-axis [m]", c = color_bar, title="$object_name, $key_label, $def_name",
         size=(600, 500), colorbar_title=" \n Defect reflectivity", right_margin=5mm)

        # Add X markers at sensor locations
        scatter!(p, [x for (x,y) in sensor_annotations], 
            [y for (x,y) in sensor_annotations], 
            color=:red, 
            markershape=:diamond, 
            markersize=5, 
            label="")

        savefig(plot_dir * "/solution_$key.png")
    end
    println("Plots saved to $plot_dir")
    println()
end



"""
Plots a greens function kernel for a given sensor and location in the time domain.
"""
function plot_greens_kernel(object_name::String, sensor::Int, location::Int)
    # Get the time domain kernel
    kernel = get_real_greens(object_name, sensor, location)
    n_t = length(kernel)
    max_freq = get_max_freq(object_name)
    Δt = 1 / (2 * max_freq)
    time_range = (0:Δt:Δt*(n_t-1)) .* 1000  # Convert to milliseconds

    # Calculate frequency domain kernel
    freq_kernel = rfft(kernel)
    freq_range = FFTW.rfftfreq(n_t, 1/Δt)./10e2  # Get frequency bins

    # Create plot directory if it doesn't exist
    plot_dir = "Plots/$object_name/greens_kernels"
    if !isdir(plot_dir)
        mkpath(plot_dir)
    end

    # Plot the time domain kernel
    p_time = plot(time_range, kernel, 
        title="Impulse Resp. kernel emitter $sensor to loc. $location", 
        xlabel="Time [ms]", 
        ylabel="Amplitude",
        linewidth=2,
        legend=false)

    # Save the time domain plot
    file_path_time = "$plot_dir/greens_kernel_s$(sensor)_loc$(location)_time.png"
    savefig(p_time, file_path_time)
    
    # Plot the absolute value of frequency domain kernel
    p_freq = plot(freq_range, abs.(freq_kernel), 
        title="Impulse Resp. Freq. Spectrum emitter $sensor to loc. $location", 
        xlabel="Frequency [kHz]", 
        ylabel="Absolute Amplitude",
        linewidth=2,
        legend=false)
    
    # Save the frequency domain plot
    file_path_freq = "$plot_dir/greens_kernel_s$(sensor)_loc$(location)_freq.png"
    savefig(p_freq, file_path_freq)
    
    println("Saved Green's kernel plots to: $file_path_time and $file_path_freq")
end

"""
Plots time-domain defect measurements for all emitter/sensor pairs.
"""
function plot_defects(object_name::String, defect_set::String, range)
    #bug fixing thing. Plots defects, raw sources and what not. requires defect measurements to
    #exist on all emitters.

    target_dir = "Plots/$object_name/defect_visualisation/$defect_set"

    n_t = (length(get_avg_source(object_name, 1, 1)) - 1) * 2
    max_freq = get_max_freq(object_name)
    Δt = 1 /(2 * max_freq)

    # Check if the provided range is within bounds
    if maximum(range) > n_t || minimum(range) < 1
        @warn "Provided range ($(first(range)) to $(last(range))) is out of bounds. Valid range is 1:$n_t. Range is set to 1_$n_t."
        # Ensure range is a valid UnitRange
        range = max(1, first(range)):min(n_t, last(range))
    end

    time_range = (0:Δt:Δt*(length(range)-1)) .* 1000
    sensors = get_sensors(object_name)

    for (e,emit) in enumerate(sensors)
        for (s,sen) in enumerate(sensors)
        def = get_avg_defect_measurement(object_name, e, s, defect_set)
        def = irfft(def, n_t)[range]

        plot(time_range, def, title = "$defect_set from emitter $e, at sensor $s", xlabel = "time [ms]",  
             ylabel = "Normalized velocity", legend = false)
        
        emitter_dir = target_dir * "/emitter_$e"
        if !isdir(emitter_dir)
            mkpath(emitter_dir)
        end
        file_path = target_dir * "/emitter_$e/sensor_$s.png"
        savefig(file_path)
        end
    end
    println("Saved defect plots to directory: $target_dir")
end


"""
Plots time-domain reference (no-defect) measurements for all emitter/sensor pairs.
"""
function plot_reference_measurements(object_name::String, range)
    #bug fixing thing. Plots defects, raw sources and what not. requires defect measurements to
    #exist on all emitters.

    target_dir = "Plots/$object_name/defect_visualisation/reference_measurements"

    n_t = (length(get_avg_source(object_name, 1, 1)) - 1) * 2
    max_freq = get_max_freq(object_name)
    Δt = 1 /(2 * max_freq)

    # Check if the provided range is within bounds
    if maximum(range) > n_t || minimum(range) < 1
        @warn "Provided range ($(first(range)) to $(last(range))) is out of bounds. Valid range is 1:$n_t. Range is set to 1_$n_t."
        # Ensure range is a valid UnitRange
        range = max(1, first(range)):min(n_t, last(range))
    end

    time_range = (0:Δt:Δt*(length(range)-1)) .* 1000
    sensors = get_sensors(object_name)

    for (e,emit) in enumerate(sensors)
        for (s,sen) in enumerate(sensors)
        def = get_avg_source(object_name, e, s)
        def = irfft(def, n_t)[range]

        plot(time_range, def, title = "reference from emitter $e, at sensor $s", 
             xlabel = "time [ms]", ylabel = "Normalized velocity", legend = false)
        
        emitter_dir = target_dir * "/emitter_$e"
        if !isdir(emitter_dir)
            mkpath(emitter_dir)
        end
        file_path = target_dir * "/emitter_$e/sensor_$s.png"
        savefig(file_path)
        end
    end
    println("Saved reference measurement plots to directory: $target_dir")
    
end

"""
Plots the difference between defect and reference measurements for all emitter/sensor pairs.
"""
function plot_source_defect_difference(object_name::String, defect_set::String, range)
    #bug fixing thing. Plots defects, raw sources and what not. requires defect measurements to
    #exist on all emitters.

    target_dir = "Plots/$object_name/defect_visualisation/difference_$defect_set"

    n_t = (length(get_avg_source(object_name, 1, 1)) - 1) * 2
    max_freq = get_max_freq(object_name)
    Δt = 1 /(2 * max_freq)

    # Check if the provided range is within bounds
    if maximum(range) > n_t || minimum(range) < 1
        @warn "Provided range ($(first(range)) to $(last(range))) is out of bounds. Valid range is 1:$n_t. Range is set to 1_$n_t."
        # Ensure range is a valid UnitRange
        range = max(1, first(range)):min(n_t, last(range))
    end

    time_range = (0:Δt:Δt*(length(range)-1)) .* 1000
    sensors = get_sensors(object_name)

    for (e,emit) in enumerate(sensors)
        for (s,sen) in enumerate(sensors)
        def = get_source_defect_diff(object_name, e, s, defect_set)
        def = irfft(def, n_t)[range]

        plot(time_range, def, title = "$defect_set difference from emitter $e, at sensor $s", 
             xlabel = "time [ms]", ylabel = "Normalized velocity", legend = false)
        
        emitter_dir = target_dir * "/emitter_$e"
        if !isdir(emitter_dir)
            mkpath(emitter_dir)
        end
        file_path = target_dir * "/emitter_$e/sensor_$s.png"
        savefig(file_path)
        end
    end
    println("Saved source-defect difference plots to directory: $target_dir")
    println("Files saved:")
    
end

"""
Plots the average Green's function magnitude for points at a certain distance from a reference sensor.
"""
function plot_greens_function_analysis(object_name::String; reference_sensor::Int=5, distance_target=10.0)

    _find_source_pulses(object_name)

    xyz = get_xyz(object_name)
    sensors = get_sensors(object_name)
    sensor_idx = reference_sensor
    sensor_coord = xyz[sensor_idx, :]
    n_t = get_n_t(object_name)
    data_freq = get_data_freq(object_name)
    freqs = FFTW.rfftfreq(n_t, data_freq)

    # Print rough boundaries of xyz coordinates
    min_xyz = minimum(xyz, dims=1)
    max_xyz = maximum(xyz, dims=1)
    println("XYZ boundaries: x ∈ [$(min_xyz[1]), $(max_xyz[1])], y ∈ [$(min_xyz[2]), $(max_xyz[2])]")

    if size(xyz, 2) == 3
        println("z ∈ [$(min_xyz[3]), $(max_xyz[3])]")
    end

    # Calculate distances
    dist = [norm(xyz[i, :] .- sensor_coord) for i in 1:size(xyz,1)]
    # Select 10 random points near distance_target ± 2
    deviation = 0.2 * distance_target  # 20% of the target distance
    near_points = findall(x -> abs(x - distance_target) < deviation, dist)
    sample_points = rand(near_points, min(10, length(near_points)))

    println("Sample points found: $(length(sample_points)) / 10")

    # Load Green's block for reference sensor
    gblock = get_greens_block(object_name, reference_sensor)
    max_freq = get_max_freq(object_name)

    # Derive frequency axis as a vector to avoid callable StepRange issues
    n_freq_bins = size(gblock, 2)
    freqs = LinRange(0, max_freq, n_freq_bins)
    avg_mag = sum(abs.(gblock[sample_points, :]), dims=1) ./ length(sample_points)

    target_dir = "Plots/$object_name/Impulse_response/"
    if !isdir(target_dir)
        mkpath(target_dir)
    end
    # Plot result

    println(size(freqs), size(avg_mag))
    plot(freqs, vec(avg_mag), xlabel="Frequency (Hz)", ylabel="Amplitude",
         title="Green's Function Analysis sensor $reference_sensor")
    file_path = target_dir * "greens_analysis.png"
    savefig(file_path)
    println("Saved Green's function analysis plot to: $file_path")
end

"""
Plots a heatmap animation of a Green's function block for a sensor.
"""
function plot_greens(object_name::String, sensor::Int, range)
    block = get_greens_block(object_name, sensor)
    heatmap_3d_plot(object_name, block, range, file_name = "greens_$sensor.gif")
end

"""
Plots a heatmap animation of Green's function block multiplied by the source pulse.
"""
function plot_greens_times_source(object_name::String, sensor::Int, range)
    block = get_greens_block(object_name, sensor)
    block =_get_field_times_source(block, get_source(object_name, sensor))
    heatmap_3d_plot(object_name, block, range, file_name = "greens_$(sensor)_quality_control.gif")
end


"""
Plots the time-domain signal for a single sensor/location/set.
"""
function plot_single_location(object_name::String, sensor::Int, location::Int, set_index::Int)

    vec = get_single_location_from_raw(object_name, sensor, location, set_index)
    
    n_t = (length(get_avg_source(object_name, 1, 1)) - 1) * 2
    max_freq = get_max_freq(object_name)
    Δt = 1 /(2 * max_freq)
    time_range = (0:Δt:Δt*(n_t-1)) .* 1000

    # Convert to time domain if needed
    if length(vec) == div(n_t, 2) + 1
        # Frequency domain, convert to time domain
        data = irfft(vec, n_t)
    else
        data = vec
    end

    plt = plot(time_range, real(data), xlabel="Time (ms)", ylabel="Amplitude", title="Single Location Signal")
    plot_dir = "Plots/$object_name"
    if !isdir(plot_dir)
        mkpath(plot_dir)
    end
    file_path = "$plot_dir/single_location_s_$(sensor)_loc_$location.png"
    savefig(file_path)
    println("Saved single location plot to: $file_path")
end

"""
Plots a heatmap animation of raw data, highlighting a specific location.
"""
function heatmap_3d_plot_raw_data_and_loc(object_name::String, emitter_Index::Int, time_range, location::Int; set_number = 1)
    data = get_raw_data(object_name, emitter_Index, set_number)
    data[location,1] = 100
    heatmap_3d_plot(object_name, data, time_range, title = "Raw_data with emitter $emitter_Index", file_name = "Raw_data_and loc with emitter $emitter_Index set $set_number.gif")
end
"""
    Calls three plotting functions, since usually, these three are called together anyways.
"""
function analyze_source_terms(object_name::String, defect_set::String, range)
    plot_source_defect_difference(object_name, defect_set, range)
    plot_defects(object_name, defect_set, range)
    plot_reference_measurements(object_name, range)
end

function plot_defect_sideways(object_name, matrix, defect_name, axis, relative_max; num_ticks = 6)
    @assert axis == "x" || axis == "y" "axis must be either \"x\" or \"y\""
    
    plot_dir = "Plots/$object_name/sideview_plots_$(axis)axis/defect_sideways_$(axis)-axis_$(defect_name)_$(matrix)"
    if !isdir(plot_dir)
        mkpath(plot_dir)
    end

    xyz = get_xyz(object_name)
    sensors = get_sensors(object_name)
    result_file = "data/$object_name/results/results_$(matrix)_$(defect_name).h5"

    if !isfile(result_file)
        error("Result file $result_file does not exist.")
    end

    h5file = h5open(result_file, "r")
    lambda_keys = filter(key -> startswith(key, "Lambda_"), keys(h5file))
    defect_projector = read(h5file, "defect_projector")
    d = size(defect_projector, 2)
    def_name = read(h5file, "defect_name")
    #   heatmap_from_data(data::Vector{<:Number}, xyz::Vector{SVector{3, <:Number}}, d::Int)

    if !isdir(plot_dir)
        mkpath(plot_dir)
    end

    for key in lambda_keys
        data = read(h5file, key)
        n = size(data, 1) ÷ d

        @assert size(xyz, 1) == n "xyz must have n entries"

        vals = [norm(data[(i-1)*d+1:i*d]) for i in 1:n]

        xs = [xyz[i,1] for i in 1:n]
        ys = [xyz[i,2] for i in 1:n]

        # Use tolerance to group coordinates. This is hard coded and not good
        tol = 0.004
        function group_coords(coords)
            sorted = sort(coords)
            groups = [[sorted[1]]]
            for x in sorted[2:end]
                if isapprox(x, groups[end][1]; atol=tol)
                    push!(groups[end], x)
                else
                    push!(groups, [x])
                end
            end
            [calculate_mean(g) for g in groups]
        end

        x_unique = group_coords(xs)
        y_unique = group_coords(ys)

        # Map coordinates to grid indices using tolerance
        function find_index(val, unique_vals)
            for (i, u) in enumerate(unique_vals)
                if isapprox(val, u; atol=tol)
                    return i
                end
            end
            return nothing
        end

        grid = fill(NaN, length(y_unique), length(x_unique))
        mapped = Dict{Tuple{Int,Int},Int}()
        
        for i in 1:n
            xi = find_index(xs[i], x_unique)
            yi = find_index(ys[i], y_unique)
            if isnothing(xi) || isnothing(yi)
                @warn "Point ($(xs[i]), $(ys[i])) not mapped to grid."
                continue
            end
            idx = (yi, xi)
            #=
            if haskey(mapped, idx)
                @warn "Grid cell $idx double-mapped by points $mapped[idx] and $i."
            end
            =#
            grid[yi, xi] = vals[i]
            mapped[idx] = i
        end

        # Mark sensor locations on the grid by setting their values to 0
        for sensor_idx in sensors
            sensor_x = xyz[sensor_idx, 1]
            sensor_y = xyz[sensor_idx, 2]
            
            xi = find_index(sensor_x, x_unique)
            yi = find_index(sensor_y, y_unique)
            
            if !isnothing(xi) && !isnothing(yi)
                grid[yi, xi] = 0.0  # Set sensor locations to 0
            end
        end

        key_label = replace(key, "Lambda" => "L")
        title="$object_name, $key_label, $def_name"
        x_axis = "axis $(axis) [m]"
        
        if axis == "y"
            slices = 1:length(y_unique)
        else
            slices = 1:length(x_unique)
        end

        max_amp_list = []
        global_max = maximum(grid)
        threshold = global_max * relative_max
        threshold_slices = []
        if axis == "y"
            for y_slice in slices
                append!(max_amp_list, maximum(grid[:, y_slice]))
            end
        else
            for x_slice in slices
               append!(max_amp_list, maximum(grid[x_slice, :]))
            end
        end
        
        p = plot()
        for (idx, slice_num) in enumerate(slices)
            # Skip slices below threshold
            slice_max = (axis == "y") ? maximum(grid[:, slice_num]) : maximum(grid[slice_num, :])
            if slice_max < threshold
                continue
            end
            append!(threshold_slices, slice_num)
            # Plot slice
            if axis == "y"
                slice_vals = grid[:, slice_num]
                plot!(x_unique, slice_vals, label="y=$(round(y_unique[slice_num]; digits=3))m")
            else
                slice_vals = grid[slice_num, :]
                plot!(y_unique, slice_vals, label="x=$(round(x_unique[slice_num]; digits=3))m")
            end
        end

        #make an outline plot if multiple lines exist
        
        if axis == "x"
            grid_threshold = grid[threshold_slices, :]
            outline = [maximum(grid_threshold[:, y]) for y in slices]
        else
            grid_threshold = grid[:, threshold_slices]
            outline = [maximum(grid_threshold[x, :]) for x in slices]
        end
        if length(threshold_slices) > 1
            if axis == "x"
                plot!(x_unique, outline, color="black", label="Max Outline")
            else
                plot!(y_unique, outline, color="black", label="Max Outline")
            end
        end


        # Adjust y-axis limit to be 20% higher to accommodate annotations
        y_max = maximum(outline) * 1.2  # 20% higher than the max value
        
        hline!([threshold], label="Threshold $(round(threshold; digits=3))", linestyle=:dash, title = "Side_view $object_name, $defect_name", 
            xaxis = "$axis-axis [m]", yaxis = "defect reflectivity", ylims=(0, y_max))
        
        
        # Find where the outline crosses the threshold
        
        crossing_points = []
        axis_vals = axis == "x" ? y_unique : x_unique
        
        # Find all crossing points
        for i in 1:(length(outline)-1)
            if (outline[i] <= threshold && outline[i+1] >= threshold) || 
               (outline[i] >= threshold && outline[i+1] <= threshold)
                # Linear interpolation to find exact crossing
                t = (threshold - outline[i]) / (outline[i+1] - outline[i])
                cross_point = axis_vals[i] + t * (axis_vals[i+1] - axis_vals[i])
                push!(crossing_points, (cross_point, i))
            end
        end
        
        # Get first and last crossing points if they exist
        if !isempty(crossing_points)
            first_cross = first(crossing_points)[1]
            last_cross = last(crossing_points)[1]
            
            # Add vertical lines at crossing points
            vline!([first_cross, last_cross], label = "Threshold crossing", linestyle=:dot)
            
            # Add width annotation directly over the crossing points
            width = abs(last_cross - first_cross)
            midpoint = (first_cross + last_cross) / 2
            peak_height = maximum(outline) * 1.1  # Position it at 110% of the peak height
            annotate!([(midpoint, peak_height, 
            text("Width: $(round(width; digits=3))m", 10))])
        end
        
        # Use more standard ticks with equal steps (6 ticks)
        
        tick_min = minimum(axis_vals)
        tick_max = maximum(axis_vals)
        tick_step = (tick_max - tick_min) / (num_ticks - 1)
        tick_positions = [tick_min + i * tick_step for i in 0:(num_ticks-1)]
        tick_labels = ["$(round(pos; digits=3))" for pos in tick_positions]
        
        plot!(xticks=(tick_positions, tick_labels))
        
        savefig("$plot_dir/section_scan_$(key_label)_$(axis).png")
    end

end

"""
Plots the maximum outline of defect data across all slices along the specified axis.
"""
function plot_defect_full_sideways(object_name, matrix, defect_name, axis)
    @assert axis == "x" || axis == "y" "axis must be either \"x\" or \"y\""
    
    plot_dir = "Plots/$object_name/sideview_plots_full_$(axis)-axis/defect_sideways_full_$(axis)_axis_$(defect_name)_$(matrix)"
    if !isdir(plot_dir)
        mkpath(plot_dir)
    end

    xyz = get_xyz(object_name)
    result_file = "data/$object_name/results/results_$(matrix)_$(defect_name).h5"

    if !isfile(result_file)
        error("Result file $result_file does not exist.")
    end

    h5file = h5open(result_file, "r")
    lambda_keys = filter(key -> startswith(key, "Lambda_"), keys(h5file))
    defect_projector = read(h5file, "defect_projector")
    d = size(defect_projector, 2)
    def_name = read(h5file, "defect_name")

    for key in lambda_keys
        data = read(h5file, key)
        n = size(data, 1) ÷ d

        @assert size(xyz, 1) == n "xyz must have n entries"

        vals = [norm(data[(i-1)*d+1:i*d]) for i in 1:n]

        xs = [xyz[i,1] for i in 1:n]
        ys = [xyz[i,2] for i in 1:n]

        # Use tolerance to group coordinates
        tol = 0.004
        function group_coords(coords)
            sorted = sort(coords)
            groups = [[sorted[1]]]
            for x in sorted[2:end]
                if isapprox(x, groups[end][1]; atol=tol)
                    push!(groups[end], x)
                else
                    push!(groups, [x])
                end
            end
            [calculate_mean(g) for g in groups]
        end

        x_unique = group_coords(xs)
        y_unique = group_coords(ys)

        # Map coordinates to grid indices using tolerance
        function find_index(val, unique_vals)
            for (i, u) in enumerate(unique_vals)
                if isapprox(val, u; atol=tol)
                    return i
                end
            end
            return nothing
        end

        grid = fill(NaN, length(y_unique), length(x_unique))
        
        for i in 1:n
            xi = find_index(xs[i], x_unique)
            yi = find_index(ys[i], y_unique)
            if isnothing(xi) || isnothing(yi)
                @warn "Point ($(xs[i]), $(ys[i])) not mapped to grid."
                continue
            end
            grid[yi, xi] = vals[i]
        end

        key_label = replace(key, "Lambda" => "L")
        title = "$object_name, $key_label, $def_name"
        
        # Create outline of maximum values
        if axis == "x"
            outline = [maximum(grid[:, y]) for y in 1:length(y_unique)]
            axis_vals = y_unique
            x_label = "X-axis [m]"
        else
            outline = [maximum(grid[x, :]) for x in 1:length(x_unique)]
            axis_vals = x_unique
            x_label = "Y-axis [m]"
        end

        # Calculate FWHM (Full Width at Half Maximum)
        max_value = maximum(outline)
        half_max = max_value / 2
        
        # Find all crossing points at half maximum
        crossing_points = []
        for i in 1:(length(outline)-1)
            if (outline[i] <= half_max && outline[i+1] >= half_max) || 
                (outline[i] >= half_max && outline[i+1] <= half_max)
                # Linear interpolation to find exact crossing
                t = (half_max - outline[i]) / (outline[i+1] - outline[i])
                cross_point = axis_vals[i] + t * (axis_vals[i+1] - axis_vals[i])
                push!(crossing_points, (cross_point, i))
            end
        end
        
        # Plot the outline
        p = plot(axis_vals, outline, 
            label="Max outline", 
            title="Side view $object_name, $defect_name", 
            xlabel=x_label, 
            ylabel="Defect reflectivity",
            lw=2,
            color="black")
            
        
        
        savefig("$plot_dir/full_outline_$(key_label)_$(axis).png")
    end
end


