include("Modules/command_center.jl")
include("Modules/custom_plotters.jl")
using HDF5
using FFTW

object_name = "Showcase"
original_n_t = 2000 # Number of time steps of data before any changes.
data_freq = 200_000   # Frequency of the data in Hz
num_noise = 1 #number of measurements per data set. used for noise compensation.
noise_amplitude = 0. # Amplitude of noise added to the data, as a percentage of the RMS signal level.

n_t = 2000 # Number of target steps of data.
min_freq = 0.0  #min frequency of bandpass filter
max_freq = 100_000.0 #max frequency of bandpass filter (data_freq/2 does not filter at all)

range = 1:100   #range of data for time axis plotting. Amount in indexes, not time steps.

debug_mode = false #if true spams print statements for every step

# Ensure n_t doesn't exceed original_n_t
n_t = min(n_t, original_n_t)
data_range = 1:n_t

data_dir = "Simulated_data_downsampled/"

if isdir("Plots/$object_name")
    rm("Plots/$object_name", recursive=true)
end

raw_data_dir = "$(data_dir)fz_homo_q/SurfaceVelocity/"
defect_A_dir = "$(data_dir)fz_defect_A/"
defect_B_dir = "$(data_dir)fz_defect_B/"
defect_C_dir = "$(data_dir)fz_defect_C/"
defect_AC_dir = "$(data_dir)fz_defect_AC/"
defect_AB_dir = "$(data_dir)fz_defect_AB/"
defect_A5_dir = "$(data_dir)fz_defect_A5/"
defect_A5B_dir = "$(data_dir)fz_defect_A5B/"
defect_D_dir = "$(data_dir)fz_defect_D/"
defect_E_dir = "$(data_dir)fz_defect_E/"
defect_AC_square_dir = "$(data_dir)fz_defect_AC_square/"

defect_A_filename_prefix = "fz_defect_A_SrcEvent"
defect_B_filename_prefix = "fz_defect_B_SrcEvent"
defect_C_filename_prefix = "fz_defect_C_SrcEvent"
defect_AC_filename_prefix = "fz_defect_AC_SrcEvent"
defect_AB_filename_prefix = "fz_defect_AB_SrcEvent"
defect_A5_filename_prefix = "fz_defect_A5_SrcEvent"
defect_A5B_filename_prefix = "fz_defect_A5B_SrcEvent"
defect_D_filename_prefix = "fz_defect_D_SrcEvent"
defect_E_filename_prefix = "fz_defect_E_SrcEvent"
defect_AC_square_prefix = "fz_defect_AC_square_SrcEvent"

# List of defect directories and their corresponding filename prefixes
defect_dirs = [
    defect_A_dir,
    defect_B_dir,
    defect_C_dir,
    defect_AC_dir,
    defect_AB_dir,
    defect_A5_dir,
    defect_A5B_dir,
    defect_D_dir,
    defect_E_dir,
    defect_AC_square_dir
]

defect_prefixes = [
    defect_A_filename_prefix,
    defect_B_filename_prefix,
    defect_C_filename_prefix,
    defect_AC_filename_prefix,
    defect_AB_filename_prefix,
    defect_A5_filename_prefix,
    defect_A5B_filename_prefix,
    defect_D_filename_prefix,
    defect_E_filename_prefix,
    defect_AC_square_prefix
]
# List of defect names
defect_names = ["defect_A", "defect_B", "defect_C", "defect_AC", "defect_AB", "defect_A5", "defect_A5B",
    "defect_D", "defect_E", "defect_AC_square"
]
# Following section remaps the x-y system of the data to a n- system.
# It is specific to the data set and should be adapted for other data sets.
x_max = 101
y_max = 101

h5file = h5open(defect_A_dir * "fz_defect_A_SrcEvent01.h5", "r")
x_arr = read(h5file, "x")
y_arr = read(h5file, "y")

src01 = [x_arr[1], y_arr[1]]
src02 = [x_arr[2], y_arr[2]]
src03 = [x_arr[3], y_arr[3]]
src04 = [x_arr[4], y_arr[4]]
src05 = [x_arr[5], y_arr[5]]

sources = [src01 , src02 , src03 , src04 , src05]

xyz = Array{Float32}(undef, y_max * x_max, 3)
idx = 1

h5file = h5open(raw_data_dir * "fz_homo_q_fz_src01.h5", "r")
x_arr = read(h5file, "x")
y_arr = read(h5file, "y")

for col in 1:x_max
    for row in 1:y_max
        xyz[idx, 1] = x_arr[col]  # x
        xyz[idx, 2] = y_arr[row]  # y
        xyz[idx, 3] = 0
        global idx += 1
    end
end

sensor_indices = Int[]

for src in sources
    # xyz is (x, y, 0) before rescaling, so match on first two columns
    local idx = findfirst(x -> isapprox(x[1], src[1]; atol=1e-6) && isapprox(x[2], src[2]; atol=1e-6), eachrow(xyz))
    if idx !== nothing
        push!(sensor_indices, idx)
    else
        println("Emitter not found in xyz: ", src)
        exit()
    end
end

for (i, idx) in enumerate(sensor_indices)
    println("Emitter $i at (n) = $(idx) maps to xyz[$idx, :] = ", xyz[idx, :])
end

"""
actual start of code
"""

if isdir("data/$object_name")
    delete_object(object_name)
end


initialize_new_object(object_name, sensor_indices, xyz, min_freq, max_freq, n_t, data_freq, debug_mode = debug_mode)

# Adds the raw data from the h5 files to the object.
# Manual implementation of mean function
function manual_mean(arr)
    sum_val = 0.0
    count = 0
    for val in arr
        sum_val += val
        count += 1
    end
    return sum_val / count
end

# Define squared mean function for RMS calculations
function mean(arr)
    sum_val = 0.0
    count = 0
    for val in arr
        sum_val += val
        count += 1
    end
    return sum_val / count
end

# Root mean square function
function rms_calc(arr)
    return sqrt(mean_squared(arr))
end

for (s, sen) in enumerate(sensor_indices)
    for n in 1:num_noise
        local h5file = h5open(raw_data_dir * "fz_homo_q_fz_src0$s.h5", "r")
        local key = "Vz"
        local measurements = read(h5file, key)[:,:,data_range]
        measurements = reshape(measurements, x_max * y_max, n_t)
        if s == 1
            global rescaler = maximum(measurements)
            global rms = sqrt(mean(measurements[sensor_indices[2], 1:div(n_t, 10)].^2))
        end
        
        # Add noise with amplitude that is noise_amplitude of the RMS signal level
        measurements ./= rescaler    
        measurements .+= rms / rescaler * noise_amplitude * randn(size(measurements)) # Add random noise
        add_raw_data(object_name, measurements, sen) 
    end
end


# Process measurements from already preprocessed data
function add_files_and_downsample(file_dir, file_prefix, defect_name)
    for src in 1:5
        local key = "data"
        local h5file = h5open(file_dir * file_prefix * "0$src.h5", "r")
        defect_data = read(h5file, key)
        for msr in 1:5
            for n in 1:num_noise
                defect = defect_data[data_range, msr]
                defect ./= rescaler 
                defect .+= rms / rescaler * noise_amplitude * randn(size(defect)) 
                add_measurements(object_name, defect_name, sensor_indices[src], sensor_indices[msr], defect)
            end
        end  
    end
end

# add measurements
for (dir, prefix, name) in zip(defect_dirs, defect_prefixes, defect_names)
    add_files_and_downsample(dir, prefix, name)
end

#calculate_greens_function(object_name, lambda = 0.000001)

#define the defect basis matrix
def_mat = zeros(Float32, n_t,2)
delta_vec = zeros(Float32, n_t)
delta_vec[1] = 1.0
delta_1 = rfft(delta_vec)
cutoff_freq = length(delta_1) ÷ 2
delta_1[(cutoff_freq+1):end] .= 0  # Low-pass filter
def_mat[:,1] = irfft(delta_1, n_t)  # Transform back to time domain
delta_2 = rfft(delta_vec)
cutoff_freq = length(delta_2) ÷ 2
delta_2[1:(cutoff_freq)] .= 0  
def_mat[:,2] = irfft(delta_2, n_t)  # Transform back to time domain

def_mat = zeros(Float32, n_t)
def_mat[1] = 1

calculate_greens_function(object_name, lambda = 0.00000001)
heatmap_3d_plot_raw_data_sensors(object_name, 2, range)

plot_greens_kernel(object_name, 2, 5106)

#analyze_source_terms(object_name, "defect_A", range)
#analyze_source_terms(object_name, "defect_A5", range)

create_reduced_system_matrix(object_name, "reduced_mat",def_mat, lambdas = [100])

color = :ice
for defect in defect_names
    evaluate_matrix(object_name, "reduced_mat", defect, color_bar = color)
end



plot_defect_sideways(object_name, "reduced_mat","defect_A5", "x", 0.5)
plot_defect_sideways(object_name, "reduced_mat","defect_A", "x", 0.5)

plot_defect_sideways(object_name, "reduced_mat","defect_A5", "y", 0.5)
plot_defect_sideways(object_name, "reduced_mat","defect_A", "y", 0.5)





exit()

#=


range = 1:500
analyze_source_terms(object_name, "defect_A", range)
analyze_source_terms(object_name, "defect_B", range)
analyze_source_terms(object_name, "defect_C", range)
analyze_source_terms(object_name, "defect_AB", range)
analyze_source_terms(object_name, "defect_AC", range)
analyze_source_terms(object_name, "defect_A5", range)


create_single_system_matrix(object_name, "single_mat_2", 2, [2], def_mat, lambdas = [1.0])
evaluate_matrix(object_name, "single_mat_2", "defect_A")


create_single_system_matrix(object_name, "single_mat_12", 1, [2], def_mat, lambdas = [1.0])
evaluate_matrix(object_name, "single_mat_12", "defect_A")


create_full_system_matrix(object_name, "full_matrix", def_mat, lambdas = [1.0])
evaluate_matrix(object_name, "full_matrix", "defect_A")
evaluate_matrix(object_name, "full_matrix", "defect_B")
evaluate_matrix(object_name, "full_matrix", "defect_C")

plot_reference_measurements(object_name,1:50)
plot_defects(object_name, "defect_A", 1:50)
plot_source_defect_difference(object_name, "defect_A", 1:50)
plot_single_location(object_name, 1, sensor_indices[1], 1)



=#




















#
#heatmap_3d_plot_raw_data(object_name, 2, range)
#heatmap_3d_plot_raw_data(object_name, 3, range)
#heatmap_3d_plot_raw_data(object_name, 4, range)
#heatmap_3d_plot_raw_data(object_name, 5, range)