include("Modules/command_center.jl")
include("Modules/custom_plotters.jl")
using HDF5
using FFTW

object_name = "Showcase"
original_n_t = 2000 # Number of time steps of original data before any changes.
data_freq = 200_000   # Frequency of the data in Hz
num_noise = 1 # Number of measurements per data set. 
# Ie it feeds num_noise sets of the same data into the algorithm with different noise patterns. 
# Set to at least 3-5 for noisy data, above 10 returns are diminishing
noise_amplitude = 0. # Amplitude of noise added to all data, as a percentage of the RMS signal over the first 10% of the data.

num_target_steps = original_n_t # Number of target steps of data. generally don't change it. It just shortens the data as [1:num_target_steps]
min_freq = 0.0  # Min frequency of bandpass filter (Due to lazy coding, it only sets frequency components to 0, so it does not increase performance to use this.)
max_freq = 100_000.0 # Max frequency of bandpass filter (data_freq/2 does not filter at all). It shortens the data in frequency space.

plotting_range = 1:100   # Range of indexes used for plotting. 1:100 on a 10ms long signal with original_n_t = 2000 = 0.5ms of data.

debug_mode = false # If true, spams print statements for every step

lambdas = [100] # List of lambda parameters used for tikhonov inversion. 
# lambas = [10,100,1000] # computes each inversion, meaning 3 times longer computation. For this dataset, this is the useful range.
# Plots the results for each in the results.
# Going higher yield diminishing returns, Going below generally breaks the results due to instability. Optimal number depends on dataset amplitudes and length. 

#Directory where the measured data is located.
data_dir = "Simulated_data/"

#Delete all Plots
if isdir("Plots/$object_name")
    rm("Plots/$object_name", recursive=true)
end

#Helper function
function mean(arr)
    sum_val = 0.0
    count = 0
    for val in arr
        sum_val += val
        count += 1
    end
    return sum_val / count
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

#read transducer locations from data
h5file = h5open(defect_A_dir * "fz_defect_A_SrcEvent01.h5", "r")
x_arr = read(h5file, "x")
y_arr = read(h5file, "y")

#x-y coordinates of transducers
src01 = [x_arr[1], y_arr[1]]
src02 = [x_arr[2], y_arr[2]]
src03 = [x_arr[3], y_arr[3]]
src04 = [x_arr[4], y_arr[4]]
src05 = [x_arr[5], y_arr[5]]

sources = [src01 , src02 , src03 , src04 , src05]

# Convert to the 1d system. Ie create a file n x 3 file that maps every n to a xyz location. 
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
        global idx += 1 # Julia for loop for some reason not sharing the same scope of this variable
    end
end

#Converting sensor xy coordinates to n- coordinates
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
Actual start of computation
"""

# Delete existing object.
if isdir("data/$object_name")
    delete_object(object_name)
end

# Initialize new object, creating a new data folder for it along with saved metadata
initialize_new_object(object_name, sensor_indices, xyz, min_freq, max_freq, num_target_steps, data_freq, debug_mode = debug_mode)


data_range = 1:num_target_steps

# Add calibration/characterization data (The scan of the full object)
for (s, sen) in enumerate(sensor_indices)
    for n in 1:num_noise
        local h5file = h5open(raw_data_dir * "fz_homo_q_fz_src0$s.h5", "r")
        local key = "Vz"
        local measurements = read(h5file, key)[:,:,data_range]
        measurements = reshape(measurements, x_max * y_max, num_target_steps)
        if s == 1
            global rescaler = maximum(measurements)
            global rms = sqrt(mean(measurements[sensor_indices[2], 1:div(num_target_steps, 10)].^2))
        end

        measurements ./= rescaler    # Rescaling helps to avoid numerical imprecision, and keeps the obtimal lambda value in the same range for different problems.
        measurements .+= rms / rescaler * noise_amplitude * randn(size(measurements)) # Add gaussian noise
        add_raw_data(object_name, measurements, sen) 
    end
end

# Function to add defect measurements
function add_files_and_downsample(file_dir, file_prefix, defect_name)
    for src in 1:5
        local key = "data"
        local h5file = h5open(file_dir * file_prefix * "0$src.h5", "r")
        defect_data = read(h5file, key)
        for msr in 1:5
            for _ in 1:num_noise
                defect = defect_data[data_range, msr]
                defect ./= rescaler 
                defect .+= rms / rescaler * noise_amplitude * randn(size(defect)) 
                add_measurements(object_name, defect_name, sensor_indices[src], sensor_indices[msr], defect)
            end
        end  
    end
end

# Add defect measurements
for (dir, prefix, name) in zip(defect_dirs, defect_prefixes, defect_names)
    add_files_and_downsample(dir, prefix, name)
end


#Define simple delta response defect matrix
def_mat = zeros(Float32, num_target_steps)
def_mat[1] = 1

#= Example on how to deine a higher dimensional defect matrix
def_mat = zeros(Float32, num_target_steps,2)
delta_vec = zeros(Float32, num_target_steps)
delta_vec[1] = 1.0
delta_1 = rfft(delta_vec)
cutoff_freq = length(delta_1) ÷ 2
delta_1[(cutoff_freq+1):end] .= 0  # Low-pass filter
def_mat[:,1] = irfft(delta_1, num_target_steps)  # Transform back to time domain
delta_2 = rfft(delta_vec)
cutoff_freq = length(delta_2) ÷ 2
delta_2[1:(cutoff_freq)] .= 0  
def_mat[:,2] = irfft(delta_2, num_target_steps)  # Transform back to time domain
=#

# Uses the calibration data to regularized fit the impulse responses.
# lambda > 0 required to assert no singularity errors
# higher lambda appears to not change the results significanty.
calculate_greens_function(object_name, lambda = 0.00000001)

#Plots the raw data as in the system (after bandpass filtering). 
#It takes a few seconds for long plotting ranges. Chosing the emitter_idx betwen 1 to 5 defines which transducer is firing
heatmap_3d_plot_raw_data_sensors(object_name, 2, plotting_range)

# Plots the impulse response function from transducer 2 to transducer 1 in real and frequency domain. 
#The target location 5106 corresponds to the trnasducer 1 location. Other trnasducer locations can be read off in the terminal output.
#plot_greens_kernel(object_name, 2, 5106)

# Plots the non-defect impulse response, defect impulse response aswell as the difference of them (The difference caused by the defect).
# Great for debugging, interesting to look at how uninteligeble the data is to the naked eye.
#analyze_source_terms(object_name, "defect_A", plotting_range)
#analyze_source_terms(object_name, "defect_A5", plotting_range)

# Builds the system matrix with a predefined sensor configuration.
# This one configures as 1 - (2-5), 2 - (3-5), 3-(4-5), 4-(5). 
# Also inverts the matrix which takes the bulk of the time.
create_reduced_system_matrix(object_name, "reduced_mat", def_mat, lambdas = lambdas)

# Evalutate the inverted matrix by feeding defect data.
color = :ice
for defect in defect_names
    evaluate_matrix(object_name, "reduced_mat", defect, color_bar = color)
end


# Sideways view plot to estimate reconstruction width
plot_defect_sideways(object_name, "reduced_mat","defect_A5", "x", 0.5)
plot_defect_sideways(object_name, "reduced_mat","defect_A", "x", 0.5)

plot_defect_sideways(object_name, "reduced_mat","defect_A5", "y", 0.5)
plot_defect_sideways(object_name, "reduced_mat","defect_A", "y", 0.5)


#=
#Pulse echo configuration computation, using sensor 2
create_single_system_matrix(object_name, "single_mat_2", 2, [2], def_mat, lambdas = [1.0])
#The results plot does not properly remove the inactive sensors.
evaluate_matrix(object_name, "single_mat_2", "defect_A")

#1-2 configuration, meaning sensor 1 fires and 2 records. 
create_single_system_matrix(object_name, "single_mat_12", 1, [2], def_mat, lambdas = [1.0])
evaluate_matrix(object_name, "single_mat_12", "defect_A")

#All sensors are active and record. Unneccesarily large due to redundant data, and for some reason less good results.
create_full_system_matrix(object_name, "full_matrix", def_mat, lambdas = [1.0])
evaluate_matrix(object_name, "full_matrix", "defect_A")
evaluate_matrix(object_name, "full_matrix", "defect_B")
evaluate_matrix(object_name, "full_matrix", "defect_C")
=#