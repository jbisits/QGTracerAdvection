#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 256, nx = 512, save_freq = 100, IC = GaussianBlob, Ensemble = true"))

## Load in the data for delay_time = Δt * 6000 with the new parameters on larger domain.
data = Array{Dict{String, Any}}(undef, 50)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end

# Not enough ram to do this all at once so need to change the way things are computed
# Will use `Glob.jl` to access all the .jld2 files, then compute a diffusivity over the 
# relevant time interval for each member and create an ensemble average field at the end.

using Glob

data_files = glob("*.jld2")
data_files = data_files[1:50]

## First moments of ensemble members and ensemble concentration

ens_conc = load(data_files[1])
save_freq = ens_conc["save_freq"]
nsteps = ens_conc["clock/nsteps"]
saved_vals = 0:save_freq:nsteps

t = time_vec(ens_conc)

first_moms = Array{Float64}(undef, length(t), 2, length(data_files))
i = 1
for file ∈ data_files

    println("$file")
    data = load(file)
    first_moms[:, :, i] = first_moment(data)
    
    if i != 1

        for j ∈ saved_vals
            @. ens_conc["snapshots/Concentration/"*string(j)] += data["snapshots/Concentration/"*string(j)]
        end

    elseif i == length(data_files)

        for j ∈ saved_vals
            ens_conc["snapshots/Concentration/"*string(j)]  ./=  length(data_files)
        end

    end

    i += 1

end

ens_first_moms = first_moment(ens_conc)

ens_avg_lf = [ones(length(t)) t] \ ens_first_moms 

plot(t, first_moms[:, :, 50], label = false)
plot(t, ens_first_moms[:, 1], legend = :topleft, label = "Ensemble average")
plot!(t, ens_avg_lf[1, 1] .+ t .* ens_avg_lf[2, 1], label = "Linear fit")

#file = "saved_data.jld2"

jldopen(file, "a+") do path
    path["First_moms/t"] = t
    path["First_moms/Ensemble_first_moms"] = first_moms
    path["First_moms/Ensemble_avg_first_moms"] =  ens_first_moms
    path["First_moms/Ensemble_avg_lf"] = ens_avg_lf
end 

# Diffusivity for ensemble members and ensemble average from linear fits 
# of the data during the 20-80 non dimensional time frame. This removes any 
# initial mixing growth and avoid s any plataueing beacuse of the domain size.

ense_diff = (ens_avg_lf[2, :] ./ (4*π)) .* 0.02 .* dims["Ld"]
t̂_16_5  =  findfirst(t .>= 16.5)
t̂_80 = findfirst(t .>= 80) - 1
linear_phase = t̂_20:t̂_80
t̂_linear = t[linear_phase]
plot(t[linear_phase], first_moms[linear_phase, :, 1], label = false)
t̂_ones = ones(length(t̂_linear))
t_fit_mat = [t̂_ones t̂_linear]
member_diffs = Array{Float64}(undef, length(first_moms[1, 1, :]), 2)
for i ∈ 1:length(first_moms[1, 1, :])

    linfit = t_fit_mat \ first_moms[linear_phase, :, i]
    member_diffs[i, :] = (linfit[2, :] ./ (4*π)) .* 0.02 .* dims["Ld"]

end
member_diffs
jldopen(file, "a+") do path
    path["Diffusivity/member_diffs"] = member_diffs
    path["Diffusivity/ens_avg_diff"] = ense_diff
end

## Bootstrap the samples
## Bootstrap the ensemble members to form error for ensemble average

N = 1000 
sample = 30
sample_vec = Array{Int64}(undef, sample)
sample_vals = 1:length(data_files)
ens_diff = Array{Float64}(undef, N, 2)
t_mat = [ones(length(t)) t]

for n ∈ 1:N

    println("Bootstrap "*string(n))
    sample_vec = StatsBase.sample(sample_vals, sample, replace = false)
    data_sample = data_files[sample_vec]
    ens_conc_sample = load(data_sample[1])
    for i ∈ 2:length(data_sample)

        data = load(data_sample[i])
        
        for j ∈ saved_vals
        
            @. ens_conc_sample["snapshots/Concentration/"*string(j)] += data["snapshots/Concentration/"*string(j)]
        
        end
        if i == length(data_sample)

            for j ∈ saved_vals
            
                ens_conc_sample["snapshots/Concentration/"*string(j)]  ./=  length(data_sample)

            end

        end

    end
    ensemble_avg_sample = first_moment(ens_conc_sample)
    ens_fit_sample = t_mat \ ensemble_avg_sample
    ens_diff[n, :] =  ens_fit_sample[2, :]   

end

ens_diff = (ens_diff ./ (4*π)) .* 0.02 .* dims["Ld"]

jldopen("saved_data.jld2", "a+") do path
    path["Bootstrap/diff_samples_v2"] = ens_diff
end

ens_diff
## Spatial subset
# Info needed
ens_conc = load(data_files[1])
save_freq = ens_conc["save_freq"]
nsteps = ens_conc["clock/nsteps"]
saved_vals = 0:save_freq:nsteps
dims = nondim2dim(ens_conc)

t = time_vec(ens_conc)
t̂_16_5 =  findfirst(t .>= 16.5)
t̂_80 = findfirst(t .>= 80) - 1
linear_phase = t̂_16_5:t̂_80
t̂_linear = t[linear_phase]
t̂_ones = ones(length(t̂_linear))
t_fit_mat = [t̂_ones t̂_linear]

## Spatial subsets, full temporal resolution
zonal_subset = 2 .^ (0:8)
meridional_subset = 2 .^ (0:8)

member_diffs_ss = Array{Float64}(undef, length(data_files), 2, length(zonal_subset) * length(meridional_subset))
k = 1
for i ∈ zonal_subset, j ∈ meridional_subset
    
    # Cannot load all at once so need to use a loop here fir tge first moments
    first_moms = Array{Float64}(undef, length(t), 2, length(data_files))
    l = 1
    for file ∈ data_files

        data = load(file)
        first_moms[:, :, l] = first_moment(data, i, j)
        l += 1

    end
    member_subset_linfit = [t_fit_mat \ first_moms[linear_phase, :, k] for k ∈ 1:length(data_files)]
    member_subset_linfit = [[member_subset_linfit[k][2, 1] for k ∈ 1:length(data_files)] [member_subset_linfit[k][2, 2] for k ∈ 1:length(data_files)]]
    K_subset_member_linfit = member_subset_linfit ./ (4π)
    member_diffs_ss[:, :, k] = @. K_subset_member_linfit * dims["Ld"] * 0.02
    k += 1

end
member_diffs_ss

# RMS error compared to diffusivity from ensemble average concentration field
member_diffs_ss_rms_err = Array{Float64}(undef, 1, 2, length(member_diffs_ss[1, 1, :]))

for i ∈ 1:length(member_diffs_ss_v2[1, 1, :])

    member_diffs_ss_rms_err[:, 1, i] = sqrt.( mean((member_diffs_ss[:, 1, i] .- ense_diff[1]).^ 2, dims = 1 ) )
    member_diffs_ss_rms_err[:, 2, i] = sqrt.( mean((member_diffs_ss[:, 2, i] .- ense_diff[2]).^ 2, dims = 1 ) )

end

member_diffs_ss_rms_err

jldopen(file, "a+") do path
    path["Spatial_subset/RMS_error"] = member_diffs_ss_rms_err
end

## Temporal subsets, full spatial resolution
# Info needed
ens_conc = load(data_files[1])
save_freq = ens_conc["save_freq"]
nsteps = ens_conc["clock/nsteps"]
saved_vals = 0:save_freq:nsteps
dims = nondim2dim(ens_conc)

t = time_vec(ens_conc)

t̂_16_5 =  findfirst(t .>= 16.5)
t̂_80 = findfirst(t .>= 80) - 1

time_inc = @. 2^(0:6)
time_subset_linfits = Array{Float64}(undef, length(data_files), 2, length(time_inc))
# Cannot load all at once so need to use a loop here fir tge first moments
first_moms = Array{Float64}(undef, length(t), 2, length(data_files))
l = 1
for file ∈ data_files

    data = load(file)
    first_moms[:, :, l] = first_moment(data)
    l += 1

end
j = 1
for i ∈ time_inc

    @info "Time subset $i"
    # Create the time subset
    linear_phase = t̂_16_5:i:t̂_80
    @info "Linear phase $linear_phase"
    t̂_linear = t[linear_phase]
    t̂_ones = ones(length(t̂_linear))
    t_fit_mat = [t̂_ones t̂_linear]
    member_first_moms_linfit = [t_fit_mat \ first_moms[linear_phase, :, k] for k ∈ 1:length(data_files)]
    member_first_moms_linfit = [[member_first_moms_linfit[k][2, 1] for k ∈ 1:length(data_files)] [member_first_moms_linfit[k][2, 2] for k ∈ 1:length(data_files)]]
    K_subset_member_lifit = member_first_moms_linfit ./ (4π)
    time_subset_linfits[:, :, j] .= K_subset_member_lifit .* dims["Ld"] .* 0.02
    j += 1

end
time_subset_linfits

ense_diff = load("saved_data.jld2")["Diffusivity/ens_avg_diff"]

time_RMS = Array{Float64}(undef, 1, 2, length(time_subset_linfits[1, 1, :]))
for i ∈ 1:length(time_subset_linfits[1, 1, :])
    time_RMS[:, 1, i] .= sqrt.( mean((time_subset_linfits[:, 1, i] .- ense_diff[1]).^ 2, dims = 1 ) )
    time_RMS[:, 2, i] .= sqrt.( mean((time_subset_linfits[:, 2, i] .- ense_diff[2]).^ 2, dims = 1 ) )
end
time_RMS

jldopen("saved_data.jld2", "a+") do path
    path["Temporal_subset/RMS_error"] = time_RMS
end

## Spatio-temporal subsets
time_inc = @. 2^(0:6)
spatial_subset = @. 2^(0:8)

time_spatial_subset_linfits = Array{Union{Missing, Float64}}(missing, length(data_files), 2, length(time_inc), length(spatial_subset))
n = 1
for m ∈ spatial_subset

    @info "Spatial subset $m"
    # Cannot load all at once so need to use a loop here fir tge first moments
    first_moms = Array{Float64}(undef, length(t), 2, length(data_files))
    l = 1
    for file ∈ data_files

        data = load(file)
        first_moms[:, :, l] = first_moment(data, m, m)
        l += 1

    end

    k = 1
    for i ∈ time_inc

        # Create the time subset
        linear_phase = t̂_16_5:i:t̂_80
        @info "Time subset $i"
        t̂_linear = t[linear_phase]
        t̂_ones = ones(length(t̂_linear))
        t_fit_mat = [t̂_ones t̂_linear]

        member_first_moms_linfit = [t_fit_mat \ first_moms[linear_phase, :, k] for k ∈ 1:length(data_files)]
        member_first_moms_linfit = [[member_first_moms_linfit[k][2, 1] for k ∈ 1:length(data_files)] [member_first_moms_linfit[k][2, 2] for k ∈ 1:length(data_files)]]
        temp = member_first_moms_linfit ./ (4π)
        time_spatial_subset_linfits[:, :, k, n] = @. temp * dims["Ld"] * 0.02

        k += 1

    end

    n += 1

end
time_spatial_subset_linfits
ts_rms_err = Array{Float64}(undef, 1, 2, length(time_inc), length(spatial_subset))

for j ∈ 1:length(spatial_subset)

    for i ∈ 1:length(time_inc)

        diffs = [collect(skipmissing(time_spatial_subset_linfits[:, 1, i, j])) collect(skipmissing(time_spatial_subset_linfits[:, 2, i, j]))]
        ts_rms_err[:, 1, i, j] = sqrt.( mean((diffs[:, 1] .-  ense_diff[1]).^ 2, dims = 1 ) )
        ts_rms_err[:, 2, i, j] = sqrt.( mean((diffs[:, 2] .-  ense_diff[2]).^ 2, dims = 1 ) )

    end

end

ts_rms_err

jldopen("saved_data.jld2", "a+") do path
    path["Spatio_temp_subset/RMS_error"] = ts_rms_err
end
