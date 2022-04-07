#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 256, nx = 512, save_freq = 100, IC = GaussianBlob, Ensemble = true"))

## Load in the data for delay_time = Δt * 6000 with the new parameters on larger domain.
data = Array{Dict{String, Any}}(undef, 25)
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