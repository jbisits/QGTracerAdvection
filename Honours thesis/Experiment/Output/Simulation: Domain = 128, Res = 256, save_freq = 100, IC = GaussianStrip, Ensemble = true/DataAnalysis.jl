cd(joinpath(SimPath, "Output/Simulation: Domain = 128, res = 256, save_freq = 100, IC = GaussianStrip, Ensemble = true"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data. This is an ensemble simulation so now have an array of dictionaries.
data = Array{Dict{String, Any}}(undef, 2)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end