#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 64, save_freq = 50, IC = GaussianStrip, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

t = time_vec(data)
mer_sec_mom = meridional_second_mom(data)

plot(t, log.(mer_sec_mom))