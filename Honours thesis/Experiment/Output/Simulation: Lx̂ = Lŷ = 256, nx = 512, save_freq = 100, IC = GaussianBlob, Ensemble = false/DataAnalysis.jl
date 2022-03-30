data_path = joinpath(pwd(), "Output/Simulation: Lx̂ = Lŷ = 256, nx = 512, save_freq = 100, IC = GaussianBlob, Ensemble = false/SimulationData.jld2")

data = load(data_path)

t = time_vec(data[1])
nsteps = data["clock/nsteps"]
save_freq = data["save_freq"]
t = [data["snapshots/t/"*string(i)] for i in 0:save_freq:nsteps]
first_moms = first_moment(data)

plot(t, first_moms,label = false)

α = (first_moms[end, :] - first_moms[40, :]) ./ (t[end] - t[40])

diff_nondim = α ./ (4π)

dims = nondim2dim(data)

diff_dim = diff_nondim .* 0.02 .* dims["Ld"]


days = (((t[end] / 0.005) * dims["Δt"]) / 3600) / 24

years = days / 365