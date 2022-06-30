#All have delay_time = Δt * 3500

cd(joinpath(SimPath, "Output/Meridional area increase blob"))
files = [joinpath(pwd(), "SimulationData_64_64.jld2"), joinpath(pwd(), "SimulationData_64_128.jld2"), 
        joinpath(pwd(), "SimulationData_64_256.jld2")]

#Load in the data
data64 = load(files[1])
data128 = load(files[2])
data256 = load(files[3])

t64 = time_vec(data64)
t128 = time_vec(data128)
t256 = time_vec(data256)

first_mom64 = first_moment(data64)
first_mom128 = first_moment(data128)
first_mom256 = first_moment(data256)

upperplot = Plots.plot(t64, first_mom64[:, 1],
    title = "(a) Upper layer",
    xlabel = "t",
    ylabel = "⟨A⟩",
    label = "Meridional length = 64",
    legend = :bottomright)
Plots.plot!(upperplot, t128, first_mom128[:, 1], label = "Meridional length = 128")
Plots.plot!(upperplot, t256, first_mom256[:, 1], label = "Meridional length = 256")

lowerplot = Plots.plot(t64, first_mom64[:, 2], 
    title = "(b) Lower layer",
    xlabel = "t",
    ylabel = "⟨A⟩",
    label = "Meridional length = 64",
    legend = :bottomright)
Plots.plot!(lowerplot, t128, first_mom128[:, 2], label = "Meridional length = 128")
Plots.plot!(lowerplot, t256, first_mom256[:, 2], label = "Meridional length = 256")

fullplot = Plots.plot(upperplot, lowerplot, layout = (2, 1), size = (700, 800))

savefig(fullplot, "meridareaincrease_blob.png")