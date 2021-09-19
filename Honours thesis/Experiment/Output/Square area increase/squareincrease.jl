cd(joinpath(SimPath, "Output/Square area increase"))
files = [joinpath(pwd(), "SimulationData_32.jld2"), joinpath(pwd(), "SimulationData_64.jld2"),
        joinpath(pwd(), "SimulationData_128.jld2"), joinpath(pwd(), "SimulationData_256.jld2")]

#Load in the data
data32 = load(files[1])
data64 = load(files[2])
data128 = load(files[3])
data256 = load(files[4])

t32 = time_vec(data32)
t64 = time_vec(data64)
t128 = time_vec(data128)
t256 = time_vec(data256)

first_mom32 = first_moment(data32)
first_mom64 = first_moment(data64)
first_mom128 = first_moment(data128)
first_mom256 = first_moment(data256)

upperplot = plot(t32, first_mom32[:, 1],
    title = "Upper layer average area on square domains",
    xlabel = "t",
    ylabel = "⟨A⟩",
    label = "Lx̂ = Lŷ = 32",
    legend = :topleft)
plot!(upperplot, t64, first_mom64[:, 1], label = "Lx̂ = Lŷ = 64")
plot!(upperplot, t128, first_mom128[:, 1], label = "Lx̂ = Lŷ = 128")
plot!(upperplot, t256[1:round(Int, 3 * end / 4)], first_mom256[1:round(Int, 3 * end / 4), 1], label = "Lx̂ = Lŷ = 256")

lowerplot = plot(t32, first_mom32[:, 2], 
    title = "Lower layer average area on square domains",
    xlabel = "t",
    ylabel = "⟨A⟩",
    label = "Lx̂ = Lŷ = 32",
    legend = :topleft)
plot!(lowerplot, t64, first_mom64[:, 2], label =  "Lx̂ = Lŷ = 64")
plot!(lowerplot, t128, first_mom128[:, 2], label =  "Lx̂ = Lŷ = 128")
plot!(lowerplot, t256[1:round(Int, 3 * end / 4)], first_mom256[1:round(Int, 3 * end / 4), 2], label =  "Lx̂ = Lŷ = 256")

fullplot = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))

savefig(fullplot, "sqaureareaincrease.png")