cd(joinpath(SimPath, "Output/Square area increase band"))
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

second_mom32 = second_moment(data32) ./ 32^2
second_mom64 = second_moment(data64) ./ 64^2
second_mom128 = second_moment(data128) ./ 128^2
second_mom256 = second_moment(data256) ./ 256^2

upperplot = plot(t32, second_mom32[:, 1],
    title = "Upper layer second moment of area on square domains",
    xlabel = "t",
    ylabel = "σ²ₐ/Lₓ²",
    label = "Lx̂ = Lŷ = 32",
    legend = :topleft)
plot!(upperplot, t64, second_mom64[:, 1], label = "Lx̂ = Lŷ = 64")
plot!(upperplot, t128, second_mom128[:, 1], label = "Lx̂ = Lŷ = 128")
plot!(upperplot, t256, second_mom256[:, 1], label = "Lx̂ = Lŷ = 256")

lowerplot = plot(t32, second_mom32[:, 2], 
    title = "Lower layer second moment of area on square domains",
    xlabel = "t",
    ylabel = "σ²ₐ/Lₓ²",
    label = "Lx̂ = Lŷ = 32",
    legend = :topleft)
plot!(lowerplot, t64, second_mom64[:, 2], label =  "Lx̂ = Lŷ = 64")
plot!(lowerplot, t128, second_mom128[:, 2], label =  "Lx̂ = Lŷ = 128")
plot!(lowerplot, t256, second_mom256[:, 2], label =  "Lx̂ = Lŷ = 256")

fullplot = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))

savefig(fullplot, "sqaureareaincrease_band.png")

## Smaller plot
lowerplot = plot(t32[1:50], second_mom32[1:50, 2], 
    title = "Lower layer second moment of area on square domains",
    xlabel = "t",
    ylabel = "σ²ₐ",
    label = "Lx̂ = Lŷ = 32",
    legend = :topleft)
plot!(lowerplot, t64[1:50], second_mom64[1:50, 2], label =  "Lx̂ = Lŷ = 64")
plot!(lowerplot, t128[1:50], second_mom128[1:50, 2], label =  "Lx̂ = Lŷ = 128")
plot!(lowerplot, t256[1:50], second_mom256[1:50, 2], label =  "Lx̂ = Lŷ = 256")

## Compare the diffusivities

(second_mom64[50, 1] - second_mom64[1, 1]) / (8 * 64^2 *(t64[50] - t64[1]))

(second_mom128[end, 1] - second_mom128[1, 1]) / (8 * 128^2 *(t128[end] - t128[1]))

(second_mom256[end, 1] - second_mom256[80, 1]) / (8 * 256^2 *(t256[end] - t256[80]))