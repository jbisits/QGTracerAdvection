#Change to the current directory
#These are done with delauy_time = Δt * 3500
cd(joinpath(SimPath, "Output/Meridional area increase band"))
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

second_mom32 = second_moment(data32)
second_mom64 = second_moment(data64)
second_mom128 = second_moment(data128)
second_mom256 = second_moment(data256)

upperplot = plot(t32, second_mom32[:, 1],
    title = "Upper layer area second moment for zonal length x̂ = 32 \n and four meridional lengths",
    xlabel = "t",
    ylabel = "σ²ₐ",
    label = "Meridional length = 32",
    legend = :bottomright)
plot!(upperplot, t64, second_mom64[:, 1], label = "Meridional length = 64")
plot!(upperplot, t128, second_mom128[:, 1], label = "Meridional length = 128")
plot!(upperplot, t256, second_mom256[:, 1], label = "Meridional length = 256")

lowerplot = plot(t32, second_mom32[:, 2], 
    title = "Lower layer area second moment for zonal length x̂ = 32 \n and four meridional lengths",
    xlabel = "t",
    ylabel = "σ²ₐ",
    label = "Meridional length = 32",
    legend = :bottomright)
plot!(lowerplot, t64, second_mom64[:, 2], label = "Meridional length = 64")
plot!(lowerplot, t128, second_mom128[:, 2], label = "Meridional length = 128")
plot!(lowerplot, t256, second_mom256[:, 2], label = "Meridional length = 256")

fullplot = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))
##
upperplot = plot(t32, second_mom32[:, 1],
    title = "Upper layer area second moment for zonal length x̂ = 32 \n and four meridional lengths",
    xlabel = "t",
    ylabel = "σ²ₐ",
    label = "Meridional length = 32",
    legend = :bottomright)
plot!(upperplot, t64[1:80], second_mom64[1:80, 1], label = "Meridional length = 64")
plot!(upperplot, t128[1:80], second_mom128[1:80, 1], label = "Meridional length = 128")
plot!(upperplot, t256[1:80], second_mom256[1:80, 1], label = "Meridional length = 256")

lowerplot = plot(t32, second_mom32[:, 2], 
    title = "Lower layer area second moment for zonal length x̂ = 32 \n and four meridional lengths",
    xlabel = "t",
    ylabel = "σ²ₐ",
    label = "Meridional length = 32",
    legend = :bottomright)
plot!(lowerplot, t64[1:80], second_mom64[1:80, 2], label = "Meridional length = 64")
plot!(lowerplot, t128[1:80], second_mom128[1:80, 2], label = "Meridional length = 128")
plot!(lowerplot, t256[1:80], second_mom256[1:80, 2], label = "Meridional length = 256")

fullplot_shortertime = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))
savefig(fullplot_shortertime, "meridareainc_band_dt3500.png")