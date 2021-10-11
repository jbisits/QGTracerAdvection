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

upperplot = plot(t64, first_mom64[:, 1],
    title = "Upper layer average area for zonal length x̂ = 64 \n and three meridional lengths",
    xlabel = "t",
    ylabel = "⟨A⟩",
    label = "Meridional length = 64",
    legend = :bottomright)
plot!(upperplot, t128, first_mom128[:, 1], label = "Meridional length = 128")
plot!(upperplot, t256, first_mom256[:, 1], label = "Meridional length = 256")

lowerplot = plot(t64, first_mom64[:, 2], 
    title = "Lower layer average area for zonal length x̂ = 64 \n and three meridional lengths",
    xlabel = "t",
    ylabel = "⟨A⟩",
    label = "Meridional length = 64",
    legend = :bottomright)
plot!(lowerplot, t128, first_mom128[:, 2], label = "Meridional length = 128")
plot!(lowerplot, t256, first_mom256[:, 2], label = "Meridional length = 256")

fullplot = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))

#savefig(fullplot, "meridareainc_blob_dt6000.png")