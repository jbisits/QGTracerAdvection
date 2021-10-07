cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 64, nx = 128, save_freq = 50, IC = GaussianStrip, Ensemble = true"))

## Load in the data for delay_time = Δt * 3000
data = Array{Dict{String, Any}}(undef, 15)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end
## Load in the data for delay_time = Δt * 4500
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData_25.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 24)*".jld2")
        data[i] = load(file)
    end
end
## Load in the data for delay_time = Δt * 5000
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData_15.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 14)*".jld2")
        data[i] = load(file)
    end
end
## New flow params for delay_time = Δt * 6000
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData_35.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i + 34)*".jld2")
        data[i] = load(file)
    end
end

##
t = time_vec(data[1])
sec_mom = second_moment(data)

upperplot = plot(t, sec_mom[:, 1, 1],
                title = "Upper layer second moment of area growth \n for Gaussian band initial condition",
                xlabel = "t",
                ylabel = "σ²ₐ",
                label = "Member 1",
                legend = :bottomright)
lowerplot = plot(t, sec_mom[:, 2, 1],
                title = "Lower layer second moment of area growth \n for Gaussian band initial condition",
                xlabel = "t",
                ylabel = "σ²ₐ",
                label = "Member 1",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(upperplot, t, sec_mom[:, 1, i], label = "Member "*string(i)) 
    plot!(lowerplot, t, sec_mom[:, 2, i], label = "Member "*string(i)) 
end

ens_conc = ensemble_concentration(data)
ens_sec_mom = second_moment(ens_conc)

plot!(upperplot, t, ens_sec_mom[:, 1], label = "Ensemble", line = (:dash, :black, 2))
plot!(lowerplot, t, ens_sec_mom[:, 2], label = "Ensemble", line = (:dash, :black, 2))

fullplot = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))
#savefig(fullplot, "Gaussianband_64dom_td4500.png")
##
ΔA² = ens_sec_mom[50, :] - ens_sec_mom[1, :]
Δt = t[50] - t[1]
K = ΔA² / (data[1]["grid/Lx"]^2 * 8 * Δt)

##
ΔA²/ Δt
slope = t[1:50] \ ens_sec_mom[1:50, :]

plot(t[1:50], slope .* t[1:50], label =  ["Upper best fit" "Lower best fit"], legend = :bottomright)
plot!(t[1:50], ens_sec_mom[1:50, :], label = ["Upper ensemble average data" "Lower ensemble average data"])
##
tracer_plots = tracer_plot(data[1]; plot_freq = 500)
plot(tracer_plots[:, 1]..., size = (1200, 1200))
plot(tracer_plots[:, 2]..., size = (1200, 1200))
#########################################################################################
##
t = time_vec(data[1])
mer_sec_mom = meridional_second_mom(data)

upperplot = plot(t, mer_sec_mom[:, 1, 1],
                title = "Upper layer",
                xlabel = "t",
                ylabel = "σ²y",
                label = "Ensemble member",
                legend = :bottomright)
lowerplot = plot(t, mer_sec_mom[:, 2, 1],
                title = "Lower layer",
                xlabel = "t",
                ylabel = "σ²y",
                label = "Ensemble member",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(upperplot, t, mer_sec_mom[:, 1, i], label = false) 
    plot!(lowerplot, t, mer_sec_mom[:, 2, i], label = false) 
end

ens_conc = ensemble_concentration(data)
ens_mer_sec_mom = meridional_second_mom(ens_conc)

plot!(upperplot, t, ens_mer_sec_mom[:, 1], label = "Ensemble", line = (:dash, :black, 2))
plot!(lowerplot, t, ens_mer_sec_mom[:, 2], label = "Ensemble", line = (:dash, :black, 2))

plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))