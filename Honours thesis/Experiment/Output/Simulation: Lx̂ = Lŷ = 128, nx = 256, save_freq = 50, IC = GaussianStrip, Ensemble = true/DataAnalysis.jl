cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianStrip, Ensemble = true"))

## New flow params for delay_time = Δt * 6000, seed = 1234
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
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
#savefig(fullplot, "Gaussianband_128dom_td6000.png")
## Diffusivity
ΔA² = ens_sec_mom[end, :] - ens_sec_mom[41, :]
Δt = t[end] - t[41]
Lₓ = data[1]["grid/Lx"]
K = ΔA² / (Lₓ^2 * 8 * Δt)

## Linear best fit
ΔA²/ Δt
slope = t \ ens_sec_mom

plot(t, slope .* t, label =  ["Upper best fit" "Lower best fit"], legend = :bottomright)
plot!(t, ens_sec_mom, label = ["Upper ensemble average data" "Lower ensemble average data"])
## Plots of the tracer
tracer_plots = tracer_plot(data[1]; plot_freq = 500)
upperlayerband = plot(tracer_plots[:, 1]..., size = (1400, 1400))
savefig(upperlayerband, "upperlayerbandtracer.png")
lowerlayerband = plot(tracer_plots[:, 2]..., size = (1400, 1400))
savefig(lowerlayerband, "lowerlayerbandtracer.png")
upperlayerbandIC = plot(tracer_plots[1, 1], size = (800, 400))
#savefig(upperlayerbandIC, "upperlayerbandIC.png")