cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianStrip, Ensemble = true"))

## New flow params for delay_time = Δt * 6000, 
# first 10 seed 1234, 10-20 seed 4321
data = Array{Dict{String, Any}}(undef, 20)
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
                label = "Ensemble member",
                legend = :bottomright)
lowerplot = plot(t, sec_mom[:, 2, 1],
                title = "Lower layer second moment of area growth \n for Gaussian band initial condition",
                xlabel = "t",
                ylabel = "σ²ₐ",
                label = "Ensemble member",
                legend = :bottomright)
for i ∈ 2:length(data)
    plot!(upperplot, t, sec_mom[:, 1, i], label = false) 
    plot!(lowerplot, t, sec_mom[:, 2, i], label = false) 
end

ens_conc = ensemble_concentration(data)
ens_sec_mom = second_moment(ens_conc)

plot!(upperplot, t, ens_sec_mom[:, 1], label = "Ensemble", line = (:dash, :black, 2))
plot!(lowerplot, t, ens_sec_mom[:, 2], label = "Ensemble", line = (:dash, :black, 2))

fullplot = plot(upperplot, lowerplot, layout = (2, 1), size = (800, 800))
#savefig(fullplot, "Gaussianband_128dom_td6000.png")

## Diffusivity
ΔA² = ens_sec_mom[61, :] - ens_sec_mom[21, :]
Δt = t[61] - t[21]
Lₓ = data[1]["grid/Lx"]
K = ΔA² / (Lₓ^2 * 8 * Δt)

nondim2dim(data[1])
## Linear best fit
fit = [ones(length(t)) t] \ ens_sec_mom

upperlinfit = plot(t, fit[1, 1] .+ fit[2, 1] .* t, 
                    title = "Upper layer linear best fit of the growth of \n the ensemble second moment of area",
                    label = "Best fit of ensemble data", 
                    xlabel = "t",
                    ylabel = "σ²ₐ",
                    legend = :bottomright, 
                    line = (:dash))
plot!(upperlinfit, t, ens_sec_mom[:, 1], 
    label = "Ensemble data")
savefig(upperlinfit, "upperlinfitband.png")

lowerlinfit = plot(t, fit[1, 2] .+ fit[2, 2] .* t, 
                    title = "Lower layer linear best fit of the growth of \n the ensemble second moment of area",
                    label = "Best fit of ensemble data", 
                    xlabel = "t",
                    ylabel = "σ²ₐ",
                    legend = :bottomright, 
                    line = (:dash))
plot!(lowerlinfit, t, ens_sec_mom[:, 2], 
    label = "Ensemble data")
savefig(lowerlinfit, "lowerlinfitband.png")

Lₓ = data[1]["grid/Lx"]
K_linfit = fit[2, :] ./ (Lₓ^2 * 8)

## Plots of the tracer
tracer_plots = tracer_plot(data[1]; plot_freq = 500)
upperlayerband = plot(tracer_plots[:, 1]..., size = (1400, 1400))
savefig(upperlayerband, "upperlayerbandtracer.png")
lowerlayerband = plot(tracer_plots[:, 2]..., size = (1400, 1400))
savefig(lowerlayerband, "lowerlayerbandtracer.png")
upperlayerbandIC = plot(tracer_plots[1, 1], size = (800, 400))
#savefig(upperlayerbandIC, "upperlayerbandIC.png")

## Plots of ensemble Concentration
ens_plots = tracer_plot(ens_conc; plot_freq = 500)
plot(ens_plots[:, 1]..., size = (1400, 1400))
plot(ens_plots[:, 2]..., size = (1400, 1400))