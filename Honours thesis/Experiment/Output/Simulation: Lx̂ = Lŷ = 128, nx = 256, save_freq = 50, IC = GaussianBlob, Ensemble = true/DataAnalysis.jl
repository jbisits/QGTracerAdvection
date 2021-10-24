#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

## Load in the data for delay_time = Δt * 6000 with the new parameters, 
#seed = 1234 for first 10, seed = 4321 for 10-20, seed = 2341 for 20-30
data = Array{Dict{String, Any}}(undef, 30)
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
first_moms = first_moment(data)
first_mom_upper = plot(t, first_moms[:, 1, 1], 
                        label = "Ensemble member", 
                        title = "Upper layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
first_mom_lower = plot(t, first_moms[:, 2, 1], 
                        label = "Ensemble member", 
                        title = "Lower layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
for i ∈ 2:length(data)
    plot!(first_mom_upper, t, first_moms[:, 1, i], label = false)
    plot!(first_mom_lower, t, first_moms[:, 2, i], label = false)
end

ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

plot!(first_mom_upper, t, ensemble_avg[:, 1], label = "Ensemble average", line = (:dash, 2, :black))
plot!(first_mom_lower, t, ensemble_avg[:, 2], label = "Ensemble average", line = (:dash, 2, :black))

fullplot = plot(first_mom_upper, first_mom_lower, layout = (2, 1), size= (800, 800))
#savefig(fullplot, "Gaussianblob_128dom_td6000.png")
## Diffusivity
Δt = t[41] - t[1]
ΔA = ensemble_avg[41, :] .- ensemble_avg[1, :]
K = ΔA ./ (4 * π * Δt)

## Linear fit
fit = [ones(length(t)) t] \ ensemble_avg

upperlinfit = plot(t, fit[1, 1] .+ fit[2, 1] .* t, 
                    title = "Upper layer linear best fit of the growth of \n the ensemble average area",
                    label = "Best fit of ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(upperlinfit, t, ensemble_avg[:, 1], 
    label = "Ensemble data")
#savefig(upperlinfit, "upperlinfitblob.png")

lowerlinfit = plot(t, fit[1, 2] .+ fit[2, 2] .* t, 
                    title = "Lower layer linear best fit of the growth of \n the ensemble average area",
                    label = "Best fit of ensemble data", 
                    xlabel = "t",
                    ylabel = "⟨A⟩",
                    legend = :bottomright, 
                    line = (:dash))
plot!(lowerlinfit, t, ensemble_avg[:, 2], 
    label = "Ensemble data")
#savefig(lowerlinfit, "lowerlinfitblob.png")

K_linfit = fit[2, :] ./ (4 * π)

dims = nondim2dim(data[1])
K_linfit_dim = @. K_linfit * dims["Ld"] * 0.02

## Tracer plots
tracer_plots = tracer_plot(data[1]; plot_freq = 500)
upperlayerblob = plot(tracer_plots[:, 1]..., size = (1400, 1400))
savefig(upperlayerblob, "upperlayertracerblob.png")
lowerlayerblob = plot(tracer_plots[:, 2]..., size = (1400, 1400))
savefig(lowerlayerblob, "lowerlayertracerblob.png")
upperlayerblobIC = plot(tracer_plots[1, 1], size = (900, 400))
savefig(upperlayerblobIC, "blobIC.png")

## Ensemble plots
ens_plots = tracer_plot(ensemble_conc; plot_freq = 500)
upperblobens = plot(ens_plots[:, 1]..., size = (1400, 1400))
savefig(upperblobens, "upperblobens.png")
lowerblobens = plot(ens_plots[:, 2]..., size = (1400, 1400))
savefig(lowerblobens, "lowerblobens.png")