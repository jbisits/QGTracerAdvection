#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

## Load in the data for delay_time = Δt * 6000 with the new parameters, seed = 1234
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
first_moms = first_moment(data)
first_mom_upper = plot(t, first_moms[:, 1, 1], 
                        label = "Member 1", 
                        title = "Upper layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
first_mom_lower = plot(t, first_moms[:, 2, 1], 
                        label = "Member 1", 
                        title = "Lower layer average area growth for Gaussian blob initial condition",
                        xlabel = "t",
                        ylabel = "⟨A⟩",
                        legend = :topleft)
for i ∈ 2:length(data)
    plot!(first_mom_upper, t, first_moms[:, 1, i], label = "Memeber "*string(i))
    plot!(first_mom_lower, t, first_moms[:, 2, i], label = "Memeber "*string(i))
end

ensemble_conc = ensemble_concentration(data)
ensemble_avg = first_moment(ensemble_conc)

plot!(first_mom_upper, t, ensemble_avg[:, 1], label = "Ensemble average", line = (:dash, 2, :black))
plot!(first_mom_lower, t, ensemble_avg[:, 2], label = "Ensemble average", line = (:dash, 2, :black))

fullplot = plot(first_mom_upper, first_mom_lower, layout = (2, 1), size= (800, 800))
#savefig(fullplot, "Gaussianblob_128dom_td6000.png")
##
Δt = t[end] - t[1]
ΔA = ensemble_avg[end, :] .- ensemble_avg[1, :]
K = ΔA ./ (4 * π * Δt)

##
slope = t \ ensemble_avg

plot(t, slope[1] .* t, label = ["Upper best fit" "Lower best fit"], legend = :bottomright)
plot!(t, ensemble_avg[:, 1], label = ["Upper ensemble average data" "Lower ensemble average data"])
##
##
tracer_plots = tracer_plot(data[1]; plot_freq = 500)
upperlayerblob = plot(tracer_plots[:, 2]..., size = (1400, 1400))
#plot(tracer_plots[1, 1], tracer_plots[2, 1], tracer_plots[3, 1], tracer_plots[4, 1], tracer_plots[5, 1], size = (1400, 1400))
savefig(upperlayerblob, "upperlayertracerblob.png")
upperlayerblobIC = plot(tracer_plots[1, 1], size = (900, 400))
savefig(upperlayerblobIC, "blobIC.png")
plot(tracer_plots[:, 2], size = (1200, 1200))