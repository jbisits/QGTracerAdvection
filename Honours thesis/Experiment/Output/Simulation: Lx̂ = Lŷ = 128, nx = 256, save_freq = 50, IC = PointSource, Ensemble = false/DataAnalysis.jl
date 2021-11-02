cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = PointSource, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

## Look at first moment.
t = time_vec(data)
first_mom = first_moment(data)
pt_release = plot(t, first_mom, 
                xlabel = "t",
                ylabel = "⟨A⟩",
                title = "Average area growth from a point release tracer",
                label = ["Upper layer" "Lower layer"], 
                legend = :topleft)

savefig(pt_release, "pt_release.png")
lin_phase = plot(t[50:100], first_mom[50:100, 1], 
                    xlabel = "t", 
                    ylabel = "⟨A⟩",
                    title = "Garrett stage 3 growth phase of \npoint release in upper layer",
                    label = false)

savefig(lin_phase, "lin_phase_pt_release.png")
ΔA = first_mom[100, :] .- first_mom[50, :]
Δt = t[100] - t[50]

K = ΔA ./ (4 * π * Δt)

dims = nondim2dim(data)

K_dim = K * dims["Ld"] * 0.02

##
tracer_plots = tracer_plot(data; plot_freq = 50)
