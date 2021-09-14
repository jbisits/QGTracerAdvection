#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Domain = 128, res = 256, save_freq = 200, IC = GaussianBlob, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

#Produce histogram plots from the saved concentration data
histplots = hist_plot(data; plot_freq = 3000)

#Produce heatmaps of tacer concentration from the saved concentration data
tracerplots = tracer_plot(data; plot_freq = 3000)

#Plot heatmaps and histograms togehter.
uppertacer = plot(tracerplots[:, 1]...)

upperhist = plot(histplots[:, 1]...)
plot(uppertacer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[:, 2]...)

lowerhist = plot(histplots[:, 2]...)
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

#Time vector for plotting
t = time_vec(data)

#Create some plots of concentration diagnostics.
ConcentrationVaricance = conc_var(data)
ConcentrationMean = conc_mean(data)

meanplot = plot(t, ConcentrationMean, 
                    label = ["Upper Layer" "Lower layer"],
                    title = "Mean concentration \n over the grid",
                    xlabel = "t",
                    ylabel = "Concentration",
                    ylims = (0, 0.001)
                )
varplot = plot(t, ConcentrationVaricance, 
                    label = ["Upper Layer" "Lower layer"],
                    title = "Variance of concentration \n over the grid",
                    xlabel = "t",
                    ylabel = "Concentration"    
                )
plot(meanplot, varplot, size = (1000, 600))

ConcVsArea = concarea_animate(data)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)

area_per = tracer_area_percentile(data; Cₚ = 0.5)
plot(t, area_per, 
    label = ["Upper layer" "Lower layer"],
    title = "Growth of area of tracer patch in both layers layer",
    legend = :topleft
    )

plot(t, area_per[:, 1],
     label = "Upper layer",
     title = "Growth of area of tracer patch in the upper layer",
     legend = :topleft
     )
scatter!([t[230]], [area_per[230, 1]],
        label = "Stage 2 -> stage 3",
        annotations = (13, area_per[230, 1], Plots.text("Stage three begins after \n approximately 2.6 years", :left, :orange))
        )
scatter!([t[315]], [area_per[315, 1]],
        label = "Tracer patch ≈ size of domain",
        annotations = ([t[310]], 3.5, Plots.text("After approx 3.7 years \n tracer patch is size  \n of domain", :left, :green))
        )


phys_params = nondim2dim(data)

steps = t[251] / data["clock/dt"]
days = (steps * phys_params["Δt̂"]) / 3600
years = days / 365

plot(t, area_per[:, 2],
     label = "Lower layer",
     title = "Growth of area of tracer patch in the Lower layer",
     legend = :topleft,
     color = :red
     )
scatter!([t[251]], [area_per[251, 2]],
        label = "Stage 2 -> stage 3",
        annotations = (13, area_per[230, 1], Plots.text("Stage three begins after \n approximately 2.8 years", :left, :orange))
        )
scatter!([t[315]], [area_per[315, 1]],
        label = "Tracer patch ≈ size of domain",
        annotations = ([t[310]], 3.5, Plots.text("After approx 3.7 years \n tracer patch is size  \n of domain", :left, :green))
        )

## Look at whether the second_mom function produces reasonable output and diffusivity

tsecs = time_vec(data; time_measure = "secs")
tdays = time_vec(data; time_measure = "days")
avg_area = tracer_avg_area(data)
second_moments = tracer_second_mom(data)

plot(t, avg_area)
plot(t, second_moments)

K = [ (second_moments[i + 1, 2] - second_moments[i, 2]) ./ (2 * (t[i + 1] - t[i])) for i ∈ 1:length(t) - 1]
K = second_moments ./ (2 .* t)
plot(tdays, K)