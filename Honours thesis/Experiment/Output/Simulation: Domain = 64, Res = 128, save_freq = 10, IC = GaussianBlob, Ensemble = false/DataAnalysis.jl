#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Domain = 64, res = 128, save_freq = 10, IC = GaussianBlob, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

#Produce histogram plots from the saved concentration data
histplots = hist_plot(data; plot_freq = 1000, xlims_same = false)

#Produce heatmaps of tacer concentration from the saved concentration data
tracerplots = tracer_plot(data; plot_freq = 1000)

#Plot heatmaps and histograms togehter.
uppertacer = plot(tracerplots[1]...)

upperhist = plot(histplots[1]...)
plot(uppertacer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[2]...)

lowerhist = plot(histplots[2]...)
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

#Time vector for plotting
t = time_vec(data)

#Create some plots of concentration diagnostics.
ConcentrationVaricance = conc_var(data)
ConcentrationMean = conc_mean(data)
SecondMoment = ConcentrationVaricance .+ ConcentrationMean.^2

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

plot(t, SecondMoment, 
        label = ["Upper Layer" "Lower layer"],
        title = "Inverse of variance of concentration \n over the lower layer grid",
        xlabel = "t",
        legend = :bottomright 
    )

plot(t, 1 ./ SecondMoment, 
        label = ["Upper Layer" "Lower layer"],
        title = "Inverse of variance of concentration \n over the upper layer grid",
        xlabel = "t",
        legend = :bottomright
    )

#Instead consider the integral ∫C²dA which is the concentration per unit area as Garrett defines
conc_int = Garrett_int(data)

plot(t, conc_int, 
        label = ["Upper layer" "Lower layer"], 
        xlabel = "t", 
        ylabel = "∑C²",
        title = "Sum of squared concentration over \n grid calculated at each time step"
    )
plot!(t, ConcentrationMean,
        label = ["Upper Layer" "Lower layer"],
        title = "Variance of concentration \n over the grid",
        xlabel = "t",
        ylabel = "Concentration"    
    )
plot(t, 1 ./ conc_int, 
        label = ["Upper layer" "Lower layer"],
        xlabel = "t", 
        ylabel = "(∑C²)⁻¹",
        title = "Inverse sum of squared concentration \n over the grid calculated at each time step", 
        legend = :bottomright
    )

ConcVsArea = concarea_animate(data; number_of_bins = 30)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)

TracerAnim = tracer_animate(data)
mp4(TracerAnim, "TracerAnim.mp4", fps = 18)

AreaVConnc = tracer_area_avg(data)

area_per = tracer_area_percentile(data; conc_min = 0.05)
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
        annotations = ([t[310]], .85, Plots.text("After approx 3.7 years \n tracer patch is size  \n of domain", :left, :green))
        )

phys_params = nondim2dim(data)

steps = t[100] / data["clock/dt"]
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
        annotations = ([t[310]], .85, Plots.text("After approx 3.7 years \n tracer patch is size  \n of domain", :left, :green))
        )