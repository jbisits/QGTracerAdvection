#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Domain = 128, res = 256, save_freq = 20, IC = GaussianStrip, Ensemble = false"))
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

avg_area = tracer_area_avg(data)
plot(t, avg_area)

t = time_vec(data)
area_per = tracer_area_percentile(data; conc_min = 0.1)
p1 = plot(t, area_per, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers; domain = 128, res = 256. \n Gaussian strip IC",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of 90% of log(area of tracer patch) \n in both layers; domain = 128, res = 256. \n Gaussian strip IC",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

expfit = exp_fit(data; conc_min = 0.1, tfitfinal = 52, tplot_length = 40)
linfit = linear_fit(data; conc_min = 0.1, tfitvals = [35 180], tplot_length = [0 0])
upper_area = plot(t, area_per[:, 1], 
                label = "Lower layer",
                title = "Growth of area of 90% tracer patch in upper layer; \n domain = 128, res = 256, Gaussian strip IC",
                legend = :topleft,
                lw =2,
                size = (700, 700)
                )
plot!(upper_area, linfit[:, 1, 1], linfit[:, 2, 1],
    label = "Linear fit for data points 40 to 185",
    line = (:dash, 2),
    color = :red
    )
annotate!((t[end], 0.4, text("Linear growth phase lasts ≈ 3.3 years. \n In this time the growth of 90% of the area is 80%. \n This gives diffusivity of 2.7e6m²/s.", 10, :right, :red)))
annotate!((t[end], 0.075, text("Exponetial growth lasts for ≈ 0.78 years.", 10, :right, :orange)))

phys_params = nondim2dim(data)

steps = t[35] / data["clock/dt"]
days = (steps * phys_params["Δt"]) / 3600
years = days / 365

area_per[180, 1] - area_per[35, 1]

diff = diffusivity(data, [35 180; 52 201]; conc_min = 0.1)
tracer_growth = plot(p1, logp1, upper_area, layout = @layout([a b; c]), size = (1200, 1200))
save("tracer_growth_90_dom128_stripIC.png", tracer_growth)