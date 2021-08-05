#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Domain = 256, res = 512, save_freq = 50, IC = GaussianBlob, Ensemble = false"))
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

t = time_vec(data; days = true)
area_per = tracer_area_percentile(data; conc_min = 0.1)
p1 = plot(t, area_per, 
        xlabel = "Days",
        label = ["Upper layer" "Lower layer"],
        title = "Growth of 90% area of tracer patch in \n both layers layers \n domain = 256, res = 512",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per),
            xlabel = "Days",
            label = ["Upper layer" "Lower layer"],
            title = "Growth of log(90% area of tracer patch) \n in both layers \n domain = 256, res = 512",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

expfit = exp_fit(data; conc_min = 0.1, tfitfinal = 40, tplot_length = 20, days = true)
linfit1 = linear_fit(data; conc_min = 0.1, tfitvals = [40 100], tplot_length = [0 0], days = true)
linfit2 = linear_fit(data; conc_min = 0.1, tfitvals = [100 170], tplot_length = [0 0], days = true)
lower_area = plot(t, area_per[:, 2], 
                label = "Lower layer",
                title = "Growth of area of tracer patch in lower layer",
                legend = :topleft,
                lw =2,
                size = (700, 700)
                )
plot!(lower_area, expfit[:, 1, 2], expfit[:, 2, 2],
    label = "Exponential fit for first 70 data points",
    line = (:dash, 2),
    color = :orange
)
plot!(lower_area, linfit1[:, 1, 2], linfit1[:, 2, 2],
    label = "Linear fit for data points 40 to 100",
    line = (:dash, 2),
    color = :red
    )  
plot!(lower_area, linfit2[:, 1, 2], linfit2[:, 2, 2],
    label = "Linear fit for data points 100 to 170",
    line = (:dash, 2),
    color = :green
    )

annotate!((t[1], 0.25, text("Exponential growth lasts \n for ≈ 33 days", 10, :left, :orange)))
annotate!((t[241], 0.25, text("Linear growth first phase \n lasts for ≈ 52 days. In this \n time the growth of 90% of \n the area is 16%. \n This gives diffusivity of 2.1e6m²/s.", 10, :right, :red)))
annotate!((t[1], 0.5, text("Linear growth second phase \n lasts for ≈ 60 days. \n In this time the growth \n of 90% of the area is 42%. \n This gives diffusivity of 4.7e6m²/s.", 10, :left, :green)))

phys_params = nondim2dim(data)
steps = t[40] / data["clock/dt"]
hours = (steps * phys_params["Δt"]) / 3600
days = hours / 24
years = days / 365

linphase1 = ((t[170] / data["clock/dt"]) * phys_params["Δt"] - (t[100] / data["clock/dt"]) * phys_params["Δt"])/ (3600 * 365)

area_per[170, 2] - area_per[100, 2]
d1 = diffusivity(data, [1 2; 40 100]; conc_min = 0.1)
d2 = diffusivity(data, [1 2; 100 170]; conc_min = 0.1)

tracer_growth = plot(p1, logp1, lower_area, layout = @layout([a b; c]), size = (1200, 1200))
save("tracer_growth_90_dom256.png", tracer_growth)