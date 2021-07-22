#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Domain = 32, res = 64, save_freq = 10, IC = GaussianBlob, Ensemble = false"))
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

t = time_vec(data)
area_per = tracer_area_percentile(data; conc_min = 0.1)
p1 = plot(t, area_per, 
    label = ["Upper layer" "Lower layer"],
    title = "Growth of 90% area of tracer patch in both layers \n domain = 32, res = 64",
    legend = :topleft
    )
logp1 =  plot(t, log.(area_per), 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of log(90% area of tracer patch) in both layers \n domain = 32, res = 64",
        legend = :topleft
        )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

expfit = exp_fit(data; conc_min = 0.1, tfitfinal = 201, tplot_length = 0)
linfit1 = linear_fit(data; conc_min = 0.1, tfitvals = [201 320], tplot_length = [0 0])
linfit2 = linear_fit(data; conc_min = 0.1, tfitvals = [320 465], tplot_length = [0 0])
lower_area = plot(t, area_per[:, 2], 
                label = "Lower layer",
                title = "Growth of area of 90% tracer patch in lower layer",
                legend = :topleft,
                lw =2,
                size = (700, 700)
                )
plot!(lower_area, expfit[:, 1, 2], expfit[:, 2, 2],
    label = "Exponential fit for first 201 data points",
    line = (:dash, 2),
    color = :orange
) 
plot!(lower_area, linfit1[:, 1, 2], linfit1[:, 2, 2],
    label = "Linear fit for data points 201 to 320",
    line = (:dash, 2),
    color = :red
    )
plot!(lower_area, linfit2[:, 1, 2], linfit2[:, 2, 2],
    label = "Linear fit for data points 320 to 465",
    line = (:dash, 2),
    color = :green
    )

annotate!((t[1], 0.25, text("Exponential growth lasts \n for ≈ 2.8 years", 10, :left, :orange)))
annotate!((t[801], 0.25, text("Linear growth first phase lasts for ≈ 1.35 years. \n In this time the growth of 90% of the area is 28%. \n This gives diffusivity of 2.3e6m²/s.", 10, :right, :red)))
annotate!((t[801], 0.75, text("Linear growth second phase \n lasts for ≈ 1.65 years. \n In this time the growth \n of 90% of the area is 60%. \n This gives diffusivity of 4e6m²/s.", 10, :right, :green)))

phys_params = nondim2dim(data)
steps = t[450] / data["clock/dt"]
days = (steps * phys_params["Δt"]) / 3600
years = days / 365

linphase1 = ((t[465] / data["clock/dt"]) * phys_params["Δt"] - (t[320] / data["clock/dt"]) * phys_params["Δt"])/ (3600 * 365)

area_per[465, 2] - area_per[320, 2]
d1 = diffusivity(data, [1 2; 201 320]; conc_min = 0.1)
d2 = diffusivity(data, [1 2; 320 465]; conc_min = 0.1)

tracer_growth = plot(p1, logp1, lower_area, layout = @layout([a b; c]), size = (1200, 1200))
save("tracer_growth_90_dom32.png", tracer_growth)