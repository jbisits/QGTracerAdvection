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

t = time_vec(data)
area_per = tracer_area_percentile(data; conc_min = 0.25)
p1 = plot(t, area_per, 
        label = ["Upper layer" "Lower layer"],
        title = "Growth of area of tracer patch in \n both layers layer",
        legend = :topleft
        )
logp1 = plot(t, log.(area_per), 
            label = ["Upper layer" "Lower layer"],
            title = "Growth of log(area of tracer patch) \n in both layers",
            legend = :bottomright
            )
plot(p1, logp1, layout = (2, 1), size = (700, 700))

expfit = exp_fit(data; conc_min = 0.05, tfitfinal = 81, tplot_length = 10)
linfit = linear_fit(data; conc_min = 0.05, tfitvals = [101 190], tplot_length = [0 0])
lower_area = plot(t, area_per[:, 2], 
                label = "Lower layer",
                title = "Growth of area of tracer patch in lower layer",
                legend = :topleft,
                lw =2,
                size = (700, 700)
                )
plot!(lower_area, expfit[:, 1, 2], expfit[:, 2, 2],
    label = "Exponential fit for first 81 data points",
    line = (:dash, 2),
    color = :orange
)  
plot!(lower_area, linfit[:, 1, 2], linfit[:, 2, 2],
    label = "Linear fit for data points 101 to 200",
    line = (:dash, 2),
    color = :red
    )

phys_params = nondim2dim(data)

steps = t[450] / data["clock/dt"]
days = (steps * phys_params["Δt̂"]) / 3600
years = days / 365

#Diffusivity calcuation
Area = phys_params["Lx̂"] * phys_params["Lŷ"]
Area_inc = area_per[190, 2] - area_per[101, 2]

no_of_seconds = (t[190] / data["clock/dt"]) * phys_params["Δt̂"] - (t[101] / data["clock/dt"]) * phys_params["Δt̂"]

diffusivity = (Area * Area_inc) / no_of_seconds