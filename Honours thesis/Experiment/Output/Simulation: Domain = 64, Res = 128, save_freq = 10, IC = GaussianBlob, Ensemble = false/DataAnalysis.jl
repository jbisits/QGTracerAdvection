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

plot(t, area_per[:, 1],
     label = "Upper layer",
     title = "Growth of area of tracer patch in the upper layer",
     legend = :topleft
     )
scatter!([t[150]], [area_per[150, 1]],
        label = "Stage 2 -> stage 3",
        annotations = ([t[150] + 1], area_per[150, 1], Plots.text("Stage three begins after \n approximately 2.6 years", :left, :orange))
        )
scatter!([t[315]], [area_per[315, 1]],
        label = "Tracer patch ≈ size of domain",
        annotations = ([t[310]], .85, Plots.text("After approx 3.7 years \n tracer patch is size  \n of domain", :left, :green))
        )

phys_params = nondim2dim(data)

steps = t[90] / data["clock/dt"]
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


#########################################################
#Calculations that may turn into functions

#Slope of each increment may be too variable with many saves
slope = [(area_per[i + 1, 1] - area_per[i, 1])/(t[i + 1] - t[i]) for i in 1:(length(area_per[:, 1]) - 1)]

plot(t[2:end], slope) #This is potentially too variable to be of use

##### This exponential and linear fitting are both now functions that will fit curves to top and bottom layer for cpecifed time lengths
#test = linear_fit(data; tvals = [100 250])
#test = lexp_fit(data; tfinal = 100)
#Fit curves via least squares and see what happens as t increases

X = [ones(length(t[1:50])) t[1:50]]
M = log.(area_per[1:50, 1])
res50 = inv(X' * X) * (X' * M)

X = [ones(length(t[1:100])) t[1:100]]
M = log.(area_per[1:100, 1])
res100 = inv(X' * X) * (X' * M)

X = [ones(length(t[1:125])) t[1:125]]
M = log.(area_per[1:125, 1])
res125 = inv(X' * X) * (X' * M)

X = [ones(length(t[1:150])) t[1:150]]
M = log.(area_per[1:150, 1])
res150 = inv(X' * X) * (X' * M)

X = [ones(length(t[1:200])) t[1:200]]
M = log.(area_per[1:200, 1])
res200 = inv(X' * X) * (X' * M)

#Have fitted ln(A) = ln(α) + tk, so A = exp(α)*exp(tk) ⟹ A = exp(res[1]) * exp(res[2] .* t)

#These are the real data as well as some fitted exponential curves to see if there is a clear indication where it departs from exp growth.
plot(t, area_per[:, 1],
     label = "Upper layer",
     title = "Growth of area of tracer patch in the upper layer",
     legend = :bottomright,
     lw = 2
     )
plot!(t[1:120], exp(res50[1]) .* exp.(res50[2] .* t[1:120]), label = "Exp fit to first 50 data points")
plot!(t[1:130], exp(res100[1]) .* exp.(res100[2] .* t[1:130]), label = "Exp fit to first 100 data points")
plot!(t[1:140], exp(res125[1]) .* exp.(res125[2] .* t[1:140]), label = "Exp fit to first 125 data points")
plot!(t[1:160], exp(res150[1]) .* exp.(res150[2] .* t[1:160]), label = "Exp fit to first 150 data points")
plot!(t[1:180], exp(res200[1]) .* exp.(res200[2] .* t[1:180]), label = "Exp fit to first 200 data points")

#On the log scale
plot(t, log.(area_per[:, 1]),
     label = "Upper layer",
     title = "Growth of log(area of tracer patch) in the upper layer",
     legend = :bottomright,
     lw = 2
     )
plot!(t[1:120], res50[1] .+ res50[2] .* t[1:120], label = "Exp fit to first 50 data points")
plot!(t[1:130], res100[1] .+ res100[2] .* t[1:130], label = "Exp fit to first 100 data points")
plot!(t[1:140], res125[1] .+ res125[2] .* t[1:140], label = "Exp fit to first 125 data points")
plot!(t[1:160], res150[1] .+ res150[2] .* t[1:160], label = "Exp fit to first 150 data points")
plot!(t[1:180], res200[1] .+ res200[2] .* t[1:180], label = "Exp fit to first 200 data points")

#From the above looks like something in between 100 and 125 fits the first part of the data best.
#The fit to the 110 steps below looks reasonable
X = [ones(length(t[1:110])) t[1:110]]
M = log.(area_per[1:110, 1])
res110 = inv(X' * X) * (X' * M)
plot(t, area_per[:, 1],
     label = "Upper layer tracer area",
     title = "Growth of area of tracer patch in the upper layer",
     legend = :bottomright,
     lw = 2
     )
plot!(t[1:130], exp(res110[1]) .* exp.(res110[2] .* t[1:130]), 
    label = "Exp fit to first 110 data points",
    line = (:dash, 2))

#Now look at linear fit to the next part of the curve
X = [ones(length(t[100:250])) t[100:250]]
M = area_per[100:250, 1]
reslin = inv(X' * X) * (X' * M)

plot!(t[100:250], reslin[1] .+ reslin[2] .* t[100:250], 
    label = "Linear fit to data points 100 to 250", 
    line = (:dash, 2), color = :orange)
scatter!([t[240]], [area_per[240,1]], color = :orange, label = false)
scatter!([t[100]], [area_per[100,1]], color = :orange, label = false)

#Calculate the diffusivity from the increase in area over time. Not exactly sure ho to do this
phys_params = nondim2dim(data)

steps = t[100] / data["clock/dt"]
days = (steps * phys_params["Δt̂"]) / 3600
years = days / 365

area_growth = area_per[240,1] - area_per[100,1]
no_of_steps = t[240] / data["clock/dt"] - t[100] / data["clock/dt"]
no_of_seconds = (no_of_steps * phys_params["Δt̂"])
days = (no_of_seconds) / 3600
years = days / 365