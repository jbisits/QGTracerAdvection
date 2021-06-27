#Change to the current directory
cd(joinpath(SimPath, "Output/Simulation: Domain = 32, res = 128, save_freq = 200"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

#Produce histogram plots from the saved concentration data
histplots = hist_plot(data)

#Produce heatmaps of tacer concentration from the saved concentration data
tracerplots = tracer_plot(data)

#Plot heatmaps and histograms togehter.
uppertacer = plot(tracerplots[1][1], tracerplots[1][2], tracerplots[1][3], tracerplots[1][4],
                    tracerplots[1][5],tracerplots[1][6])

upperhist = plot(histplots[1][1], histplots[1][2], histplots[1][3], histplots[1][4], histplots[1][5], histplots[1][6])
plot(uppertacer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[2][1], tracerplots[2][2], tracerplots[2][3], tracerplots[2][4],
                    tracerplots[2][5],tracerplots[2][6])

lowerhist = plot(histplots[2][1], histplots[2][2], histplots[2][3], histplots[2][4], histplots[2][5], histplots[2][6])
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

#Time vector for plotting
t = time_vec(data)

#Create some plots of concentration diagnostics.
ConcentrationMean = conc_mean(data)
ConcentrationVaricance = conc_var(data)

meanplot = plot(t, ConcentrationMean, 
                    label = ["Upper Layer" "Lower layer"],
                    title = "Mean concentration \n over the grid",
                    xlabel = "t",
                    ylabel = "Concentration",
                    ylims = (0, 0.01)
                )
varplot = plot(t, ConcentrationVaricance, 
                    label = ["Upper Layer" "Lower layer"],
                    title = "Variance of concentration \n over the grid",
                    xlabel = "t",
                    ylabel = "Concentration"    
                )
plot(meanplot, varplot, size = (1000, 600))

uppersecondmoment = plot(t, 1 ./ ConcentrationVaricance[:, 1], 
                            label = "Upper layer",
                            title = "Inverse of variance of concentration \n over the upper layer grid",
                            xlabel = "t",
                            legend = :topleft 
                        )
lowersecondmoment = plot(t, 1 ./ ConcentrationVaricance[:, 2], 
                            label = "Lower layer",
                            title = "Inverse of variance of concentration \n over the lower layer grid",
                            xlabel = "t",
                            legend = :topleft    
                        )
plot(uppersecondmoment, lowersecondmoment, size = (1000, 600))

#Instead consider the integral ∫C²dA which is the concentration per unit area as Garrett defines
conc_int = Garrett_int(data)

plot(t, conc_int, 
        label = ["Upper layer" "Lower layer"], 
        xlabel = "t", 
        ylabel = "∫C²dA",
        title = "Concentration per unit area \n calculated at each time step"
    )
plot(t, 1 ./ conc_int, 
        label = ["Upper layer" "Lower layer"],
        xlabel = "t", 
        ylabel = "(∫C²dA)⁻¹",
        title = "Inverse concentration per unit area \n calculated at each time step", 
        legend = :topleft
    )

ConcVsArea = concarea_animate(data)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)