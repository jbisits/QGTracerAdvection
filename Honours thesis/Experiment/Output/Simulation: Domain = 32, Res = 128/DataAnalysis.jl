#Load in the data

file = joinpath(SimPath, "Output/Simulation: Domain = 32, res = 128/SimulationData.jld2")
data = load(file)

#Produce histogram plots from the saved concentration data
histplots = hist_plot(data, nsteps)

#Produce heatmaps of tacer concentration from the saved concentration data
tracerplots = tracer_plot(data, ADProb, nsteps)

#Plot heatmaps and histograms togehter.
uppertacer = plot(tracerplots[1][1], tracerplots[1][2], tracerplots[1][3], tracerplots[1][4],
                    tracerplots[1][5],tracerplots[1][6])

upperhist = plot(histplots[1][1], histplots[1][2], histplots[1][3], histplots[1][4], histplots[1][5], histplots[1][6])
plot(uppertacer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[2][1], tracerplots[2][2], tracerplots[2][3], tracerplots[2][4],
                    tracerplots[2][5],tracerplots[2][6])
                    
lowerhist = plot(histplots[2][1], histplots[2][2], histplots[2][3], histplots[2][4], histplots[2][5], histplots[2][6])
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

ConcentrationMean = conc_mean(data, nsteps)
ConcentrationVaricance = conc_var(data, nsteps)
t = 0:1:nsteps
meanplot = plot(t, ConcentrationMean, label = ["Upper Layer" "Lower layer"],
                    title = "Mean concentration \n over the grid",
                    ylims = (0, 0.01)
            )
varplot = plot(t, ConcentrationVaricance, label = ["Upper Layer" "Lower layer"],
                title = "Variance of concentration \n over the grid"    
            )
plot(meanplot, varplot, sixe = (1000, 600))

plot(t, 1 ./ ConcentrationVaricance, label = ["Upper Layer" "Lower layer"],
title = "Inverse of variance of concentration over the grid"    
    )

ConcVsArea = concarea_animate(data, nsteps)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)

testfit = fit(InverseGaussian, reshape(data["snapshots/Concentration/0"], :))