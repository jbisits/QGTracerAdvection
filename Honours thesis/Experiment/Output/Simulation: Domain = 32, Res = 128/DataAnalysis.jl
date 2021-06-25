#Load in the data
data = load(filename)

#Produce histogram plots from the histogram objects that are saved
histplots = hist_plot(data)

uppertracer = plot(tracerplots[1][1], tracerplots[1][2], tracerplots[1][3], tracerplots[1][4],
                    tracerplots[1][5],tracerplots[1][6])
upperhist = plot(histplots[1][1], histplots[1][2], histplots[1][3], histplots[1][4], histplots[1][5], histplots[1][6])
plot(uppertracer, upperhist, layout=(2, 1), size = (1200, 1200))

lowertracer = plot(tracerplots[2][1], tracerplots[2][2], tracerplots[2][3], tracerplots[2][4],
                    tracerplots[2][5],tracerplots[2][6])
lowerhist = plot(histplots[2][1], histplots[2][2], histplots[2][3], histplots[2][4], histplots[2][5], histplots[2][6])
plot(lowertracer, lowerhist, layout=(2, 1), size = (1200, 1200))

ConcVsArea = concarea_animate(data, nsteps)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)