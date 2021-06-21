#Load in the data
data = load(filename)

#Produce histogram plots from the histogram objects that are saved
histplots = MeasureMixing.hist_plot(data)

#Extract tracer plots from the simulation
uppertracerplots = data["TracerPlots"]

uppertracer = plot(uppertracerplots[1][1], uppertracerplots[1][2], uppertracerplots[1][3], uppertracerplots[1][4],
                    uppertracerplots[1][5],uppertracerplots[1][6])
upperhist = plot(histplots[1][1], histplots[1][2], histplots[1][3], histplots[1][4], histplots[1][5], histplots[1][6])
plot(uppertracer, upperhist, layout=(2, 1), size = (1200, 1200))

ConcVsArea = MeasureMixing.concarea_animate(data, nsteps)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)