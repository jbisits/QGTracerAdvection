#Load in the data
data = load(filename)

#Produce histogram plots from the histogram objects that are saved
histplots = MeasureMixing.hist_plot(data)

#Extract tracer plots from the simulation
uppertracerplots = data["TracerPlots"]

uppertracer = plot(uppertracerplots[1]...)
upperhist = plot(histplots[1]...)
plot(uppertracer, upperhist, layout=(2, 1), size = (1200, 1200))

ConcVsArea = MeasureMixing.concarea_animate(data, nsteps)
mp4(ConcVsArea, "ConcVsArea.mp4", fps=18)