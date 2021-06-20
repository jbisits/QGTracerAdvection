#Tracer advcetion diffusion on large domain

#Experiment setup script
include("ExperimentSetup.jl")

#Import a flow that has already been set up.
include("Flows/FlowSetup_500domain_res128.jl")
#include("Flows/FlowSetup_1000domain_res128.jl")
#include("Flows/FlowSetup_1000domain_res256.jl")
#include("Flows/FlowSetup_1500domain_res128.jl")
#include("Flows/FlowSetup_1500domain_res256.jl")

κ = 0.01
#Set delay time (that is flow for some number of days, then drop tracer in)
delay_time = 0
#Set the tracer advection probelm by passing in the QG problem 
ADProb = TracerAdvDiff_QG.Problem(;prob = QGProb, delay_time = delay_time, nsubs = nsubs, κ = κ)
ADSol, ADClock, ADVars, ADParams, ADGrid = ADProb.sol, ADProb.clock, ADProb.vars, ADProb.params, ADProb.grid
ADx, ADy = gridpoints(ADGrid)
x, y = ADGrid.x, ADGrid.y

#Set the Gaussian blob initial condition
μIC = [0, 0]
Σ = [1 0; 0 1]
IC = GaussianBlobIC(μIC, Σ, ADGrid)

TracerAdvDiff_QG.QGset_c!(ADProb, IC.C₀)

max_conc = [findmax(ADVars.c[:, :, 1])[1], findmax(ADVars.c[:, :, 2])[1]]

#Define frequency at which to save a plot.
#plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
#Setting them the same gives plots at equal time increments. (Might be a better work around)
plot_time_AD, plot_time_inc = 0.2, 0.2
#Blank array to save the step number so can plot histogram corresponding to the tracer advection plot.
step_nums = []

#Create a file to save data to
filename = CreateFile(ADProb)

#Set arguments for the plots from the ADProb
plotargs = Set_plotargs(ADProb)
#Step the tracer advection problem forward and plot at the desired time step.
while ADClock.step <= nsteps
    if ADClock.step == 0
        tp_u = heatmap(x, y, ADVars.c[:, :, 1]'; plotargs...)
        push!(UpperLayerTracerPlots, tp_u)
        tp_l = heatmap(x, y, ADVars.c[:, :, 2]'; plotargs...)
        push!(LowerLayerTracerPlots, tp_l)
        push!(step_nums, ADClock.step)
    elseif round(Int64, ADClock.step) == round(Int64, plot_time_AD*nsteps)
        tp_u = heatmap(x, y, ADVars.c[:, :, 1]'; plotargs)
        push!(UpperLayerTracerPlots, tp_u)
        tp_l = heatmap(x, y, ADVars.c[:, :, 2]'; plotargs)
        push!(LowerLayerTracerPlots, tp_l)
        push!(step_nums, ADClock.step)
        global plot_time_AD += plot_time_inc
    end
    MeasureMixing.fit_hist!(filename, ADProb, number_of_bins = 30)
    stepforward!(ADProb, nsubs)
    TracerAdvDiff_QG.QGupdatevars!(ADProb)
    TracerAdvDiff_QG.vel_field_update!(ADProb, QGProb, nsubs)
end

UpperTracerPlots = plot(UpperLayerTracerPlots[1], UpperLayerTracerPlots[2], 
                        UpperLayerTracerPlots[3], UpperLayerTracerPlots[4],
                        UpperLayerTracerPlots[5], UpperLayerTracerPlots[6])
     

LowerTracerPlots = plot(LowerLayerTracerPlots[1], LowerLayerTracerPlots[2], 
                    LowerLayerTracerPlots[3], LowerLayerTracerPlots[4],
                    LowerLayerTracerPlots[5], LowerLayerTracerPlots[6])

#Now can load the output of the .jld2 file created to create histograms and plots of Concentration ~ normalised area.

data = load(filename)

histargs1 = (
    label = false, 
    xlabel = "Concentration", 
    ylabel = "Normalised area",
    xlims = (0, max_conc[1] + 0.01)
)
histargs2 = (
    label = false, 
    xlabel = "Concentration", 
    ylabel = "Normalised area",
    xlims = (0, max_conc[2] + 0.01)
)

for i ∈ step_nums
    push!(UpperConcentrationHistograma, plot(data["Histograms/step"*string(i)][1]; histargs1))
    push!(LowerConcentrationHistograma, plot(data["Histograms/step"*string(i)][2]; histargs2))
end

UpperConcentrationHistograms = plot(UpperConcentrationHistograms[1], UpperConcentrationHistograms[2], 
                                    UpperConcentrationHistograms[3], UpperConcentrationHistograms[4], 
                                    UpperConcentrationHistograms[5], UpperConcentrationHistograms[6])
LowerConcentrationHistograms = plot(LowerConcentrationHistograms[1], LowerConcentrationHistograms[2], 
                                    LowerConcentrationHistograms[3], LowerConcentrationHistograms[4], 
                                    LowerConcentrationHistograms[5], LowerConcentrationHistograms[6])

plot(UpperTracerPlots, UpperConcentrationHistograms, layout=(2, 1), size=(1200, 1200))
plot(LowerTracerPlots, LowerConcentrationHistograms, layout=(2, 1), size=(1200, 1200))

#Concentration ~ area example plot
plot(data["ConcentrationData/step400"][1], data["Histograms/step400"][1].edges,
         label = false,
        xlabel = "Normalised area",
        ylabel = "Concentration",
        xlims = (0, max_conc[1] + 0.01)
        )

#Animation of #Concentration ~ area 
#This animation is already saved so is commented out
#=
ConcVsArea = @animate for i in 1:10:nsteps
    p1 = plot(data["ConcentrationData/step"*string(i)][1], data["Histograms/step"*string(i)][1].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 xlims = (0, max_conc[1] + 0.01),
                 title = "Top layer"
                )
    p2 = plot(data["ConcentrationData/step"*string(i)][2], data["Histograms/step"*string(i)][2].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 xlims = (0, max_conc[2] + 0.01),
                 title = "Bottom layer"
                )
    plot(p1, p2)
end

mp4(ConcVsArea, "Movies/ConcVsArea.mp4", fps=18)
=#