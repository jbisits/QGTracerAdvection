#' # PV projection into the tracer simulations

#' ## Snapshot of PV field
#' First need a snapshot of a PV field. 
#' The PV field is from an advection-diffusion simulation using the QGPV snapshot as the initial condition.

using JLD2, Plots

SimPath = joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection", "Honours thesis/Experiment")
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = QGPV, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

data = load(file)

QGPV_init = data["snapshots/Concentration/"*string(2000)]

x = data["grid/x"]
y = data["grid/y"]
Lₓ = data["grid/Lx"]
nx = data["grid/nx"]

heatmap(x, y, QGPV_init[:, :, 1]', 
        title = "Snapshot of PV field",
        xlabel = "x", ylabel = "y", colorbar_title = "PV")


#+ echo = false
conc_path = joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true/SimulationData.jld2")
conc = load(conc_path);
#+

#' ## PV ($q$) as quasi meridional coordinate $\tilde{y}$
#' This is a shortcut to see if what I am doing is correct.
#' Take the first column of PV values (that is all meridional values at the first zonal gridpoint) and use those for the transformation.
#' This is instead of the histogram in terms of area and dividing by zonal width (apologies if this shortcut is too short!).
#' Then reshape and sort these values into an array of smallest PV to largest PV and we get the plot of $y - q(y)$,

sorted_PV = sort(QGPV_init[1, :, 1])
plot(y, sorted_PV, xlabel = "y", ylabel = "ỹ = q(y)", label = false)

#' ## Tracer concentration as a function of PV using $\tilde{y}$
#' The above gives the function $\tilde{y} = q(y)$ which can be used as the meridional cooridinate?
#' The tracer concentration with the usual meridional gridpoints $y$ is initially the Gaussian blob 

conc_init = conc["snapshots/Concentration/"*string(0)]
heatmap(x, y, conc_init[:, :, 1]', color = :deep, 
        xlabel = "x", ylabel = "y", colorbar_title = "Concentration")

#' Instead of the meridional gridpoints use the PV values that area a function of meridional gridpoints

heatmap(x, sorted_PV, conc_init[:, :, 1]', color = :deep, 
        title = "Concentration as a function of PV (or ỹ)",
        xlabel = "x", ylabel = "ỹ", colorbar_title = "Concentration", size = (600, 400))

#' The tracer concentration is a function of $x$ and $\tilde{y}$.
#' Is the right idea?
#' Then can take some values of PV say $(q_{1}, q_{2}) = (140, 160)$ and look at how the tracer evolves?
#' I think I might have missed something though as all I have done is rescale the meridional axis!