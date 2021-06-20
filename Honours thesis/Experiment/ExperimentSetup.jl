#Import the required packages to use for the simulation
using .TracerAdvDiff_QG
using GeophysicalFlows.MultiLayerQG, Plots, Distributions, StatsBase, LinearAlgebra, JLD2

#Initial conditions
abstract type InitialCondition end

struct GaussianBlob{T, U, V} <: InitialCondition
    μ  :: T
    Σ  :: U
    C₀ :: V
end

struct GaussianStrip{T, U} <: InitialCondition
    μ  :: T
    σ² :: T
    C₀ :: U
end

function GaussianBlobIC(μ, Σ, grid)
    x, y = gridpoints(grid)
    blob = MvNormal(μ, Σ)
    blob_IC(x, y) = pdf(blob, [x, y])
    C₀ = @. blob_IC(ADx, ADy)
    return GaussianBlob(μ, Σ, C₀)
end

function GaussianStripIC(μ, σ², grid)
    x, y = gridpoints(grid)
    strip = Normal(μ, σ²)
    strip_IC(x) = pdf(strip, x)
    C₀ = Array{Float64}(undef, grid.nx, grid.ny)
    for i in 1:grid.nx
        C₀[i, :] = strip_IC(y[i, :])
    end
    return GaussianStrip(μ, σ², C₀)
end

#Create file name for a given run that will be appended with flow characteristics and .jld2
function CreateFile(ADProb)
    Lx, nx = ADProb.grid.Lx, ADProb.grid.nx
    filepath = pwd()
    filename = joinpath(filepath, "Honours thesis/Experiment/Output/SavedData: domain = "*string(round(Int, Lx/1e3))*", res = "*string(nx)*".jld2")
    return filename
end

#Blank arrays in which to store the plots of tracer diffusion in each layer.
UpperLayerTracerPlots = Plots.Plot{Plots.GRBackend}[]
LowerLayerTracerPlots = Plots.Plot{Plots.GRBackend}[]

#Key arguments for plots
plotargs = (
            aspectratio = 1,
            c = :deep,
            xlabel = "x",
            ylabel = "y",
            colorbar = true,
            xlims = (-ADGrid.Lx/2, ADGrid.Lx/2),
            ylims = (-ADGrid.Ly/2, ADGrid.Ly/2),
            xticks = (-ADGrid.Lx/2:round(Int, ADGrid.Lx/6):ADGrid.Lx/2, string.(-Int(ADGrid.Lx/2e3):round(Int, ADGrid.Lx/6e3):Int(ADGrid.Lx/2e3))),
            yticks = (-ADGrid.Lx/2:round(Int, ADGrid.Lx/6):ADGrid.Lx/2, string.(-Int(ADGrid.Lx/2e3):round(Int, ADGrid.Lx/6e3):Int(ADGrid.Lx/2e3))),
            title = "C(x,y,t), day = "*string(round(Int, ADClock.t/3600))
)

#Histograms for concentration ~ normalised area
UpperConcentrationHistograms = Plots.Plot{Plots.GRBackend}[]
LowerConcentrationHistograms = Plots.Plot{Plots.GRBackend}[]

histargs = (
    label = false, 
    xlabel = "Concentration", 
    ylabel = "Normalised area",
    xlims = (0, max_conc[2] + 0.01)
)