#Import the required packages to use for the simulation

include(joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection/Modules/TracerAdvDiff_QG.jl"))
include(joinpath("/Users/Joey/Desktop/ThesisCode/QG_tracer_advection/Modules/MeasureMixing.jl"))
using .TracerAdvDiff_QG, .MeasureMixing #local modules
using GeophysicalFlows.MultiLayerQG, Plots, Distributions, StatsBase, LinearAlgebra, JLD2 #packages
using Random: seed!

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

struct PointSource{T, U} <: InitialCondition
    x  :: T
    y  :: T
    C₀ :: U
end

struct Square{T, U} <: InitialCondition
    x  :: T
    y  :: T
    C₀ :: U
end

function GaussianBlobIC(μ, Σ, grid)
    x, y = grid.x, grid.y
    xgrid, ygrid = gridpoints(grid)
    blob = MvNormal(μ, Σ)
    blob_IC(x, y) = pdf(blob, [x, y])
    C₀ = @. blob_IC(xgrid, ygrid)
    return GaussianBlob(μ, Σ, C₀)
end

function GaussianStripIC(μ, σ², grid)
    x = grid.x
    xgrid, ygrid = gridpoints(grid)
    strip = Normal(μ, σ²)
    strip_IC(x) = pdf(strip, x)
    C₀ = Array{Float64}(undef, grid.nx, grid.ny)
    for i in 1:grid.nx
        C₀[i, :] = strip_IC(ygrid[i, :])
    end
    return GaussianStrip(μ, σ², C₀)
end

function PointSourceIC(ConcentrationPoint, grid)
    xconcpt = ConcentrationPoint[1]
    yconcpt = ConcentrationPoint[2]
    C₀ = Array{Float64}(undef, grid.nx, grid.ny)
    for i in 1:grid.nx
        for j in 1:grid.ny
            if [i, j] != [xconcpt, yconcpt]
                C₀[i, j] = 0
            else
                C₀[i, j] = 1
            end
        end
    end
    return PointSource(xconcpt, yconcpt, C₀)
end

function SquareIC(x, grid)
    C₀ = zeros(grid.nx, grid.ny)
    for i ∈ x
        for j ∈ x
            C₀[i, j] = 1
        end
    end
    return Square(x, x, C₀)
end

#Create directory and file for a given run that will be appended with flow characteristics and .jld2. If already exists uses directory and removes the file.
function CreateFile(ADProb)
    Lx, nx = ADProb.grid.Lx, ADProb.grid.nx
    filepath = joinpath(pwd(), "Honours thesis/Experiment/Output/Simulation: Domain = "*string(round(Int, Lx))*", Res = "*string(nx))
    if !isdir(filepath)
        mkdir(filepath)
    end
    filename = joinpath(filepath, "SimulationData.jld2")
    if isfile(filename)
        rm(filename)
    end
    return filename
end

#Blank arrays in which to store the plots of tracer diffusion in each layer.
UpperLayerTracerPlots = Plots.Plot{Plots.GRBackend}[]
LowerLayerTracerPlots = Plots.Plot{Plots.GRBackend}[]

#Key arguments for plots
function Set_plotargs(ADProb)
    ADGrid = ADProb.grid
    plotargs = (
                aspectratio = 1,
                c = :deep,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-ADGrid.Lx/2, ADGrid.Lx/2),
                ylims = (-ADGrid.Ly/2, ADGrid.Ly/2),
                #=
                These arguments are for plots on a larger domain.
                xticks = (-ADGrid.Lx/2:round(Int, ADGrid.Lx/6):ADGrid.Lx/2, string.(-Int(ADGrid.Lx/2e3):round(Int, ADGrid.Lx/6e3):Int(ADGrid.Lx/2e3))),
                yticks = (-ADGrid.Lx/2:round(Int, ADGrid.Lx/6):ADGrid.Lx/2, string.(-Int(ADGrid.Lx/2e3):round(Int, ADGrid.Lx/6e3):Int(ADGrid.Lx/2e3)))
                =#
    )
    return plotargs
end