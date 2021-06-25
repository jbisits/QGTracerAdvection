module ExperimentSetup

export 
    GaussianBlobIC,
    GaussianStripIC,
    PointSourceIC,
    CreateFile

using GeophysicalFlows, Distributions, StatsBase, LinearAlgebra, JLD2
"""
Abstract uper type for initial conditions
"""
abstract type InitialCondition end
"""
    Struct for Gaussian blob initial condition
"""
struct GaussianBlob{T, U, V} <: InitialCondition
    μ  :: T
    Σ  :: U
    C₀ :: V
end
"""
    Struct for Gaussian strip initial condition
"""
struct GaussianStrip{T, U} <: InitialCondition
    μ  :: T
    σ² :: T
    C₀ :: U
end
"""
    Strcut for point source initial condition
"""
struct PointSource{T, U} <: InitialCondition
    x  :: T
    y  :: T
    C₀ :: U
end
"""
    function GaussianBlobIC(μ, Σ, grid)
Create a Gaussian blob initial condition on a advection diffusion problem grid from a given mean 
μ (as vector) and covariance matrx Σ.
"""
function GaussianBlobIC(μ::Vector, Σ::Matrix, grid)
    x, y = grid.x, grid.y
    xgrid, ygrid = gridpoints(grid)
    blob = MvNormal(μ, Σ)
    blob_IC(x, y) = pdf(blob, [x, y])
    C₀ = @. blob_IC(xgrid, ygrid)
    return GaussianBlob(μ, Σ, C₀)
end
"""
    function GaussianStripIC(μ, σ², grid)
Create a Gaussian strip initial condition on a advection diffusion problem grid from a given mean 
μ and variance σ².
"""
function GaussianStripIC(μ::Union{Int, Float64}, σ²::Union{Int, Float64}, grid)
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
"""
    function PointSourceIC(ConcentrationPoint, grid)
Create a point source initial condition at `ConcentrationPoint` on a advection diffusion problem grid.
"""
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

"""
    Struct to hold output from a tracer advection diffusion simulation.
"""
struct AdvectionOutput
    ADProb :: FourierFlows.Problem
    filepath :: String
    fields :: Dict{Symbol, Function}
end
"""
    CreateFile(ADProb)
Create directory and file for a given run that will be appended with 
flow characteristics and .jld2. If already exists uses directory and removes the file.
"""
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

end #module