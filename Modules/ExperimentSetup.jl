module ExperimentSetup

export 
    GaussianBlobIC,
    GaussianStripIC,
    PointSourceIC,
    QGPVIC,
    CreateFile, 
    GetConcentration

using GeophysicalFlows, Distributions, StatsBase, LinearAlgebra, JLD2
"""
    Abstract super type for initial conditions
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
struct PointSource{T, U, V} <: InitialCondition
    x  :: T
    y  :: T
    C₀ :: U
    ConcentrationAmount :: V
end
"""
    Struct for QGPV initial condition
"""
struct QGPV{T, U} <: InitialCondition
    q :: T
    plan_vort :: U
    C₀:: T
end
"""
    function GaussianBlobIC(μ, Σ, grid)
Create a Gaussian blob initial condition on a advection diffusion problem grid from a given mean 
μ (as vector) and covariance matrx Σ.
"""
function GaussianBlobIC(μ::Vector, Σ::Matrix, grid)

    xgrid, ygrid = gridpoints(grid)
    nx, ny = grid.nx, grid.ny
    C₀ = Array{Float64}(undef, nx, ny)
    blob = MvNormal(μ, Σ)
    blob_IC(x, y) = pdf(blob, [x,y])
   
    @. C₀ = blob_IC(xgrid, ygrid)

    return GaussianBlob(μ, Σ, C₀)
end
"""
    function GaussianStripIC(μ, σ², grid)
Create a Gaussian strip initial condition on a advection diffusion problem grid from a given mean 
μ and variance σ².
"""
function GaussianStripIC(μ::Union{Int, Float64}, σ²::Union{Int, Float64}, grid)

    x, y = grid.x, grid.y
    strip = Normal(μ, σ²)
    strip_IC(x) = pdf(strip, x)
    nx, ny = grid.nx, grid.ny
    C₀ = Array{Float64}(undef, nx, ny)
    for i in 1:nx
        C₀[i, :] = strip_IC(y)
    end

    return GaussianStrip(μ, σ², C₀)
end
"""
    function PointSourceIC(ConcentrationPoint, grid)
Create a point source initial condition at `ConcentrationPoint` on a advection diffusion problem grid.
"""
function PointSourceIC(ConcentrationPoint::Vector, ConcentrationAmount, grid)

    xconcpt = ConcentrationPoint[1]
    yconcpt = ConcentrationPoint[2]
    C₀ = Array{Float64}(undef, grid.nx, grid.ny)
    for i in 1:grid.nx
        for j in 1:grid.ny
            if [i, j] != [xconcpt, yconcpt]
                C₀[i, j] = 0
            else
                C₀[i, j] = ConcentrationAmount
            end
        end
    end

    return PointSource(xconcpt, yconcpt, C₀, ConcentrationAmount)
end
"""
    function QGPVIC(QGProb::FourierFlows.Problem)
Create an initial condition that is the full QGPV.
"""
function QGPVIC(QGProb::FourierFlows.Problem)

    q = QGProb.vars.q 
    f₀ = QGProb.params.f₀
    β = QGProb.params.β
    y = QGProb.grid.y
    plan_vort = f₀ .+ β .* y

    Q = q .+ plan_vort'

    return QGPV(q, plan_vort, Q)
end
"""
    CreateFile(ADProb)
Create directory and file for a given run that will be appended with 
flow characteristics and .jld2. If already exists uses directory and removes the file.
"""
function CreateFile(ADProb::FourierFlows.Problem, IC::InitialCondition, save_freq::Int64, SimPath::String; Ensemble = false)

    IC = string(typeof(IC))
    first = length("Main.ExperimentSetup.")
    last = findfirst('{', IC)
    IC = IC[(first + 1):(last - 1)]
    Lx, Ly = ADProb.grid.Lx, ADProb.grid.Ly
    nx = ADProb.grid.nx
    if Lx == Ly
        filepath = joinpath(SimPath, 
                            "Output/Simulation: Lx̂ = Lŷ = "*string(round(Int, Lx))*", nx = "*string(nx)*", save_freq = "*string(save_freq)*", IC = "*IC*", Ensemble = "*string(Ensemble)
                            )
    else
        filepath = joinpath(SimPath, 
                            "Output/Simulation: Lx̂ = "*string(round(Int, Lx))*", Lŷ = "*string(round(Int, Ly))*", save_freq = "*string(save_freq)*", IC = "*IC*", Ensemble = "*string(Ensemble)
        )
    end

    if !isdir(filepath)
        mkdir(filepath)
    end
    filename = joinpath(filepath, "SimulationData.jld2")
    if Ensemble == true & isfile(filename)
        filename = FourierFlows.uniquepath(filename)
    elseif isfile(filename)
        rm(filename)
    end

    return filename
end
"""
    GetConcentration(ADProb)
Function to extract the concentration data from the `ADProb` so it can be saved into .jld2 
that has been created for the simulation
"""
function GetConcentration(ADProb)
    
    Concentration = @. ADProb.vars.c
    return Concentration
end

end #module