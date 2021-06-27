module ExperimentSetup

export 
    GaussianBlobIC,
    GaussianStripIC,
    PointSourceIC,
    CreateFile, 
    GetConcentration,
    nondim2dim

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
    CreateFile(ADProb)
Create directory and file for a given run that will be appended with 
flow characteristics and .jld2. If already exists uses directory and removes the file.
"""
function CreateFile(ADProb, IC::InitialCondition, SimPath)
    Lx, nx = ADProb.grid.Lx, ADProb.grid.nx
    filepath = joinpath(SimPath, 
                        "Output/Simulation: Domain = "*string(round(Int, Lx))*", Res = "*string(nx)*", IC = "*string(typeof(IC))
                        )
    if !isdir(filepath)
        mkdir(filepath)
    end
    filename = joinpath(filepath, "SimulationData.jld2")
    if isfile(filename)
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
"""
    function nondim2dim(prob)
Translates parameters that have been set in the non-dimensional space (as I use in my thesis) for a QG problem
to the phsyical space based off mid-latitude values in metres and seconds. The values have defaults set.
"""
function nondim2dim(prob;
                    Ω = 7.29e-5,     # Earth"s rotation
                    ϕ = π/3,         # Latitude
                    a = 6378e3,      # Earth's radius
                    g = 9.81,        # Gravity
                    H = 1500,        # Total depth (in metres)
                    ρ₁ = 1034,       # Density of top layer
                    ρ₂ = 1035,       # Density of bottom layer
                    )

    f₀ = 2*Ω*sin(ϕ)             # Coriolis computed from above values
    gprime = g*((ρ₂ - ρ₁)/ρ₂)   # Reduced gravity
    
    Ld = sqrt(gprime*H)/(f₀)    #Rossby deformation radius
    U = 0.1

    #Domain
    Lx̂, Lŷ = prob.grid.Lx, prob.grid.Ly
    Lx = Ld * Lx̂
    Ly = Ld * Lŷ

    #Parameters
    f̂₀, β̂, μ̂, ν̂ = prob.params.f₀, prob.params.β, prob.params.μ, prob.params.ν'
    f₀ = (U/Ld) * f̂₀
    β = (U/Ld^2) * β̂
    μ = (U/Ld) * μ̂
    ν = (U/Ld) * ν̂

    #Time
    Δt̂ = prob.clock.dt
    Δt = (Ld/U) * Δt̂

    return Dict("f̂₀" => f₀,
                "β̂"  => β,
                "μ̂"  => μ,
                "ν̂"  => ν,
                "Lx̂" => Lx,
                "Lŷ" => Ly,
                "Δt̂" => Δt
                )

end

end #module