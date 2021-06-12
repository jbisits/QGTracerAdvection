#=
    Diagnostics to measure mixing.
    These are:
     - variance of concentration over the grid
     - evolution of the isopycnal second moment
=#
module MeasureMixing

export
    conc_var!,
    area_tracer_patch!,
    fit_normal!

using Distributions, GeophysicalFlows, StatsBase, LinearAlgebra

"""
    function conc_var!(concentration_variance, AD_prob)
Calculate the variance of the tracer concentration in each layer for advection-diffusion problem `prob` and store the
result at each timestep in the array concentration_variance. 
"""
function conc_var!(concentration_variance, AD_prob) 

    nlayers = AD_prob.params.nlayers
    step = AD_prob.clock.step + 1
    for i in 1:nlayers
        concentration_variance[step, i] = var(AD_prob.vars.c[:, :, i])
    end

end

"""
    function fit_normal!(σ², AD_prob)
Fit a normal distribution to the concentration of tracer at each time step to look at how σ² changes.
"""
function fit_normal!(σ², AD_prob)

    nlayers = AD_prob.params.nlayers
    step = AD_prob.clock.step + 1
    for i in 1:nlayers
        conc_data = reshape(AD_prob.vars.c[:, :, i], :, 1) 
        fit_norm = fit_mle(Normal, conc_data)
        σ = params(fit_norm)[2]
        σ²[step, i] = σ^2
    end

end
"""
    function isopycnal_second_moment(AD_prob)
Fit a histogram to the concentration data, then use cumulative sum to numerically integrate so that a picture
of normalised area ~ concentration. The axis are then swapped to concentration ~ normalised area.
A plot of this is then saved to a .jld2 file (and maybe also the fitted histogram?).
"""
function isopycnal_second_moment(AD_prob)

    nlayers, C = AD_prob.params.nlayers, AD_prob.vars.c
    hist_layer = Array{Histogram}(undef, nlayers)

    for i in 1:nlayers
        conc_data = reshape(C[:, :, i], :)
        temp_fit = fit(Histogram, conc_data)
        hist_layer[i] = normalize(temp_fit, mode = :probability)
        
    end

end

"""
    function area_tracer_patch!(area_vals, AD_prob, QG_prob, Kₛ)
Calculate the evolution of the area of the tracer patch in each layer that is advected by `AD_prob` and store the
result at each timestep in the array area_vals. 
"""
function area_tracer_patch!(second_moment_conc, AD_prob, QG_prob, Kₛ)

    α = 0.5
    nlayers = QG_prob.params.nlayers
    step = AD_prob.clock.step + 1
    t = AD_prob.clock.t

    Q = Array{Float64}(undef, nlayers)
    for i in 1:nlayers
        Q[i] = mean(AD_prob.vars.c[:, :, i]) #The total amount of tracer added.
    end

    # ux = ∂u/∂x and vy = ∂v/∂y; they are computed here in Fourier space and then transformed back to physical space
    uh, vh = QG_prob.vars.uh, QG_prob.vars.vh
    uxh = @. im * QG_prob.grid.kr * uh
    vyh = @. im * QG_prob.grid.l * vh

    ux, vy = QG_prob.vars.u, QG_prob.vars.v #Use these as dummy variables
    MultiLayerQG.invtransform!(ux, uxh, QG_prob.params)
    MultiLayerQG.invtransform!(vy, vyh, QG_prob.params)

    γ = Array{Float64}(undef, nlayers)
    for i in 1:nlayers
        ms =  mean(ux[:, :, i].^2 .+ vy[:, :, i].^2) 
        γ[i] = sqrt(ms)
    end
    
    Aₜ = @. π * (Kₛ / γ) * exp(α * γ * (t - 0.25 * γ^(-1)))

    @. second_moment_conc[step, :] = Q^2 * (2 * Aₜ)^(-1)
    
end

end #module
