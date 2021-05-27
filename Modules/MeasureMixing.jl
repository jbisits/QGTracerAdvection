#=
    Diagnostics to measure mixing.
    These are:
     - variance of concentration over the grid
     - evolution of the second moment of the tracer patch.
=#
module MeasureMixing

export
    conc_var!,
    area_tracer_patch!,
    fit_normal!

using Distributions, GeophysicalFlows

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
Fit a normal distribution to the concentration of tracer at each time step to look at how σ² grows.
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
        ms =  mean(ux[:, :, i].^2 .+ vy[:, :, i].^2) #Maybe calculating this at each step is not the way to go. Could instead calculate mean for the whole flow then do the calcualtions?
        γ[i] = sqrt(ms)
    end
    
    Aₜ = @. π * (Kₛ / γ) * exp(α * γ * (t - 0.25 * γ^(-1)))

    @. second_moment_conc[step, :] = Q^2 * (2 * Aₜ)^(-1)
    
end

end #module
