#=
    Diagnostics to measure mixing.
    These are:
     - variance of concentration over the grid
     - evolution of the second moment of the tracer patch.
=#
module MeasureMixing

export
    conc_var!
    area_tracer_patch

using Distributions, GeophysicalFlows

"""
    function conc_var!(concentration_variance, prob)
Calculate the variance of the tracer concentration for advection-diffusion problem `prob`.
"""
function conc_var!(concentration_variance, prob) 

    nlayers = prob.params.nlayers
    step = prob.clock.step + 1
    for i in 1:nlayers
        @. concentration_variance[step, i] = var(prob.vars.c[:, :, i])
    end

end

"""
    function area_tracer_patch(AD_prob, QG_prob, Kₛ)
Calculate the evolution of the area of the tracer patch that is advected by `AD_prob`.
"""
function area_tracer_patch(AD_prob, QG_prob, Kₛ)

    α = 0.5
    params = QG_prob.params
    t = AD_prob.clock.t

    # ux = ∂u/∂x and vy = ∂v/∂y; they are computed here in Fourier space and then transformed back to physical space
    uh, vh = QG_prob.vars.uh, QG_prob.vars.vh
    uxh = @. im * QG_prob.grid.kr * uh
    vyh = @. im * QG_prob.grid.l * vh

    ux, vy = QG_prob.vars.u, QG_prob.vars.v #Use these as dummy variables
    MultiLayerQG.invtransform!(ux, uxh, params)
    MultiLayerQG.invtransform!(vy, vyh, params)

    γ = mean(@. sqrt(ux^2 + vy^2))

    Aₜ = @. π * (Kₛ / γ) * exp(α * γ * (t - 0.25 / γ))

    var_patch = Q^2 * (2 * Aₜ)^(-1)
end

end #module
