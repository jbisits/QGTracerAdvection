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

#Calculate the variance of the tracer concentration.
function conc_var!(conc_var, prob) 

    nlayers = prob.params.nlayers
    step = prob.clock.step + 1
    for i in 1:nlayers
        @. conc_var[step, i] = var(prob.vars.c[:, :, i])
    end

end

#Second moment evolution of tracer patch

function area_tracer_patch(AD_prob, QG_prob, Kₛ)

    α = 0.5
    params = QG_prob.params
    t = AD_prob.clock.t
    #Q = total amount of tracer thought not sure how to define this.

    # ux = ∂u/∂x and vy = ∂v/∂y; they are computed here in Fourier space and then transformed back to physical space
    uh, vh = QG_prob.vars.uh, QG_prob.vars.vh
    uhx = im * QG_prob.grid.kr * uh
    vhy = im * QG_prob.grid.l * vh

    ux, vy = QG_prob.vars.u, QG_prob.vars.v #Use these as dummy variables
    MultiLayerQG.invtransform!(ux, uhx, params)
    MultiLayerQG.invtransform!(vy, vhy, params)

    γ = @. sqrt(ux^2 + vy^2) #Presumably want this as a single number so take the mean?

    Aₜ = @. π * (Kₛ / γ) * exp(α * γ * (t - 0.25 / γ))

    var_patch = Q^2 * (2 * Aₜ)^(-1)
end

end #module
