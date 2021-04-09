#=
    Diagnostics to measure mixing. My idea is to have all the mixing measuring diagnostics in one module.
    Eventually would like to set it up so that only one thing need to be called similar to the way MultiLayerQG
    has the argument `diags` which updates all relevant diagnostics. 
=#
module MeasureMixing

export
    var_vector!

using Distributions
#=
    Calculate the variance of tracer distribution and save to vector
=#
function var_vector!(var_vals, prob) 

    nlayers = prob.params.nlayers
    step = prob.clock.step + 1
    temp = zeros(nlayers)
    for i in 1:nlayers
        temp[i] = var(prob.vars.c[:, :, i])
    end
    @. var_vals[step, :] = temp

end

#effective diffusivity (eventually)

end #module