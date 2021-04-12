#=
    Diagnostics to measure mixing. My idea is to have all the mixing measuring diagnostics in one module.
    Eventually would like to set it up so that only one thing need to be called similar to the way MultiLayerQG
    has the argument `diags` which updates all relevant diagnostics. 
=#
module MeasureMixing

export
    conc_var!

using Distributions

#Calculate the variance of the tracer concentration.
function conc_var!(conc_var, prob) 

    nlayers = prob.params.nlayers
    step = prob.clock.step + 1
    temp = zeros(nlayers)
    for i in 1:nlayers
        temp[i] = var(prob.vars.c[:, :, i])
    end
    @. conc_var[step, :] = temp

end

#effective diffusivity (eventually)

end #module