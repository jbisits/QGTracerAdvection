#=
    Diagnostics to measure mixing. My idea is to have all the mixing measuring diagnostics in one module.
    Eventually would like to set it up so that only one thing need to be called similar to the way MultiLayerQG
    has the argument `diags` which updates all relevant diagnostics. 
=#
module MeasureMixing

export
    second_moment!

using Distributions
#=
    Calculate the second moment of the tracer distribution and save to vector
=#
function second_moment!(second_moment, prob) 

    nlayers = prob.params.nlayers
    step = prob.clock.step + 1
    temp = zeros(nlayers)
    for i in 1:nlayers
        #Use variance or second moment here? Ask Jan and Geoff
        temp[i] = sum(prob.vars.c[:, :, i].^2)
        #temp[i] = var(prob.vars.c[:, :, i])
    end
    @. second_moment[step, :] = temp

end

#effective diffusivity (eventually)

end #module