#=
    Diagnostics to measure mixing.
    These are:
     - variance of concentration over the grid
     - evolution of the second moment of the tracer patch.
    
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

#Second moment evolution of tracer patch

end #module