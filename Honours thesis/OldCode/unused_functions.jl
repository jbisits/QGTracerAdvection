# These are functions that were initially in the MeasureMixing.jl module but 
# are now not used. This is either because they apply to passing in blank arrays
# during the simulation which I now do not do (it is better to save the data)
# or they are not useful for my thesis.

"""
    function conc_var!(concentration_variance, AD_prob)
Calculate the variance of the tracer concentration in each layer for advection-diffusion problem `prob` and store the
result at each timestep where the data was saved in the array concentration_variance. 
"""
function conc_var!(concentration_variance::Array, AD_prob::FourierFlows.Problem) 

    nlayers = AD_prob.params.nlayers
    step = AD_prob.clock.step + 1
    for i in 1:nlayers
        concentration_variance[step, i] = var(AD_prob.vars.c[:, :, i])
    end

end
"""
    function Garrett_int(data::Dict{String, Any})
Compute the diagnostic for tracer concentration ∫C²dA at each saved data timestep (Garrett 1983).
"""
function Garrett_int(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    save_freq = data["save_freq"]
    saved_steps = length(0:save_freq:nsteps)
    conc_int = Array{Float64}(undef, saved_steps, nlayers)
    for i in 1:saved_steps
        conc_int[i, :] = [sum(data["snapshots/Concentration/"*string( (i-1) * save_freq )][:, :, j].^2) for j ∈ 1:nlayers]
    end

    return conc_int
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
    function fit_hist!(filename, AD_prob, max_conc)
Fits a histogram to the concentration data at each time step. From the histogram the concentration data 
and area data can be extracted. This is for use in a simulation like in `QG_hist.jl`.
"""
function fit_hist!(filename, AD_prob; number_of_bins = 0)

    nlayers, C = AD_prob.params.nlayers, AD_prob.vars.c
    hist_layer = Array{Histogram}(undef, nlayers)
    conc_data = []

    for i in 1:nlayers
        if number_of_bins == 0
            temp_fit = fit(Histogram, reshape(C[:, :, i], :))
        else
            temp_fit = fit(Histogram, reshape(C[:, :, i], :), nbins = number_of_bins)
        end
        hist_layer[i] = normalize(temp_fit, mode = :probability)
        temp_conc_data = cumsum(reverse(hist_layer[i].weights))
        push!(conc_data, reverse!(vcat(0, temp_conc_data)))
    end
    jldopen(filename, "a+") do file
        file["Histograms/step"*string(AD_prob.clock.step)] = hist_layer
        file["ConcentrationData/step"*string(AD_prob.clock.step)] = conc_data
    end

end