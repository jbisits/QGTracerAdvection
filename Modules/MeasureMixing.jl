#=
    Diagnostics to measure mixing and functions for data extraction from a .jld2 file created by a tracer 
    advection diffusion simulation.
    The diagnostics are:
     - variance of concentration over the grid
     - evolution of the isopycnal second moment
    From the .jld2 file can create plots and extract relevant information.
=#
module MeasureMixing

export
    conc_var!,
    area_tracer_patch!,
    fit_normal!, 
    fit_hist!,
    hist_plot,
    concarea_plot,
    concarea_animate

using Distributions, GeophysicalFlows, StatsBase, LinearAlgebra, JLD2, Plots

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
    function fit_hist!(filename, AD_prob, max_conc)
Fits a histogram to the concentration data at each time step. From the histogram the concentration data 
and area data can be extracted.
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
"""
    function hist_plot(data)
Create plots of histograms at the same timesteps as the tracer plots from the saved data
in the output file. The input `data` is the loaded .jld2 file.
"""
function hist_plot(data)
    step_nums = data["SaveStepsforTracerPlots"]
    max_conc = data["MaxConcentration"]
    UpperConcentrationHistograms = Plots.Plot{Plots.GRBackend}[]
    LowerConcentrationHistograms = Plots.Plot{Plots.GRBackend}[]
    for i ∈ step_nums
        push!(UpperConcentrationHistograms, plot(data["Histograms/step"*string(i)][1],
                                                    label = false, 
                                                    xlabel = "Concentration", 
                                                    ylabel = "Normalised area",
                                                    xlims = (0, max_conc[1] + 0.01)))
        push!(LowerConcentrationHistograms, plot(data["Histograms/step"*string(i)][2],
                                                    label = false, 
                                                    xlabel = "Concentration", 
                                                    ylabel = "Normalised area",
                                                    xlims = (0, max_conc[2] + 0.01)))
    end
    return [UpperConcentrationHistograms, LowerConcentrationHistograms]
end
"""
    function concarea_plot(data)
Create plots of Concetration ~ normalised area at the same time steps as the tracer plots from the 
saved data in the output file. The input `data` is the loaded .jld2 file.
"""
function concarea_plot(data)
    step_nums = data["SaveStepsforTracerPlots"]
    max_conc = data["MaxConcentration"]
    UpperConcentrationArea = Plots.Plot{Plots.GRBackend}[]
    LowerConcentrationArea = Plots.Plot{Plots.GRBackend}[]
    for i ∈ step_nums
        push!(UpperConcentrationArea, plot(data["ConcentrationData/step"*string(i)][1], data["Histograms/step"*string(i)][1].edges,
                label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                xlims = (0, max_conc[1] + 0.01)
                ))
        push!(LowerConcentrationArea, plot(data["ConcentrationData/step"*string(i)][2], data["Histograms/step"*string(i)][2].edges,
                label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                xlims = (0, max_conc[2] + 0.01)
                ))
    end
    return [UpperConcentrationArea, LowerConcentrationArea]
end
"""
    function concarea_animate(data, nsteps)
Create an animation of Concetration ~ normalised area from the saved data in the output file.
"""
function concarea_animate(data, nsteps)

    max_conc = data["MaxConcentration"]
    ConcVsArea = @animate for i in 0:10:nsteps
    p1 = plot(data["ConcentrationData/step"*string(i)][1], data["Histograms/step"*string(i)][1].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 ylims = (0, max_conc[1] + 0.01),
                 title = "Top layer"
                )
    p2 = plot(data["ConcentrationData/step"*string(i)][2], data["Histograms/step"*string(i)][2].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 ylims = (0, max_conc[2] + 0.01),
                 title = "Bottom layer"
                )
    plot(p1, p2)
    end
    return ConcVsArea
end

end #module
