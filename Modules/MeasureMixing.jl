#=
    Diagnostics to measure mixing and functions for data extraction from a .jld2 file created by a tracer 
    advection diffusion simulation.
=#
module MeasureMixing

export
    conc_mean,
    conc_var!,
    conc_var,
    Garrett_int,
    area_tracer_patch!,
    fit_normal!, 
    fit_hist!,
    hist_plot,
    concarea_plot,
    concarea_animate,
    tracer_plot,
    tracer_animate,
    time_vec,
    tracer_area_avg,
    tracer_area_percentile,
    avg_ensemble_tracer_area,
    exp_fit,
    linear_fit,
    diffusivity,
    nondim2dim

using Distributions, GeophysicalFlows, StatsBase, LinearAlgebra, JLD2, Plots
"""
    conc_mean(data::Dict{String, Any})
Calculate the mean concentration at each timestep where the data was saved.
"""
function conc_mean(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    save_freq = data["save_freq"]
    saved_steps = length(0:save_freq:nsteps)
    concentration_mean = Array{Float64}(undef, saved_steps, nlayers)
    for i ∈ 1:saved_steps
        concentration_mean[i, :] = [mean(data["snapshots/Concentration/"*string( (i-1) * save_freq )][:, :, j]) for j ∈ 1:nlayers]
    end

    return concentration_mean
end
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
    function conc_var(data::Dict{String, Any})
Compute the concentration variance from saved output for an advection-diffusion problem.
"""
function conc_var(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    save_freq = data["save_freq"]
    saved_steps = length(0:save_freq:nsteps)
    concentration_variance = Array{Float64}(undef, saved_steps, nlayers)
    for i ∈ 1:saved_steps
        concentration_variance[i, :] = [var(data["snapshots/Concentration/"*string( (i-1) * save_freq )][:, :, j]) for j ∈ 1:nlayers]
    end

    return concentration_variance
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
"""
    function hist_plot(data)
Create plots of histograms at the same timesteps as the tracer plots from the saved data
in the output file. The input `data` is the loaded .jld2 file.
"""
function hist_plot(data::Dict{String, Any}; plot_freq = 1000, number_of_bins = 0, xlims_same = false)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    plot_steps = 0:plot_freq:nsteps
    max_conc = [findmax(data["snapshots/Concentration/0"][:, :, i])[1] for i ∈ 1:nlayers]
    UpperConcentrationHistograms = Plots.Plot{Plots.GRBackend}[]
    LowerConcentrationHistograms = Plots.Plot{Plots.GRBackend}[]
    if xlims_same == true
        upperxlims = (0, max_conc[1])
        lowerxlims = (0, max_conc[2])
    else
        upperxlims = nothing
        lowerxlims = nothing
    end
    for i ∈ plot_steps
        upperdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 1], :)
        lowerdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 2], :)
        if number_of_bins == 0
            upperhist = fit(Histogram, upperdata)
            lowerhist = fit(Histogram, lowerdata)
        else
            upperhist = fit(Histogram, upperdata, nbins = number_of_bins)
            lowerhist = fit(Histogram, lowerdata, nbins = number_of_bins)
        end
        upperhist = normalize(upperhist, mode = :probability)
        lowerhist = normalize(lowerhist, mode = :probability)
        push!(UpperConcentrationHistograms, plot(upperhist,
                                                    label = false, 
                                                    xlabel = "Concentration", 
                                                    ylabel = "Normalised area",
                                                    xlims = upperxlims
                                                )
            )
        push!(LowerConcentrationHistograms, plot(lowerhist,
                                                    label = false, 
                                                    xlabel = "Concentration", 
                                                    ylabel = "Normalised area",
                                                    xlims = lowerxlims
                                                )
            )
    end

    return [UpperConcentrationHistograms, LowerConcentrationHistograms]
end
"""
    function concarea_plot(data)
Create plots of Concetration ~ normalised area at the same time steps as the tracer plots from the 
saved data in the output file. The input `data` is the loaded .jld2 file.
"""
function concarea_plot(data::Dict{String, Any}; plot_freq = 1000, number_of_bins = 0)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    plot_steps = 0:plot_freq:nsteps
    max_conc = [findmax(data["snapshots/Concentration/0"][:, :, i])[1] for i ∈ 1:nlayers]
    UpperConcentrationArea = Plots.Plot{Plots.GRBackend}[]
    LowerConcentrationArea = Plots.Plot{Plots.GRBackend}[]
    for i ∈ plot_steps
        upperdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 1], :)
        lowerdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 2], :)
        if number_of_bins == 0
            upperhist = fit(Histogram, upperdata)
            lowerhist = fit(Histogram, lowerdata)
        else
            upperhist = fit(Histogram, upperdata, nbins = number_of_bins)
            lowerhist = fit(Histogram, lowerdata, nbins = number_of_bins)
        end
        upperhist = normalize(upperhist, mode = :probability)
        lowerhist = normalize(lowerhist, mode = :probability)
        upperarea = vcat(0, cumsum(reverse(upperhist.weights)))
        lowerarea = vcat(0, cumsum(reverse(lowerhist.weights)))
        upperconc = reverse(Vector(upperhist.edges[1]))
        lowerconc = reverse(Vector(lowerhist.edges[1]))
        push!(UpperConcentrationArea, plot(upperarea, upperconc,
                                                label = false,
                                                xlabel = "Normalised area",
                                                ylabel = "Concentration",
                                                ylims = (0, max_conc[1])
                                            )
                )
        push!(LowerConcentrationArea, plot(lowerarea, lowerconc,
                                                label = false,
                                                xlabel = "Normalised area",
                                                ylabel = "Concentration",
                                                ylims = (0, max_conc[2])
                                            )
            )
    end
    return [UpperConcentrationArea, LowerConcentrationArea]
end
"""
    function concarea_animate(data)
Create an animation of Concetration ~ normalised area from the saved data in the output file.
"""
function concarea_animate(data::Dict{String, Any}; number_of_bins = 0)

    save_freq = data["save_freq"]
    if save_freq <  10
        plot_freq = 10 * save_freq
    else
        plot_freq = save_freq
    end
    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    max_conc = [findmax(data["snapshots/Concentration/0"][:, :, i])[1] for i ∈ 1:nlayers]
    ConcVsArea = @animate for i in 0:plot_freq:nsteps
        upperdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 1], :)
        lowerdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 2], :)
        if number_of_bins == 0
            upperhist = fit(Histogram, upperdata)
            lowerhist = fit(Histogram, lowerdata)
        else
            upperhist = fit(Histogram, upperdata, nbins = number_of_bins)
            lowerhist = fit(Histogram, lowerdata, nbins = number_of_bins)
        end
        upperhist = normalize(upperhist, mode = :probability)
        lowerhist = normalize(lowerhist, mode = :probability)
        upperarea = vcat(0, cumsum(reverse(upperhist.weights)))
        lowerarea = vcat(0, cumsum(reverse(lowerhist.weights)))
        upperconc = reverse(Vector(upperhist.edges[1]))
        lowerconc = reverse(Vector(lowerhist.edges[1]))
        p1 = plot(upperarea , upperconc,
                    label = false,
                    xlabel = "Normalised area",
                    ylabel = "Concentration",
                    ylims = (0, max_conc[1]),
                    title = "Upper layer"
                )
        p2 = plot(lowerarea, lowerconc,
                    label = false,
                    xlabel = "Normalised area",
                    ylabel = "Concentration",
                    ylims = (0, max_conc[2]),
                    title = "Lower layer"
                )
    plot(p1, p2)
    end

    return ConcVsArea
end
"""
    function tracer_plot(data)
Plot a heatmap of the concentration field at specified time steps from a tracer advection
diffusion simulation. The input is a loaded .jld2 output file.
"""
function tracer_plot(data::Dict{String, Any}; plot_freq = 1000)

    nsteps = data["clock/nsteps"]
    Lx, Ly = data["grid/Lx"], data["grid/Ly"]
    if Lx >= 500e3
        #This just makes the domain a little easier to read if a really large domain is being used
        set_xticks = (-Lx/2:round(Int, Lx/6):Lx/2, string.(-Int(Lx/2e3):round(Int, Lx/6e3):Int(Lx/2e3)))
        set_yticks = (-Ly/2:round(Int, Ly/6):Ly/2, string.(-Int(Ly/2e3):round(Int, Ly/6e3):Int(Ly/2e3)))
        plotargs = (
                    aspectratio = 1,
                    color = :deep,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlims = (-Lx/2, Lx/2),
                    ylims = (-Ly/2, Ly/2),
                    xticks = set_xticks,
                    yticks = set_yticks
                    )  
    else
        plotargs = (
                    aspectratio = 1,
                    color = :deep,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlims = (-Lx/2, Lx/2),
                    ylims = (-Ly/2, Ly/2)
                    ) 
    end
    
    plot_steps = 0:plot_freq:nsteps
    x, y = data["grid/x"], data["grid/y"]
    UpperTracerPlots = Plots.Plot{Plots.GRBackend}[]
    LowerTracerPlots = Plots.Plot{Plots.GRBackend}[]
    for i ∈ plot_steps
        uppertracer = heatmap(x, y, data["snapshots/Concentration/"*string(i)][:, :, 1],
                                title = "C(x,y,t) step = "*string(i); 
                                plotargs...
                            )
        push!(UpperTracerPlots, uppertracer)
        lowertracer = heatmap(x, y, data["snapshots/Concentration/"*string(i)][:, :, 2],
                                title = "C(x,y,t) step = "*string(i); 
                                plotargs...
                            )
        push!(LowerTracerPlots, lowertracer)
    end

    return [UpperTracerPlots, LowerTracerPlots]                   
end
"""
    function tracer_animate(data)
Turn the saved concentration data into an animation.
"""
function tracer_animate(data::Dict{String, Any})

    nsteps = data["clock/nsteps"]
    Lx, Ly = data["grid/Lx"], data["grid/Ly"]
    x, y = data["grid/x"], data["grid/y"]
    if Lx >= 500e3
        #This just makes the domain a little easier to read if a really large domain is being used
        set_xticks = (-Lx/2:round(Int, Lx/6):Lx/2, string.(-Int(Lx/2e3):round(Int, Lx/6e3):Int(Lx/2e3)))
        set_yticks = (-Ly/2:round(Int, Ly/6):Ly/2, string.(-Int(Ly/2e3):round(Int, Ly/6e3):Int(Ly/2e3)))
        plotargs = (
                    aspectratio = 1,
                    color = :deep,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlims = (-Lx/2, Lx/2),
                    ylims = (-Ly/2, Ly/2),
                    xticks = set_xticks,
                    yticks = set_yticks
                    )  
    else
        plotargs = (
                    aspectratio = 1,
                    color = :deep,
                    xlabel = "x",
                    ylabel = "y",
                    colorbar = true,
                    xlims = (-Lx/2, Lx/2),
                    ylims = (-Ly/2, Ly/2)
                    ) 
    end

    save_freq = data["save_freq"]
    if save_freq <=  10
        plot_freq = 10 * save_freq
    else
        plot_freq = save_freq
    end
    TracerAnimation = @animate for i ∈ 0:plot_freq:nsteps
        uppertracer = heatmap(x, y, data["snapshots/Concentration/"*string(i)][:, :, 1]',
                                title = "Upper layer, C(x,y,t) step = "*string(i); 
                                plotargs...
                            )
        lowertracer = heatmap(x, y, data["snapshots/Concentration/"*string(i)][:, :, 2]',
                                title = "Lower layer, C(x,y,t) step = "*string(i); 
                                plotargs...
                            )   

        plot(uppertracer, lowertracer, size = (900, 600))
    end

    return TracerAnimation
end
"""
    function time_vec(data::Dict{String, Any})
Create a time vector for plotting from the saved .jld2 output.
"""
function  time_vec(data::Dict{String, Any})

    nsteps = data["clock/nsteps"]
    save_freq = data["save_freq"]
    t = [data["snapshots/t/"*string(i)] for i in 0:save_freq:nsteps]

    return t
end
"""
    function tracer_area_avg(data::Dict{String, Any})
Calculate the average area of the tracer patch. This is done by transfroming the 
concentration into a function of area then computing 
    A = (∫Ad * C(Ad) dAd)/∫∫C dAd.
This may change but for now I will work with this and see where I get to.
"""
function tracer_area_avg(data::Dict{String, Any}; number_of_bins = 0)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    AreaVConcentration = Array{Float64}(undef, length(plot_steps), nlayers)
    for i ∈ plot_steps
        upperdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 1], :)
        lowerdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 2], :)
        if number_of_bins == 0
            upperhist = fit(Histogram, upperdata)
            upperhist = normalize(upperhist; mode = :pdf)
            lowerhist = fit(Histogram, lowerdata)
            lowerhist = normalize(lowerhist; mode = :pdf)
        else
            upperhist = fit(Histogram, upperdata, nbins = number_of_bins)
            upperhist = normalize(upperhist; mode = :pdf)
            lowerhist = fit(Histogram, lowerdata, nbins = number_of_bins)
            lowerhist = normalize(lowerhist; mode = :pdf)
        end
        
        Cupper = reverse(Vector(upperhist.edges[1]))
        Clower = reverse(Vector(lowerhist.edges[1]))
        Aupper = cumsum(reverse(upperhist.weights))
        Alower = cumsum(reverse(lowerhist.weights))
        dAupper = reverse(upperhist.weights)
        dAlower = reverse(lowerhist.weights)
        upperarea = sum( [ Aupper[i] * 0.5 * (Cupper[i+1] + Cupper[i]) * dAupper[i] for i in 1:length(Aupper) ] )
        lowerarea = sum( [ Alower[i] * 0.5 * (Clower[i+1] + Clower[i]) * dAlower[i] for i in 1:length(Alower) ] )
        
        #=
        uppertraceramount = sum( [0.5 * (Cupper[i+1] + Cupper[i]) * dAupper[i] for i in 1:length(Aupper) ])
        lowertraceramount = sum( [0.5 * (Clower[i+1] + Clower[i]) * dAlower[i] for i in 1:length(Alower) ])
        =#
        j = round(Int, i/saved_steps)
        #=AreaVConcentration[j + 1, :] .= [upperarea/uppertraceramount, 
                                         lowerarea/lowertraceramount]=#
        AreaVConcentration[j + 1, :] .= [upperarea, 
                                         lowerarea]

    end

    return AreaVConcentration

end
 """
    function tracer_area_percentile(data::Dict{String, Any})
Compute the percentile of area from a concentration field using 
        ∫A(C)dC over interval (Cmax, C₁)
where C₁ is some chosen value of concentration. By default the 
standard deviation of concentration at each time step is used for C₁
but this can also be set to false and some other quantile can be entered.
"""
function tracer_area_percentile(data::Dict{String, Any}; conc_min = 0.1, standard_dev = false, sd_multiple = 1)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    grid_area = data["grid/nx"] * data["grid/ny"]
    plot_steps = 0:saved_steps:nsteps
    area_percentiles = Array{Float64}(undef, length(plot_steps), nlayers)
    for i ∈ plot_steps
        upperdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 1], :)
        sort!(upperdata, rev = true)
        lowerdata = reshape(data["snapshots/Concentration/"*string(i)][:, :, 2], :)
        sort!(lowerdata, rev = true)
        max_conc = [findmax(data["snapshots/Concentration/"*string(i)][:, :, j])[1] for j ∈ 1:nlayers]
        
        if standard_dev == true && sd_multiple == 1
            upperarea = findlast(upperdata .> std(upperdata))
            lowerarea = findlast(lowerdata .> std(lowerdata))
        elseif standard_dev == true && sd_multiple != 1
            upperarea = findlast(upperdata .> std(upperdata) * sd_multiple)
            lowerarea = findlast(lowerdata .> std(lowerdata) * sd_multiple)
        else
            upperarea = findlast(upperdata .> conc_min * max_conc[1])
            lowerarea = findlast(lowerdata .> conc_min * max_conc[2])
        end

        k = round(Int, i/saved_steps)
        area_percentiles[k + 1, :] .= [upperarea/grid_area, lowerarea/grid_area]
    end

    return area_percentiles

end
"""
Calculate the average tracer_area_percentile from an ensemble simulation. Here the saved data is an array of dictionaries.
"""
function tracer_area_percentile(data::Array{Dict{String, Any}}; conc_min = 0.1, standard_dev = false, sd_multiple = 1)

    nlayers = data[1]["params/nlayers"] #These are the same value over all simuations
    nsteps = data[1]["clock/nsteps"]
    saved_steps = data[1]["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    no_of_sims = length(data)
    area_percentiles = Array{Float64}(undef, length(plot_steps), nlayers, no_of_sims)
    for i ∈ 1:no_of_sims
        temp = tracer_area_percentile(data[i]; conc_min)
        @. area_percentiles[:, :, i] = temp
    end
    return area_percentiles

end
"""
    function avg_ensemble_tracer_area(data::Array{Dict{String, Any}}; conc_min = 0.1) 
Calculate average growth of tracer area patch from an enemble simulation using `tracer_area_percentile`.
"""
function avg_ensemble_tracer_area(data::Array{Dict{String, Any}}; conc_min = 0.1)

    area_per = tracer_area_percentile(data; conc_min)
    no_of_sims = length(data)
    avg_area = area_per[:, :, 1]
    for i ∈ 2:no_of_sims
        @. avg_area += area_per[:, :, i]
    end
    @. avg_area /= no_of_sims
    return avg_area

end
"""
    function ensemble_average(data::Array{Dict{String, Any}}; conc_min = 0.1)
Calculate average tracer field then find growth of this ensemble average.
"""
function ensemble_average_area(data::Array{Dict{String, Any}}; conc_min = 0.1)

    ensemble_average = data[:, :, 1]
    no_of_sims = length(data)
    for i ∈ 2:no_of_sims
        @. ensemble_average += data[:, :, i]
    end
    @. ensemble_average /= no_of_sims
    ensemble_average_area = tracer_area_percentile(ensemble_average; conc_min)

    return ensemble_average_area

end
"""
    function exp_fit(data::Dict{String, Any}; conc_min = 0.1, tfitfinal = 100, tplot_length = 10)
Fit an exponential curve via least squares to the second stage of the growth of the area of 
the tracer patch as calculated from the `tracer_area_percentile` function.
The data is fitted from 1:tfitfinal which can be specified then the plot is calculated for length tplot_length.
"""
function exp_fit(data::Dict{String, Any}; conc_min = 0.1, tfitfinal = 100, tplot_length = 10)

    area_per = tracer_area_percentile(data; conc_min)
    t = time_vec(data)
    tplot = t[1:tfitfinal + tplot_length]
    t = t[1:tfitfinal]
    nlayers = data["params/nlayers"]
    A = Array{Float64}(undef, length(tplot), 2, nlayers)
    for i in 1:nlayers
        X = [ones(length(t)) t]
        M = log.(area_per[1:tfitfinal, i])
        fit = inv(X' * X) * (X' * M)
        @. A[:, :, i] = [tplot exp(fit[1]) * exp(fit[2] * tplot)]
    end

    return A
end
"""
    function linear_fit(data::Dict{String, Any}; conc_min = 0.1, tfitvals = [100, 250], tplot_length = [10 0])
Fit a linear curve via least squares to the third stage of the growth of the area of 
the tracer patch as calculated from the `tracer_area_percentile` function. The data is fitted over the interval
tfitvals and the length of the out plotting vectors can be specified by extra argument tplot_length.
"""
function linear_fit(data::Dict{String, Any}; conc_min = 0.1, tfitvals = [100, 250], tplot_length = [10 0])
    
    area_per = tracer_area_percentile(data; conc_min)
    t = time_vec(data)
    tplot = t[tfitvals[1] - tplot_length[1] : tfitvals[2] + tplot_length[2]]
    t = t[tfitvals[1] : tfitvals[2]]
    nlayers = data["params/nlayers"]
    A = Array{Float64}(undef, length(tplot), 2, nlayers)
    for i in 1:nlayers
        X = [ones(length(t)) t]
        M = area_per[tfitvals[1] : tfitvals[2], i]
        fit = inv(X' * X) * (X' * M)
        @. A[:, :, i] = [tplot fit[1] + fit[2] * tplot]
    end

    return A
end
"""
    function diffusivity(data::Dict{String, Any})
Calculate the diffusivity in physical space using the growth of area of tracer patch
"""
function diffusivity(data::Dict{String, Any}, time_vals::Matrix{Int64}; conc_min = 0.1)

    t = time_vec(data)
    area_per = tracer_area_percentile(data; conc_min)
    phys_params = nondim2dim(data)
    nlayers = data["params/nlayers"]
    diff = Array{Float64}(undef, 1, nlayers)

    Area = phys_params["Lx"] * phys_params["Ly"]
    for i ∈ 1:nlayers
        Area_inc = area_per[time_vals[i, 2], i] - area_per[time_vals[i, 1], i]
        no_of_seconds = (t[time_vals[i, 2]] / data["clock/dt"]) * phys_params["Δt"] - (t[time_vals[i, 1]] / data["clock/dt"]) * phys_params["Δt"]
        diff[i] = (Area * Area_inc) / no_of_seconds
    end

    return diff
end
"""
Cacluate diffusivity from average of ensemble simulation.
"""
function diffusivity(data::Array{Dict{String, Any}}, time_vals::Matrix{Int64}; conc_min = 0.1)

    t = time_vec(data[1])
    avg_area = avg_ensemble_tracer_area(data; conc_min)
    phys_params = nondim2dim(data[1])
    nlayers = data[1]["params/nlayers"]
    diff = Array{Float64}(undef, 1, nlayers)

    Area = phys_params["Lx"] * phys_params["Ly"]
    for i ∈ 1:nlayers
        Area_inc = avg_area[time_vals[i, 2], i] - avg_area[time_vals[i, 1], i]
        no_of_seconds = (t[time_vals[i, 2]] / data[1]["clock/dt"]) * phys_params["Δt"] - (t[time_vals[i, 1]] / data[1]["clock/dt"]) * phys_params["Δt"]
        diff[i] = (Area * Area_inc) / no_of_seconds
    end

    return diff
end
"""
    function nondim2dim(prob)
Translates parameters that have been set in the non-dimensional space (as I use in my thesis) for a QG problem
to the phsyical space based off mid-latitude values in metres and seconds. The values have defaults set.
"""
function nondim2dim(prob::FourierFlows.Problem;
                    Ω = 7.29e-5,     # Earth"s rotation
                    ϕ = π/3,         # Latitude
                    a = 6378e3,      # Earth's radius
                    g = 9.81,        # Gravity
                    H = 1500,        # Total depth (in metres)
                    ρ₁ = 1034,       # Density of top layer
                    ρ₂ = 1035,       # Density of bottom layer
                    )

    f₀ = 2*Ω*sin(ϕ)             # Coriolis computed from above values
    gprime = g*((ρ₂ - ρ₁)/ρ₂)   # Reduced gravity
    
    Ld = sqrt(gprime*H)/(f₀)    #Rossby deformation radius
    U = 0.1

    #Domain
    Lx̂, Lŷ = prob.grid.Lx, prob.grid.Ly
    Lx = Ld * Lx̂
    Ly = Ld * Lŷ

    #Parameters
    f̂₀, β̂, μ̂, ν̂ = prob.params.f₀, prob.params.β, prob.params.μ, prob.params.ν'
    f₀ = (U/Ld) * f̂₀
    β = (U/Ld^2) * β̂
    μ = (U/Ld) * μ̂
    ν = (U*Ld) * ν̂

    #Time
    Δt̂ = prob.clock.dt
    Δt = (Ld/U) * Δt̂

    return Dict("f₀" => f₀,
                "β"  => β,
                "μ"  => μ,
                "ν"  => ν,
                "Lx" => Lx,
                "Ly" => Ly,
                "Δt" => Δt
                )

end
"""
Compute the nondimensionalised time and length from the saved data of a advection diffusion simulation
"""
function nondim2dim(data::Dict{String, Any};
                    Ω = 7.29e-5,     # Earth"s rotation
                    ϕ = π/3,         # Latitude
                    a = 6378e3,      # Earth's radius
                    g = 9.81,        # Gravity
                    H = 1500,        # Total depth (in metres)
                    ρ₁ = 1034,       # Density of top layer
                    ρ₂ = 1035,       # Density of bottom layer
                    )
    
    f₀ = 2*Ω*sin(ϕ)             # Coriolis computed from above values
    gprime = g*((ρ₂ - ρ₁)/ρ₂)   # Reduced gravity
    
    Ld = sqrt(gprime*H)/(f₀)    #Rossby deformation radius
    U = 0.1

    #Domain
    Lx̂, Lŷ = data["grid/Lx"], data["grid/Ly"]
    Lx = Ld * Lx̂
    Ly = Ld * Lŷ

    #Time
    Δt̂ = data["clock/dt"]
    Δt = (Ld/U) * Δt̂

    #Diffusivity
    κ̂ = data["params/κ"]
    κ = (U * Ld) * κ̂

    return Dict("Lx" => Lx,
                "Ly" => Ly,
                "Δt" => Δt,
                "κ"  => κ
                )
end

end #module
