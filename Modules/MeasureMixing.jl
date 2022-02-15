#=
    Diagnostics to measure mixing and functions for data extraction from a .jld2 file created by a tracer 
    advection diffusion simulation.
=#
module MeasureMixing

export
    conc_mean,
    conc_var,
    hist_plot,
    concarea_plot,
    concarea_animate,
    tracer_plot,
    tracer_animate,
    time_vec,
    first_moment,
    second_moment,
    meridional_second_mom,
    tracer_area_percentile,
    avg_ensemble_tracer_area,
    ensemble_average_area,
    ensemble_concentration,
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
    function hist_plot(data)
Create plots of histograms at the same timesteps as the tracer plots from the saved data
in the output file. The input `data` is the loaded .jld2 file.
"""
function hist_plot(data::Dict{String, Any}; plot_freq = 1000, number_of_bins = 0, xlims_same = false)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    plot_steps = 0:plot_freq:nsteps
    max_conc = [findmax(data["snapshots/Concentration/0"][:, :, i])[1] for i ∈ 1:nlayers]
    histograms = Array{Plots.Plot{Plots.GRBackend}}(undef, length(plot_steps), nlayers)

    if xlims_same == true
        upperxlims = (0, max_conc[1])
        lowerxlims = (0, max_conc[2])
    else
        upperxlims = nothing
        lowerxlims = nothing
    end
    for i ∈ plot_steps

        for j ∈ 1:nlayers
            concentration = reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :) 
            if number_of_bins == 0
                hist = fit(Histogram, concentration)
            else
                hist = fit(Histogram, concentration, nbins = number_of_bins)
            end
        hist = normalize(hist; mode = :probability)
        k = round(Int, i / plot_freq) + 1
        histograms[k, j] = plot(hist,
                            label = false, 
                            xlabel = "Concentration", 
                            ylabel = "Normalised area",
                            xlims = upperxlims)
        end
    end

    return histograms
end
"""
    function concarea_plot(data)
Create plots of Concetration ~ grid cells at the same time steps as the tracer plots from the 
saved data in the output file. The input `data` is the loaded .jld2 file.
"""
function concarea_plot(data::Dict{String, Any}; plot_freq = 1000)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    plot_steps = 0:plot_freq:nsteps
    max_conc = [findmax(data["snapshots/Concentration/0"][:, :, i])[1] for i ∈ 1:nlayers]
    ConcentrationArea =  Array{Plots.Plot{Plots.GRBackend}}(undef, length(plot_steps), nlayers)

    for i ∈ plot_steps
        
        for j ∈ 1:nlayers

            concentration = reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :)
            sort!(concentration, rev = true)
            grid_cells = 1:length(concentration)
            k = round(Int, i / plot_freq) + 1
            ConcentrationArea[k, j] = plot(grid_cells , concentration,
                                                label = false,
                                                xlabel = "Grid cells",
                                                ylabel = "Concentration",
                                                ylims = (0, max_conc[j]),
                                                title = "Layer "*string(j))
        end
    end
    return ConcentrationArea 
end
"""
    function concarea_animate(data)
Create an animation of Concetration ~ grid cells from the saved data in the output file.
"""
function concarea_animate(data::Dict{String, Any})

    save_freq = data["save_freq"]
    if save_freq <  10
        plot_freq = 10 * save_freq
    else
        plot_freq = save_freq
    end
    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    max_conc = [findmax(data["snapshots/Concentration/0"][:, :, i])[1] for i ∈ 1:nlayers]
    cum_area_plot = Array{Plots.Plot{Plots.GRBackend}}(undef, nlayers)

    ConcVsArea = @animate for i in 0:plot_freq:nsteps

        for j ∈ 1:nlayers
            concentration = reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :)
            concentration = reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :)
            sort!(concentration, rev = true)
            grid_cells = 1:length(concentration)
            cum_area_plot[j] = plot(grid_cells, concentration,
                                label = false,
                                xlabel = "Grid cells",
                                ylabel = "Concentration",
                                ylims = (0, max_conc[j]),
                                title = "Layer "*string(j))
        end
        if nlayers % 2 == 0
            plot_layout = (Int(nlayers / 2), 2)
        else
            plot_layout = (nlayers, 3)
        end

        plot(cum_area_plot..., layout = plot_layout)
        
    end

    return ConcVsArea
end
"""
    function tracer_plot(data)
Plot a heatmap of the concentration field at specified time steps from a tracer advection
diffusion simulation. The input is a loaded .jld2 output file.
"""
function tracer_plot(data::Dict{String, Any}; plot_freq = 1000)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    Lx, Ly = data["grid/Lx"], data["grid/Ly"]
    plotargs = (
            aspectratio = Lx/Ly,
                color = :deep,
                xlabel = "x̂",
                ylabel = "ŷ",
                colorbar = true,
                colorbar_title = "Concentration",
                xlims = (-Lx/2, Lx/2),
                ylims = (-Ly/2, Ly/2),
                framestyle = :rectangle
                ) 
    
    plot_steps = 0:plot_freq:nsteps
    x, y = data["grid/x"], data["grid/y"]
    tracer_plots = Array{Plots.Plot{Plots.GRBackend}}(undef, length(plot_steps), nlayers)

    for i ∈ plot_steps

        for j ∈ 1:nlayers
            k = round(Int, i / plot_freq) + 1
            tracer_plots[k, j] = heatmap(x, y, abs.(data["snapshots/Concentration/"*string(i)][:, :, j])',
                                        title = "C(x,y,t) step = "*string(i); 
                                        plotargs...)
        end

    end

    return tracer_plots                  
end
"""
    function tracer_animate(data)
Turn the saved concentration data into an animation.
"""
function tracer_animate(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    Lx, Ly = data["grid/Lx"], data["grid/Ly"]
    x, y = data["grid/x"], data["grid/y"] 
    plotargs = (
                color = :deep,
                xlabel = "x̂",
                ylabel = "ŷ",
                colorbar = true,
                colorbar_title = " \nConcentration",
                xlims = (-Lx/2, Lx/2),
                ylims = (-Ly/2, Ly/2),
                aspectratio = Lx/Ly
                ) 

    save_freq = data["save_freq"]
    saved_conc = 0:save_freq:nsteps
    tracer_plots = Array{Plots.Plot{Plots.GRBackend}}(undef, nlayers)

    TracerAnimation = @animate for i ∈ saved_conc

        for j ∈ 1:nlayers
            tracer_plots[j] = heatmap(x, y, abs.(data["snapshots/Concentration/"*string(i)][:, :, j]'),
                                title = "Layer "*string(j)*" C(x,y,t)\nstep = "*string(i); 
                                plotargs...)
        end

        plot(tracer_plots..., size = (1000, 600))

    end

    return TracerAnimation
end
"""
    function time_vec(data::Dict{String, Any}; days = false)
Create a time vector for plotting from the saved .jld2 output.
By defualt the simulation time is what is used for this time vector.
To return the time vector as seconds use argument time_measure = secs
To return the time vector as days (in real time) use argrument time_measure = days.
"""
function  time_vec(data::Dict{String, Any}; time_measure = nothing)

    nsteps = data["clock/nsteps"]
    save_freq = data["save_freq"]
    t = [data["snapshots/t/"*string(i)] for i in 0:save_freq:nsteps]

    if time_measure == "days"
        phys_params = nondim2dim(data)
        days = 3600 * 24
        t = ((t ./ data["clock/dt"]) .* phys_params["Δt"]) ./ days
    elseif time_measure == "secs"
        phys_params = nondim2dim(data)
        t = (t ./ data["clock/dt"]) .* phys_params["Δt"]
    end

    return t
end
"""
    function first_moment(data::Dict{String, Any})
Calculate the average area of the tracer patch. This is done by  
    ``\frac{(∫A * C(A) dA)}{∫∫C dA} ≈ (ΔxΔy)\frac{ΣkCₖ}{ΣCₖ}``  for k ∈ 1:No. of grid cells
and is used to calculate the average area growth of the Gaussian blob.
"""
function first_moment(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    first_mom = Array{Float64}(undef, length(plot_steps), nlayers)
    Δx = data["grid/Lx"] / data["grid/nx"]
    Δy = data["grid/Ly"] / data["grid/ny"]
    ΔA = Δx * Δy
 
    for i ∈ plot_steps

        for j ∈ 1:nlayers

            C = abs.(reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :)) #Absolute value avoids the negative values
            sort!(C, rev = true)
            N = length(C)
            ΣkCₖ = ΔA * sum( [k * C[k] for k ∈ 1:N] )
            ΣCₖ = sum(C)
            l = round(Int, i/saved_steps) + 1
            first_mom[l, j] = ΣkCₖ / ΣCₖ

        end
                    
    end

    return first_mom

end

function first_moment(data::Array{Dict{String, Any}})

    nlayers = data[1]["params/nlayers"]
    nsteps = data[1]["clock/nsteps"]
    saved_steps = data[1]["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    first_mom = Array{Float64}(undef, length(plot_steps), nlayers, length(data))
    Δx = data[1]["grid/Lx"] / data[1]["grid/nx"]
    Δy = data[1]["grid/Ly"] / data[1]["grid/ny"]
    ΔA = Δx * Δy

    for i ∈ 1:length(data)

        for j ∈ plot_steps

            for l ∈ 1:nlayers

                C = abs.(reshape(data[i]["snapshots/Concentration/"*string(j)][:, :, l], :))#Absolute value avoids the negative values
                sort!(C, rev = true)
                N = length(C)
                ΣkCₖ = ΔA * sum( [k * C[k] for k ∈ 1:N] )
                ΣCₖ = sum(C)
                m = round(Int, j/saved_steps) + 1
                first_mom[m, l, i] = ΣkCₖ / ΣCₖ

            end

        end

    end

    return  first_mom
end
"""
Adds argument so the data is subset meridionally and zonally
"""
function first_moment(data::Dict{String, Any}, zonal_subset::Int64, meridional_subset::Int64)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    first_mom = Array{Float64}(undef, length(plot_steps), nlayers)
    Δx = data["grid/Lx"] / data["grid/nx"]
    Δy = data["grid/Ly"] / data["grid/ny"]
    ΔA = (Δx * zonal_subset) * (Δy * meridional_subset)
    x_shift = round(Int64, meridional_subset / 2)
    y_shift = round(Int64, zonal_subset / 2)
    x_length = length(data["snapshots/Concentration/0"][:, 1, 1])
    y_length = length(data["snapshots/Concentration/0"][1, :, 1])
 
    for i ∈ plot_steps

        for j ∈ 1:nlayers

            data_subset = [data["snapshots/Concentration/"*string(i)][x + x_shift, y + y_shift, j] 
                                for x ∈ 1:meridional_subset:x_length, y ∈ 1:zonal_subset:y_length]
            C = abs.(reshape(data_subset, :)) #Absolute value avoids the negative values
            sort!(C, rev = true)
            N = length(C)
            ΣkCₖ = ΔA * sum( [k * C[k] for k ∈ 1:N] )
            ΣCₖ = sum(C)
            l = round(Int, i/saved_steps) + 1
            first_mom[l, j] = ΣkCₖ / ΣCₖ

        end
                    
    end

    return first_mom

end

function first_moment(data::Array{Dict{String, Any}}, zonal_subset::Int64, meridional_subset::Int64)

    nlayers = data[1]["params/nlayers"]
    nsteps = data[1]["clock/nsteps"]
    saved_steps = data[1]["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    first_mom = Array{Float64}(undef, length(plot_steps), nlayers, length(data))
    Δx = data[1]["grid/Lx"] / data[1]["grid/nx"]
    Δy = data[1]["grid/Ly"] / data[1]["grid/ny"]
    x_shift = round(Int64, zonal_subset / 2)
    y_shift = round(Int64, meridional_subset / 2)
    x_length = length(data[1]["snapshots/Concentration/0"][:, 1, 1])
    y_length = length(data[1]["snapshots/Concentration/0"][1, :, 1])

    zonal_vec = []
    merid_vec = []
    ΔA = 0
    if zonal_subset == 0 && meridional_subset == 0
        zonal_vec = 1:x_length
        merid_vec = 1:y_length
        ΔA = Δx * Δy
    elseif zonal_subset == 0 && meridional_subset != 0
        zonal_vec = 1:x_length
        merid_vec = 1:meridional_subset:y_length
        ΔA = Δx * (Δy * meridional_subset)
    elseif zonal_subset != 0 && meridional_subset == 0
        zonal_vec = 1:zonal_subset:x_length
        merid_vec = 1:y_length
        ΔA = (Δx * zonal_subset) * Δy
    else
        zonal_vec = 1:zonal_subset:x_length
        merid_vec = 1:meridional_subset:y_length
        ΔA = (Δx * zonal_subset) * (Δy * meridional_subset)
    end

    for i ∈ 1:length(data)

        for j ∈ plot_steps

            for l ∈ 1:nlayers

                data_subset = [data[i]["snapshots/Concentration/"*string(j)][x + x_shift, y + y_shift, l] 
                                for x ∈ zonal_vec, y ∈ merid_vec]
                C = abs.(reshape(data_subset, :)) # Absolute value avoids the negative values
                sort!(C, rev = true)
                N = length(C)
                ΣkCₖ = ΔA * sum( [k * C[k] for k ∈ 1:N] )
                ΣCₖ = sum(C)
                m = round(Int, j/saved_steps) + 1
                first_mom[m, l, i] = ΣkCₖ / ΣCₖ

            end

        end

    end

    return  first_mom
end
"""
    function second_moment(data::Dict{Union{String, Any}}
Caclutate the second moment of the area from tracer data.
"""
function second_moment(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    second_mom = Array{Float64}(undef, length(plot_steps), nlayers)
    Δx = data["grid/Lx"] / data["grid/nx"]
    Δy = data["grid/Ly"] / data["grid/ny"]
    ΔA = Δx * Δy

    for i ∈ plot_steps

        for j ∈ 1:nlayers

            C = abs.(reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :)) #Absolute value avoids the negative values
            sort!(C, rev = true)
            N = length(C)
            l = round(Int, i/saved_steps) + 1
            Σk²Cₖ =  (ΔA)^2 * sum( [k^2 * C[k] for k ∈ 1:N] )
            ΣCₖ = sum(C)
            second_mom[l, j] = Σk²Cₖ / ΣCₖ

        end

    end

    return second_mom

end
function second_moment(data::Array{Dict{String, Any}})

    nlayers = data[1]["params/nlayers"]
    nsteps = data[1]["clock/nsteps"]
    saved_steps = data[1]["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    second_mom = Array{Float64}(undef, length(plot_steps), nlayers, length(data))
    Δx = data[1]["grid/Lx"] / data[1]["grid/nx"]
    Δy = data[1]["grid/Ly"] / data[1]["grid/ny"]
    ΔA = Δx * Δy

    for i ∈ 1:length(data)

        for j ∈ plot_steps

            for l ∈ 1:nlayers

                C = abs.(reshape(data[i]["snapshots/Concentration/"*string(j)][:, :, l], :)) #Absolute value avoids the negative values
                sort!(C, rev = true)
                N = length(C)
                m = round(Int, j/saved_steps) + 1
                Σk²Cₖ =  (ΔA)^2 * sum( [k^2 * C[k] for k ∈ 1:N] )
                ΣCₖ = sum(C)
                second_mom[m, l, i] = Σk²Cₖ / ΣCₖ

            end

        end

    end

    return  second_mom
end
"""
    function meridiondal_second_mom
Calculate merional second moment from the Gaussian strip condition.
"""
function meridional_second_mom(data::Dict{String, Any})

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    nx = data["grid/nx"]
    meridional_second_mom = Array{Float64}(undef, nx)
    second_mom = Array{Float64}(undef, length(plot_steps), nlayers)
    Δy = data["grid/Ly"] / data["grid/ny"]

    for j ∈ plot_steps

        for l ∈ 1:nlayers

            for n ∈ 1:nx

                C = abs.(reshape(data["snapshots/Concentration/"*string(j)][n, :, l], :)) #Absolute value avoids the negative values
                sort!(C, rev = true)
                N = length(C)
                Σk²Cₖ = (Δy / 2)^2 * sum( [k^2 * C[k] for k ∈ 1:N] )
                meridional_second_mom[n] = Σk²Cₖ

            end
            C = abs.(reshape(data["snapshots/Concentration/"*string(j)][:, :, l], :)) #Absolute value avoids the negative values
            ΣCₖ = sum(C)
            m = round(Int, j/saved_steps) + 1
            second_mom[m, l] = sum(meridional_second_mom) / ΣCₖ
        end

    end

    return  second_mom
end
function meridional_second_mom(data::Array{Dict{String, Any}})

    nlayers = data[1]["params/nlayers"]
    nsteps = data[1]["clock/nsteps"]
    saved_steps = data[1]["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    nx = data[1]["grid/nx"]
    meridional_second_mom = Array{Float64}(undef, nx)
    second_mom = Array{Float64}(undef, length(plot_steps), nlayers, length(data))
    Δy = data[1]["grid/Ly"] / data[1]["grid/ny"]

    for i ∈ 1:length(data)

        for j ∈ plot_steps

            for l ∈ 1:nlayers

                for n ∈ 1:nx

                    C = abs.(reshape(data[i]["snapshots/Concentration/"*string(j)][n, :, l], :)) #Absolute value avoids the negative values
                    sort!(C, rev = true)
                    N = length(C)
                    Σk²Cₖ =  (Δy / 2)^2 * sum( [k^2 * C[k] for k ∈ 1:N] )
                    meridional_second_mom[n] = Σk²Cₖ 

                end
                C = abs.(reshape(data[i]["snapshots/Concentration/"*string(j)][:, :, l], :))
                ΣCₖ = sum(C)
                m = round(Int, j/saved_steps) + 1
                second_mom[m, l, i] = sum(meridional_second_mom) / ΣCₖ
            end

        end

    end

    return  second_mom
end
 """
    function tracer_area_percentile(data::Dict{String, Any})
Compute the percentile of area from a concentration field using 
        ∫A(C)dC over interval (Cmax, Cₚ)
where Cₚ is some chosen value of concentration.
"""
function tracer_area_percentile(data::Dict{String, Any}; Cₚ = 0.5)

    nlayers = data["params/nlayers"]
    nsteps = data["clock/nsteps"]
    saved_steps = data["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    area_percentiles = Array{Float64}(undef, length(plot_steps), nlayers)
    for i ∈ plot_steps

        for j ∈ 1:nlayers
            
            C = reshape(data["snapshots/Concentration/"*string(i)][:, :, j], :)
            sort!(C, rev = true)
            cumsum_C = cumsum(C)
            C_total = sum(C)
            Nₚ = findfirst(cumsum_C .> Cₚ * C_total)
            l = round(Int, i/saved_steps) + 1
            area_percentiles[l, j] = Nₚ / ( data["grid/nx"] * data["grid/ny"] )
            
        end
        
    end

    return area_percentiles

end
"""
Calculate the tracer_area_percentile for each member of an ensemble simulation. Here the saved data is an array of dictionaries.
"""
function tracer_area_percentile(data::Array{Dict{String, Any}}; Cₚ = 0.5)

    nlayers = data[1]["params/nlayers"] #These are the same value over all simuations
    nsteps = data[1]["clock/nsteps"]
    saved_steps = data[1]["save_freq"]
    plot_steps = 0:saved_steps:nsteps
    no_of_sims = length(data)
    area_percentiles = Array{Float64}(undef, length(plot_steps), nlayers, no_of_sims)
    for i ∈ 1:no_of_sims
        temp = tracer_area_percentile(data[i]; Cₚ)
        @. area_percentiles[:, :, i] = temp
    end
    return area_percentiles

end
"""
    function avg_ensemble_tracer_area(data::Array{Dict{String, Any}}; Cₚ = 0.5) 
Calculate average growth of tracer area patch from an enemble simulation using `tracer_area_percentile`.
"""
function avg_ensemble_tracer_area(data::Array{Dict{String, Any}}; Cₚ = 0.5)

    area_per = tracer_area_percentile(data; Cₚ)
    no_of_sims = length(data)
    avg_area = area_per[:, :, 1]
    for i ∈ 2:no_of_sims
        @. avg_area += area_per[:, :, i]
    end
    @. avg_area /= no_of_sims
    return avg_area

end
"""
    function ensemble_concentration(data::Array{Dict{String, Any}})
Calculate average concentration from the tracer field of the ensemble simulation.
"""
function ensemble_concentration(data::Array{Dict{String, Any}})

    save_freq = data[1]["save_freq"]
    nsteps = data[1]["clock/nsteps"]
    saved_vals = 0:save_freq:nsteps
    ensemble_concentration = deepcopy(data[1])
    no_of_sims = length(data)

    for i ∈ 2:no_of_sims

        for j ∈ saved_vals
          @. ensemble_concentration["snapshots/Concentration/"*string(j)] += data[i]["snapshots/Concentration/"*string(j)]
        end

    end

    for j ∈ saved_vals
        @.  ensemble_concentration["snapshots/Concentration/"*string(j)]  /=  no_of_sims
    end

    return ensemble_concentration

end
"""
    function exp_fit(data::Dict{String, Any}; Cₚ = 0.5, tfitfinal = 100, tplot_length = 10)
Fit an exponential curve via least squares to the second stage of the growth of the area of 
the tracer patch as calculated from the `tracer_area_percentile` function.
The data is fitted from 1:tfitfinal which can be specified then the plot is calculated for length tplot_length.
"""
function exp_fit(data::Dict{String, Any}; Cₚ = 0.5, tfitfinal = 100, tplot_length = 10, time_measure = nothing)

    area_per = tracer_area_percentile(data; Cₚ = 0.5)
    t = time_vec(data; time_measure)
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
    function linear_fit(data::Dict{String, Any}; Cₚ = 0.5, tfitvals = [100, 250], tplot_length = [10 0])
Fit a linear curve via least squares to the third stage of the growth of the area of 
the tracer patch as calculated from the `tracer_area_percentile` function. The data is fitted over the interval
tfitvals and the length of the out plotting vectors can be specified by extra argument tplot_length.
"""
function linear_fit(data::Dict{String, Any}; Cₚ = 0.5, tfitvals = [100, 250], tplot_length = [10 0], time_measure = nothing)
    
    area_per = tracer_area_percentile(data; Cₚ)
    t = time_vec(data; time_measure)
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
function diffusivity(data::Dict{String, Any}, time_vals::Matrix{Int64}; Cₚ = 0.5)

    t = time_vec(data)
    area_per = tracer_area_percentile(data; Cₚ)
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
Cacluate diffusivity from average of ensemble simulation. First calls `ensemble_concentration`
to calcuate the ensemble concentration from an array of saved data, then calls `tracer_area_percentile`.
"""
function diffusivity(data::Array{Dict{String, Any}}, time_vals::Matrix{Int64}; Cₚ = 0.5)

    t = time_vec(data[1])
    ensemble_conc = ensemble_concentration(data)
    ensemble_area_per = tracer_area_percentile(ensemble_conc; Cₚ)
    phys_params = nondim2dim(data[1])
    nlayers = data[1]["params/nlayers"]
    diff = Array{Float64}(undef, 1, nlayers)

    Area = phys_params["Lx"] * phys_params["Ly"]
    for i ∈ 1:nlayers
        Area_inc = ensemble_area_per[time_vals[i, 2], i] - ensemble_area_per[time_vals[i, 1], i]
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
                    U = 0.02,         # Background current
                    Ω = 7.29e-5,     # Earth"s rotation
                    ϕ = π/3,         # Latitude
                    a = 6378e3,      # Earth's radius
                    g = 9.81,        # Gravity
                    H = 1500,        # Total depth (in metres)
                    ρ₁ = 1034,       # Density of top layer
                    ρ₂ = 1035        # Density of bottom layer
                    )

    f₀ = 2*Ω*sin(ϕ)             # Coriolis computed from above values
    gprime = g*((ρ₂ - ρ₁)/ρ₂)   # Reduced gravity
    
    Ld = sqrt(gprime*H)/(f₀)    #Rossby deformation radius

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
                "Δt" => Δt,
                "Ld" => Ld
                )

end
"""
Compute the nondimensionalised time and length from the saved data of a advection diffusion simulation
"""
function nondim2dim(data::Dict{String, Any};
                    U = 0.02,         # Background current
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
                "κ"  => κ,
                "Ld" => Ld
                )
end

end #module
