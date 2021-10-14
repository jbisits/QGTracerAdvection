cd(joinpath(SimPath, "Output/Simulation: Domain = 128, res = 256, save_freq = 50, IC = GaussianStrip, Ensemble = true"))

#Load in the data. This is an ensemble simulation so now have an array of dictionaries.
data = Array{Dict{String, Any}}(undef, 10)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end

t = time_vec(data[1])
sec_mom = second_moment(data)

enss_conc = ensemble_concentration(data)
ense_second_mom = second_moment(enss_conc)

plot(t, ense_second_mom)

#Stuff below here may not be useful
## Average of area percentage growth
area_per = tracer_area_percentile(data; Cₚ = 0.5)
avg_area_per = avg_ensemble_tracer_area(data; Cₚ = 0.5)

upper_area = plot(t, area_per[:, 1, 1], 
                label = "Ensemble member 1", 
                title = "Growth of 90% of area of tracer patch (upper layer) \n from ensemble simulation with 10 memebers",
                legend = :topleft)
for i ∈ 2:length(data)
    plot!(upper_area, t, area_per[:, 1, i], 
        label = "Ensemble member "*string(i)
        )
end
plot!(upper_area, t, avg_area_per[:, 1], line = (:dash, 2, :black), label = "Average")

lower_area = plot(t, area_per[:, 2, 1], 
                label = "Ensemble member 1", 
                title = "Growth of 90% of area of tracer patch (lower layer) \n from ensemble simulation with 10 memebers",
                legend = :topleft)
for i ∈ 2:length(data)
    plot!(lower_area, t, area_per[:, 2, i], 
        label = "Ensemble member "*string(i)
        )
end  
plot!(lower_area, t, avg_area_per[:, 2], line = (:dash, 2, :black), label = "Average")

diff = diffusivity(data, [21 70; 21 70]; Cₚ = 0.5)

ensemble_area_per = plot(upper_area, lower_area, layout = (2, 1), size = (1200, 1200))
save("tenensemble_area_per.png", ensemble_area_per)
##

## Ensemble area (Growth of average of concentration field)
t = time_vec(data[1])
area_per = tracer_area_percentile(data; Cₚ = 0.999)
ensemble_conc = ensemble_concentration(data)
ensemble_area = tracer_area_percentile(ensemble_conc; Cₚ = 0.05)

ensemble_upper = plot(t, ensemble_area[:, 1], 
                    title = "Growth of 50% of ensemble area of tracer patch in upper layer",
                    legend = :topleft,
                    label = "Ensemble area",
                    line = (:dash, 2, :black)
                    )
for i ∈ 1:length(data)
    plot!(ensemble_upper, t, area_per[:, 1, i], 
        label = "Ensemble member "*string(i)
        )
end
ensemble_lower = plot(t, ensemble_area[:, 2], 
                    title = "Growth of 50% of ensemble area of tracer patch in lower layer",
                    legend = :topleft,
                    label = "Ensemble area",
                    line = (:dash, 2, :black)
                    )
for i ∈ 1:length(data)
    plot!(ensemble_lower, t, area_per[:, 2, i], 
        label = "Ensemble member "*string(i)
        )
end
ens_plot = plot(ensemble_upper, ensemble_lower, layout = (2,1), size = (1200, 1200))
save("ensembleavg.png", ens_plot)

ens_diff = diffusivity(data, [101 1; 101 1]; Cₚ = 0.05)

##
t = time_vec(data[1]; time_measure = "secs")
avg_area = tracer_avg_area(data[1])
plot(t, avg_area, label = false)

## Look at average area 

t = time_vec(data[1])
avg_area_members = tracer_avg_area(data)
plot(t, avg_area_members[:, 1, :], 
    title = "Average area of ensemble memebers and \n ensemble average from Gaussian strip (top layer)",
    legend = :bottomright,
    size = (800, 600))

ens_conc = ensemble_concentration(data)
avg_area_ens = tracer_avg_area(ens_conc)
plot!(t, avg_area_ens[:, 1], line = (:dash, 2, :black), label = "Ensemble")

## Now calculate diffusivity from the slope of the line
scatter!([t[70]], [avg_area_ens[70, 1]], color = :red, label = false)

K = avg_area_ens[70, 1] / (4 * π * t[70])
#K is a diffusivity and is translated into dimensions by (U * Ld) * K where U = 0.1m/s and Ld = 29862m
Kdim = (0.1 * 29862) * K

## 
t = time_vec(data[1])
second_moms = second_moment(data)
second_mom_upper = plot(t, second_moms[:, 1, 1], label = "Member 1", legend = :bottomright)
for i ∈ 2:length(data)
    plot!(second_mom_upper, t, second_moms[:, 1, i], label = "Memeber "*string(i))
end
second_mom_upper

ensemble_conc = ensemble_concentration(data)
ensemble_avg = second_moment(ensemble_conc)

plot!(second_mom_upper, t, ensemble_avg[:, 1], label = "Ensemble average", line = (:dash, 2, :black))