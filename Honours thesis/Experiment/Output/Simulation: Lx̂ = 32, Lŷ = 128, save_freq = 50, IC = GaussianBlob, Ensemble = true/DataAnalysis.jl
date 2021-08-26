#Change to correct directory
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = 32, Lŷ = 128, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

## Load in the data. This is an ensemble simulation so now have an array of dictionaries.
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


## Calculate the average area for each member of the ensemble and the ensemble average area
t = time_vec(data[1])
avg_area_members = tracer_avg_area(data)
plot(t, avg_area_members[:, 1, :], 
    title = "Average area of ensemble memebers and \n ensemble average from Gaussian blob (top layer)",
    legend = :bottomright,
    size = (800, 600))

ens_conc = ensemble_concentration(data)
avg_area_ens = tracer_avg_area(ens_conc)
plot!(t, avg_area_ens[:, 1], line = (:dash, 2, :black), label = "Ensemble")

## Now calculate diffusivity from the slope of the line
scatter!([t[70]], [avg_area_ens[70, 1]], color = :red, label = false)

K = avg_area_ens[70, 1] / (4 * π * t[70])
#K is a diffusivity and is translated into dimensions by U * Ld where U = 0.1m/s and Ld = 29862m
Kdim = (0.1 * 29862) * K