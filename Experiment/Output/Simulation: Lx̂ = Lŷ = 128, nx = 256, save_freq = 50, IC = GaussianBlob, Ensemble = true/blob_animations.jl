# Animations for the Gaussian blob

cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true"))

data = Array{Dict{String, Any}}(undef, 50)
for i ∈ 1:length(data)
    if i == 1
        file = joinpath(pwd(), "SimulationData.jld2")
        data[i] = load(file)
    else
        file = joinpath(pwd(), "SimulationData_"*string(i - 1)*".jld2")
        data[i] = load(file)
    end
end

## 

ensemble_conc = ensemble_concentration(data)

# The layout of the plots still needs to be fixed! 

## Single member

member_anim = tracer_animate(data[1])
mp4(member_anim, "member_adv_diff.mp4", fps = 18)

## Ensemble average concentration

ens_anim = tracer_animate(ensemble_conc)
mp4(ens_anim, "ensemble_adv_diff.mp4", fps = 18)
