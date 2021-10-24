# Project the PV of a simulation into the tracer field.
# Take a snapshot of the PV then apply some function to it map it into the concentration field. 

## Load in the PV. This is only run for 10 steps. Just wanted the initial condition
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = QGPV, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

QGPV_init = data["snapshots/Concentration/"*string(0)]

x = data["grid/x"]
y = data["grid/y"]

heatmap(x, y, QGPV_init[:, :, 1]')

histogram(reshape(QGPV_init[:, :, 1], :), xlabel = "PV", ylabel = "ΔA")

histogram(reshape(QGPV_init[:, :, 1], :), xlabel = "PV", ylabel = "ΔA", normalize = :pdf)

## The histogram above looks quite Gaussian so could fit a Gaussian then normalise the PV data using the fitted values.

PV_1dGaussfit = fit(Normal, reshape(QGPV_init[:, :, 1], :))
μ = PV_1dGaussfit.μ
σ = PV_1dGaussfit.σ

#Normalise the QGPV_init with these values

QGPV_init_norm = @. (QGPV_init - μ) / σ
heatmap(x, y, QGPV_init_norm[:, :, 1]')

hist_1d_norm = histogram(reshape(QGPV_init_norm[:, :, 1], :), xlabel = "PV", ylabel = "ΔA", normalize = :pdf)

## Try meridional strips of Gaussians on the normalised data
QGPV_bandfit = similar(QGPV_init)
for i ∈ 1:length(x)
    tempfit = fit(Normal, QGPV_init_norm[i, :, 1])
    μ = tempfit.μ
    σ = tempfit.σ
    @. QGPV_bandfit[i, :, 1] = pdf(Normal(μ, σ), y)
end

heatmap(x, y, QGPV_bandfit[:, :, 1]')
# Does look like a band..
# Then a case of unnormalising to see how PV has been advected?

## Or could try covariance matrix?