# Project the PV of a simulation into the tracer field.
# Take a snapshot of the PV then apply some function to it map it into the concentration field. 
# Method two the one I will first try to get going is to project the concentration field into PV

## Load in the PV. This is only run for 10 steps. Just wanted the initial condition
cd(joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = QGPV, Ensemble = false"))
file = joinpath(pwd(), "SimulationData.jld2")

#Load in the data
data = load(file)

QGPV_init = data["snapshots/Concentration/"*string(2000)]

x = data["grid/x"]
y = data["grid/y"]
Lₓ = data["grid/Lx"]
nx = data["grid/nx"]

heatmap(x, y, QGPV_init[:, :, 1]')

histogram(reshape(QGPV_init[:, :, 1], :), xlabel = "PV", ylabel = "ΔA")

## PV has a background gradient so is a quasi meridional coordinate. Want to generate PV()

## Try first just on a single meridional strip, so histogram is Δy ~ PV
PV_hist = fit(Histogram, QGPV_init[1, :, 1])
plot(PV_hist, xlabel = "PV", ylabel = "Δy", label = "Count of Δy")

## Cumulatively sum the weights
Δy_sum = similar(PV_hist.weights)
cumsum!(Δy_sum, PV_hist.weights)
Δy_sum = vcat(0, Δy_sum)
plot(Δy_sum, PV_hist.edges[1], xlabel = "Δy", ylabel = "ỹ = PV(y)", label = false)

# This could also be done using the sort method I think. Looks a little better?
sorted_PV = sort(QGPV_init[1, :, 1])
plot(y, sorted_PV, xlabel = "y", ylabel = "ỹ = q(y)", label = false)

# Now have a function PV = q(y) which can be thought of as the quasi-meridional coordinate ỹ.
# Want to connect this to tracer concentration experiments by C(ỹ). 
# This is done by finding some values q₁ < q₂ and summing all the concentration over these values.
# Can then link this to C by ∫Cdy over (q₁, q₂), so that eventually get C(q) all the concentration over some PV values.
# But have ỹ(PV) = ∫dy for q < q_*.

# Take PV between (q₁, q₂) = (120, 180).
q₁, q₂ = 120, 180
PV_vals = findall((sorted_PV .>= q₁ ) .& (sorted_PV .<= q₂))

# Now want to relate these values of PV to concentration from an advection-diffusion simulation.
## Load in some concentration data from a blob run
conc_path = joinpath(SimPath, "Output/Simulation: Lx̂ = Lŷ = 128, nx = 256, save_freq = 50, IC = GaussianBlob, Ensemble = true/SimulationData.jld2")
conc = load(conc_path)

conc_init = conc["snapshots/Concentration/"*string(0)]
heatmap(x, y, conc_init[:, :, 1]', color = :deep)

conc_PV = conc_init[:, PV_vals, 1]
heatmap(x, sorted_PV[PV_vals], conc_PV')

conc_PV = reshape(conc_init[:, PV_vals, 1], :)
Δx = conc["grid/Lx"] / conc["grid/nx"]
Δy = conc["grid/Ly"] / conc["grid/ny"]
ΔA = Δx * Δy

C_q = ΔA * sum( [conc_PV[i] for i in 1:length(conc_PV)] ) / ΔA

ỹ_q = sum( 1:length(PV_vals) )



## 
heatmap(x, y, QGPV_init[:, :, 1]') #Heatmap shows the meridional gradient of PV.
PV_hist = fit(Histogram, reshape(QGPV_init[:, :, 1], :))
plot(PV_hist, xlabel = "PV", ylabel = "ΔA")

# Want histogram of Δy ~ PV and have ΔA ~ PV. Divide by zonal resolution
PVy_weights = PV_hist.weights ./ nx
PVy_hist = fit(Histogram, PVy_weights, PV_hist.edges[1])
plot(PVy_hist)



##########################################################################################################################################
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


## By normalising the area is symmetrically about the origin with width 1 (or very close to..)
PV_normfitGaussian = fit(Normal, reshape(QGPV_init_norm[:, :, 1], :))

#Can I now use this information to generate the MV Gaussian? Below is not much use..

μ = [0, 0]
Σ = [500 0; 0 500]
testfit = MvNormal(μ, Σ)
generated_vals = [pdf(testfit, [X, Y]) for X in x, Y in y]

heatmap(x, y, generated_vals')

histogram(reshape(generated_vals, :))