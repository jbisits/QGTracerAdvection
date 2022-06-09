using JLD2, Plots, Printf

## Testing output from `PassiveTracerFlows.jl`
cd(joinpath(pwd(), "PTF test"))

## One dimensional Gaussian

file = jldopen("advection-diffusion1D_timedep.jld2")

iterations = parse.(Int, keys(file["snapshots/t"]))
t = [file["snapshots/t/$i"] for i ∈ iterations]

c = [file["snapshots/concentration/$i"] for i ∈ iterations]
nothing # hide

x  = file["grid/x"]
Lx = file["grid/Lx"]

plot_args = (xlabel = "x",
             ylabel = "c",
             framestyle = :box,
             xlims = (-Lx/2, Lx/2),
             ylims = (0, maximum(c[1])),
             legend = :false,
             color = :balance)

p = plot(x, Array(c[1]), title = "concentration, t = " * @sprintf("%.2f", t[1]); plot_args...)

nothing # hide

# Create a movie of the tracer with the streamlines.

anim = @animate for i ∈ 1:length(t)
  p[1][:title] = "concentration, t = " * @sprintf("%.2f", t[i])
  p[1][1][:y] = Array(c[i])
end

mp4(anim, "1D_advection-diffusion.mp4", fps = 12)
close(file)
## Two dimensional turbulent

file = load("advection-diffusion_1.jld2")

iterations = #=parse.(Int, keys(file["snapshots/t"]))=# string.(0:25:4000)
t = [file["snapshots/t/$i"] for i ∈ iterations]

# Concentration and streamfunction time series in the lower layer
cₗ = [file["snapshots/concentration/$i"][:, :, 2] for i ∈ iterations]
ψₗ = [file["snapshots/streamfunction/$i"][:, :, 2] for i ∈ iterations]

# We normalize all streamfunctions to have maximum absolute value `amplitude/5`.
amplitude = 10
for i in 1:length(ψₗ)
  ψₗ[i] = amplitude/5 * ψₗ[i] / maximum(abs, ψₗ[i])
end

x,  y  = file["grid/x"],  file["grid/y"]
Lx, Ly = file["grid/Lx"], file["grid/Ly"]

plot_args = (xlabel = "x",
             ylabel = "y",
             aspectratio = 1,
             framestyle = :box,
             xlims = (-Lx/2, Lx/2),
             ylims = (-Ly/2, Ly/2),
             legend = :false,
             clims = (-amplitude/5, amplitude/5),
             colorbar_title = "\n concentration",
             color = :balance)

p = heatmap(x, y, Array(cₗ[1]'), title = "concentration, t = " * @sprintf("%.2f", t[1]); plot_args...)

contour!(p, x, y, Array(ψₗ[1]'), lw=2, c=:black, ls=:solid, alpha=0.7)

nothing # hide

# Create a movie of the tracer

anim = @animate for i ∈ 1:length(t)
  p[1][1][:z] = Array(cₗ[i])
  p[1][:title] = "concentration, t = " * @sprintf("%.2f", t[i])
  p[1][2][:z] = Array(ψₗ[i])
end

mp4(anim, "turbulentflow_advection-diffusion.mp4", fps = 12)

