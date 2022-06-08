using JLD2, Plots, Printf

cd(joinpath(pwd(), "PTF test"))

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

