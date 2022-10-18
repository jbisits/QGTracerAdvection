using JLD2, CairoMakie, Printf

cd(joinpath(pwd(), "Paper plots"))

ens_file = jldopen("ens_mean_conc.jld2")
iterations = sort(parse.(Int, keys(ens_file["snapshots/t"])))
t = [ens_file["snapshots/t/$i"] for i ∈ iterations]

c₁ = [abs.(ens_file["snapshots/Concentration/$i"][:, :, 1]) for i ∈ iterations]
c₂ = [abs.(ens_file["snapshots/Concentration/$i"][:, :, 2]) for i ∈ iterations]

close(ens_file)
## Animation of ensemble mean concentration
x,  y  = ens_file["grid/x"], ens_file["grid/y"]
Lx, Ly = ens_file["grid/Lx"], ens_file["grid/Ly"]

n = Observable(1)

c_anim = @lift Array(c₁[$n])
title = @lift @sprintf("concentration, t = %.2f", t[$n])

fig = Figure(resolution = (600, 600))
ax = Axis(fig[1, 1],
          xlabel = "x",
          ylabel = "y",
          aspect = 1,
          title = title,
          limits = ((-Lx/2, Lx/2), (-Ly/2, Ly/2)))

hm = heatmap!(ax, x, y, c_anim;
              colormap = :deep)

frames = 1:length(t)
record(fig, "ensemble_mean_conc.mp4", frames, framerate = 12) do i
    n[] = i
end

cn = contour!(ax, x, y, c_anim;
            colormap = :deep)
record(fig, "ensemble_mean_conc_contour.mp4", frames, framerate = 12) do i
    n[] = i
end

## Same plots as the individual member from the ensemble mean concentration field
plot_font = "CMU Modern Serif"
plot_fs = 20
latex_fs = 29
tracer_plots = Figure(resolution = (1200, 1400), fontsize = plot_fs, font = plot_font)
plot_steps = 0:3000:15000
plot_steps_mat = reshape(plot_steps, (2, 3))
plot_times = round.(Int, [ens_file["snapshots/t/"*string(i)] for i ∈ plot_steps])
plot_times = reshape(plot_times, (3, 2))'
x̂ = ens_file["grid/x"]
ŷ = ens_file["grid/y"]
conc_plot_data = [abs.(ens_file["snapshots/Concentration/"*string(plot_steps_mat[j, i])][:, :, 1]) for j ∈ 1:2, i ∈ 1:3]
plot_letters = ["(a)" "(b)" "(c)"; "(d)" "(e)" "(f)"]

newline = "\n"
ax = [Axis(tracer_plots[i, j],
        xlabel = L"\hat{x}",
        xlabelsize = latex_fs,
        ylabel = L"\hat{y}",
        ylabelsize = latex_fs,
        title = L"%$(plot_letters[i, j]) \quad \hat{t} = %$(string(plot_times[i, j]))",
        titlesize = latex_fs,
        aspect = 1
        ) for j ∈ 1:3, i ∈ 1:2]

for (i, axis) in enumerate(ax)

    plot_data = conc_plot_data[i]
    CairoMakie.heatmap!(axis, x̂, ŷ, plot_data, colormap = :deep)

end

for i ∈ 1:2, j ∈ 1:3

    plot_data = conc_plot_data[i, j]
    clims = (minimum(plot_data), maximum(plot_data))
    cticks = range(clims[1], clims[2], length = 4)
    ctickstyle = [@sprintf("%.2E", ticks) for ticks ∈ cticks]
    println(ctickstyle)
    Colorbar(tracer_plots[i, j][2, 1], label = L"Concentration ($\hat{C}$)", labelsize = latex_fs,
            limits = clims, ticks = cticks, tickformat = ct -> [@sprintf("%.2E", ticks) for ticks ∈ cticks],
            vertical = false, flipaxis = false, colormap = :deep, 
            ticklabelsize = plot_fs, ticklabelrotation = 45.0)

end
tracer_plots
#save("ens_tracer_plots.png", tracer_plots)

## Fit a Gaussian to the concentration

## Can use the package `LsqFit.jl`

using LsqFit, GLM

## Fit the ensemble mean concentration for a single time step

x,  y  = ens_file["grid/x"], ens_file["grid/y"]
xy = hcat(x, y)

# Target twoD Gaussian. The parameters fitted are the mean μ = (x₀, y₀),
# the covarance matrix Σ = (1/σ²)[1 0; 0 1], where we assume isotropic width
# and the amplitude.
function twoD_Gaussian(xy, p)
    amplitude, x₀, y₀, σ = p

    # creating linear meshgrid from xy
    x = repeat(xy, outer = (512, 1))[:, 1]
    y = repeat(xy, inner = (512, 1))[:, 1]
    g =  amplitude .* exp.( - (1 / (2 * σ^2)) .* (((x .- x₀).^2) + ((y .- y₀).^2)))

    return g[:]
end

# Initial guess for parameters
p₀ = [0.01, 3, 2, 0.2]
# Data to fit parameters for
conc_vals = reshape(c₁[1], :)
# Fit and check values
fit_Lsq = LsqFit.curve_fit(twoD_Gaussian, xy, conc_vals, p₀)
fit_Lsq.param

# Plotting to see if matches the data.
MvNormal(fit_Lsq.param[2:3], fit_Lsq.param[4])
test_fit(x, y) = pdf(MvNormal(fit_Lsq.param[2:3], fit_Lsq.param[4]), [x, y])

xgrid, ygrid = ones(length(x)) * x', y * ones(length(y))'
heatmap(x, y, test_fit.(xgrid, ygrid)')

## Fit the Gaussian for all time steps and save the parameters for upper layer

params_t = Array{Float64}(undef, length(p₀), length(c₁))

for i ∈ eachindex(c₁)

    @info "Snapshot $i"
    conc_vals = reshape(c₁[i], :)
    fit_Lsq = LsqFit.curve_fit(twoD_Gaussian, xy, conc_vals, p₀)
    @. params_t[:, i] = fit_Lsq.param

end

params_t
σ²_t = params_t[4, :].^2

## Lower layer
params_t_lower = Array{Float64}(undef, length(p₀), length(c₂))

for i ∈ eachindex(c₂)

    @info "Snapshot $i"
    conc_vals = reshape(c₂[i], :)
    fit_Lsq = LsqFit.curve_fit(twoD_Gaussian, xy, conc_vals, p₀)
    @. params_t_lower[:, i] = fit_Lsq.param

end

params_t_lower
σ²_t = params_t_lower[4, :].^2

# Save fitted parameters from upper and lower layer
#jldsave("fitted_params.jld2"; upper_layer=params_t, lower_layer=params_t_lower)

## Load in saved data

layer = "lower"
params_t = load("fitted_params.jld2")[layer*"_layer"]
σ²_t = params_t[end, :].^2

## Check linear growth, looks like best bet use from beginning to 
# around half way through
lines(t, σ²_t)

lines(t[1:90], σ²_t[1:90])

## Linear model for the linear growth (from this we get diffusivity)
t_half = t[1:90]
σ²_t_half = σ²_t[1:90]

t_mat = [ones(length(t_half)) t_half]
sec_mom_model = lm(t_mat, σ²_t_half)

int, slope = coef(sec_mom_model)
r2(sec_mom_model)

fig_sec_mom = Figure()
ax = Axis(fig_sec_mom[1, 1], 
          xlabel = "t̂", 
          ylabel = "σ̂²(t)",
          title = "Growth of second moment fitted to the concentration data")

lines!(ax, t_half, σ²_t_half, label = "Second moment (fitted to the data)")
lines!(ax, t_half, int .+ t_half .* slope, label = "Linear fit")
axislegend(ax; position = :lt)
fig_sec_mom

## Compute diffusivity from slope of linear fit
(slope / 2) * 0.02 * 29862 # = 5559.36 upper layer, 5666.378 lower layer

# slope / 2 because we are interested in the growth of the width of the Gaussian
# as in Abernathy (2013), 0.5 * (dσ² / dt) = κ

## Same plots as ensemble mean above but using the fitted mean and variance
# Clearly does not match exactly as this is isotropic but not bad
plot_font = "CMU Modern Serif"
plot_fs = 20
latex_fs = 29
tracer_plots = Figure(resolution = (1200, 1400), fontsize = plot_fs, font = plot_font)
plot_steps = 1:31:181
x̂ = ens_file["grid/x"]
ŷ = ens_file["grid/y"]
plot_letters = ["(a)" "(b)" "(c)"; "(d)" "(e)" "(f)"]

newline = "\n"
ax = [Axis(tracer_plots[i, j],
        xlabel = L"\hat{x}",
        xlabelsize = latex_fs,
        ylabel = L"\hat{y}",
        ylabelsize = latex_fs,
        title = L"%$(plot_letters[i, j]) \quad \hat{t} = %$(string(plot_times[i, j]))",
        titlesize = latex_fs,
        aspect = 1
        ) for j ∈ 1:3, i ∈ 1:2]

for (i, axis) in enumerate(ax)

    plot_data = [pdf(MvNormal(params_t[2:3, plot_steps[i]], params_t[4, plot_steps[i]]), [xgrid, ygrid]) for xgrid ∈ x, ygrid ∈ y]
    CairoMakie.heatmap!(axis, x̂, ŷ, plot_data, colormap = :balance)

end

tracer_plots