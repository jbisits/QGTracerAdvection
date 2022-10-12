using JLD2, CairoMakie, Printf

cd(joinpath(pwd(), "Paper plots"))

ens_file = jldopen("ens_mean_conc.jld2")
iterations = sort(parse.(Int, keys(ens_file["snapshots/t"])))
t = [ens_file["snapshots/t/$i"] for i ∈ iterations]

c₁ = [abs.(ens_file["snapshots/Concentration/$i"][:, :, 1]) for i ∈ iterations]
c₂ = [abs.(ens_file["snapshots/Concentration/$i"][:, :, 2]) for i ∈ iterations]

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
save("ens_tracer_plots.png", tracer_plots)