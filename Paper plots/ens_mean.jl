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

