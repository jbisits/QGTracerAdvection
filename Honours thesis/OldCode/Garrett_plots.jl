using Distributions: rand
#Garrett paper plots

using Plots

@. Ad(t) = 4 * π * Kₕ * t
@. Ap(t) = π * (Kₛ / γ) * exp(α₁ * γ * (t - (4 * γ)^(-1)))
@. At(t) = π * (Kₛ / γ) * exp(α₂ * γ * (t - (4 * γ)^(-1)))

Kₛ = 10e-2
Kₕ = 10e3
γ = 10e-6
α₁ = 0.9
α₂ = 0.65
t = range(0, 2*365*3600, step = 50)
Aprange = Ap(t) .< Ad(t)
Atrange = At(t) .< Ad(t)

area_plot = plot(t, Ad(t), label = "Ad", legend = :topleft, 
                xticks = ([0, 1*365*3600, 2*365*3600], ["0", "1", "2"]),
                xlabel = "Time in years")
plot!(area_plot, t[Aprange], Ap(t[Aprange]), label = "Ap")
plot!(area_plot, t[Atrange], At(t[Atrange]), label = "At")
scatter!(area_plot, [(4 * γ)^(-1)], [0], label = "t = (4γ)^(-1)")