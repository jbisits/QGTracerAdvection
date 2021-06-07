using Distributions: rand
#Garrett paper plots

using Plots

@. Ad(t) = 4 * π * Kₕ * t
@. Ap(t) = π * (Kₛ / γ) * exp(α₁ * γ * (t - (4 * γ)^(-1)))
@. At(t) = π * (Kₛ / γ) * exp(α₂ * γ * (t - (4 * γ)^(-1)))

Kₛ = 10e-2
Kₕ = 10e3
γ = 10e-6
α₁ = 500
α₂ = 100
t = range(0, 40000, step = 1)
Aprange = Ap(t) .< Ad(t)
Atrange = At(t) .< Ad(t)

area_plot = plot(t, Ad(t), label = "Ad", legend = :topleft)
plot!(area_plot, t[Aprange], Ap(t[Aprange]), label = "Ap")
plot!(area_plot, t[Atrange], At(t[Atrange]), label = "At")
scatter!(area_plot, [(4 * γ)^(-1)], [0], label = "t = (4γ)^(-1)")