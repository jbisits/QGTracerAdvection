#QGflow and create a histogram of the concnetration at each time step that the plot is created at

#Passive tracer advection using a two layer QG flow (from the geophysical flows package).

using .TracerAdvDiff_QG

using GeophysicalFlows.MultiLayerQG, Plots, Distributions, StatsBase, LinearAlgebra, JLD2

#Set up the MultiLayerQG.Problem to test with the modified module.

#Choose CPU or GPU
dev = CPU()

#Numerical and time-stepping parameters
nx = 64        # 2D resolution = nx^2
ny = nx

stepper = "FilteredRK4";  # timestepper
Δt = 0.01                 # timestep
nsubs  = 1                # number of time-steps for plotting (nsteps must be multiple of nsubs)
nsteps = 6000nsubs        # total number of time-steps


#Physical parameters for a two layer QG_problem
Lx = 2π        # domain size
μ = 5e-2       # bottom drag
β = 5          # the y-gradient of planetary PV

nlayers = 2     # number of layers
f0, g = 1, 1    # Coriolis parameter and gravitational constant
H = [0.2, 0.8]  # the rest depths of each layer
ρ = [4.0, 5.0]  # the density of each layer

U = zeros(nlayers) # the imposed mean zonal flow in each layer
U[1] = 1.0
U[2] = 0.0

#Setup QG_problem and make easier to access the parts of the struct
QG_prob = MultiLayerQG.Problem(nlayers, dev; nx=nx, Lx=Lx, f₀=f0, g=g, H=H, ρ=ρ, U=U, dt=Δt, stepper=stepper, μ=μ, β=β)

sol_QG, cl_QG, pr_QG, vs_QG = QG_prob.sol, QG_prob.clock, QG_prob.params, QG_prob.vars
x_QG, y_QG = QG_prob.grid.x, QG_prob.grid.y

#Set initial conditions.
ϵ = 0.3
x, y = gridpoints(QG_prob.grid)

q_1_i = @.  ϵ * cos(4π / Lx * x_QG) * exp(-(x^2 + y^2) / 8)
q_2_i = @. -ϵ * cos(4π / Lx * x_QG) * exp(-(x^2 + y^2) / 8)

q_i = zeros(nx,ny,2)

q_i[:, :, 1] = q_1_i
q_i[:, :, 2] = q_2_i
qh_i = QG_prob.timestepper.filter .* rfft(q_i, (1, 2))         # only apply rfft in dims=1, 2
q_i  = irfft(qh_i, QG_prob.grid.nx, (1, 2))                    # only apply irfft in dims=1, 2

MultiLayerQG.set_q!(QG_prob, q_i)

#Set diffusivity
κ = 0.01
#Set delay time (that is flow for t seconds, then drop tracer in)
delay_time = 0
#Set the tracer advection probelm by passing in the QG problem 
AD_prob = TracerAdvDiff_QG.Problem(;prob = QG_prob, delay_time = delay_time, nsubs = nsubs, κ = κ)
sol_AD, cl_AD, v_AD, p_AD, g_AD = AD_prob.sol, AD_prob.clock, AD_prob.vars, AD_prob.params, AD_prob.grid
x_AD, y_AD = gridpoints(g_AD)
x, y = g_AD.x, g_AD.y
#Set the (same) initial condition in both layers.

#A Gaussian blob centred at μIC 

μIC = [0, 0]
Σ = [1 0; 0 1]
blob = MvNormal(μIC, Σ)
blob_IC(x, y) = pdf(blob, [x, y])
C₀ = @. blob_IC(x_AD, y_AD)


#A Gaussian strip around centred at μIC.
#=
μIC = 0
σ² = 0.5
strip = Normal(μIC, σ²)
strip_IC(x) = pdf(strip, x)
C₀ = Array{Float64}(undef, g_AD.nx, g_AD.ny)
for i in 1:g_AD.nx
    C₀[i, :] = strip_IC(y_AD[i, :])
end
=#

#If using strip_IC use C₀' for a vertical strip
TracerAdvDiff_QG.QGset_c!(AD_prob, C₀)

max_conc = [findmax(AD_prob.vars.c[:, :, 1])[1], findmax(AD_prob.vars.c[:, :, 1])[1]]

#Define blank arrays in which to store the plots of tracer diffusion in each layer.
lower_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]
upper_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]

#Define frequency at which to save a plot.
#plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
#Setting them the same gives plots at equal time increments. (Might be a better work around)
plot_time_AD, plot_time_inc = 0.2, 0.2
#Blank array to save the step number so can plot histogram corresponding to the tracer advection plot.
step_nums = []

#Define a file for .jld2 to save into
filepath = "."
filename = joinpath(filepath, "Honours thesis/hist_conc.jld2")

#Step the tracer advection problem forward and plot at the desired time step.
while cl_AD.step <= nsteps
    if cl_AD.step == 0
        tp_u = heatmap(x, y, v_AD.c[:, :, 1]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-g_AD.Lx/2, g_AD.Lx/2),
                ylims = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)));
        push!(upper_layer_tracer_plots_AD, tp_u)
        tp_l = heatmap(x, y, v_AD.c[:, :, 2]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-g_AD.Lx/2, g_AD.Lx/2),
                ylims = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_layer_tracer_plots_AD, tp_l)
        push!(step_nums, AD_prob.clock.step)
    elseif round(Int64, cl_AD.step) == round(Int64, plot_time_AD*nsteps)
        tp_u = heatmap(x, y, v_AD.c[:, :, 1]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-g_AD.Lx/2, g_AD.Lx/2),
                ylims = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(upper_layer_tracer_plots_AD, tp_u)
        tp_l = heatmap(x, y, v_AD.c[:, :, 2]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlims = (-g_AD.Lx/2, g_AD.Lx/2),
                ylims = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_layer_tracer_plots_AD, tp_l)
        push!(step_nums, AD_prob.clock.step)
        global plot_time_AD += plot_time_inc
    end
    #I have this file saved so it is commented out for now.
    #MeasureMixing.fit_hist!(filename, AD_prob, number_of_bins = 30)
    stepforward!(AD_prob, nsubs)
    TracerAdvDiff_QG.QGupdatevars!(AD_prob)
    #Updates the velocity field in advection problem to the velocity field in the MultiLayerQG.Problem at each timestep.
    TracerAdvDiff_QG.vel_field_update!(AD_prob, QG_prob, nsubs)
end
#Need to set this up so this does not need to be hardcoded.
#Display the tracer advection in the upper layer.
plot_top = plot(upper_layer_tracer_plots_AD[1], upper_layer_tracer_plots_AD[2], 
                upper_layer_tracer_plots_AD[3], upper_layer_tracer_plots_AD[4],
                upper_layer_tracer_plots_AD[5], upper_layer_tracer_plots_AD[6])
     
#Display the tracer advection in the lower layer.
plot_bottom = plot(lower_layer_tracer_plots_AD[1], lower_layer_tracer_plots_AD[2], 
                   lower_layer_tracer_plots_AD[3], lower_layer_tracer_plots_AD[4],
                   lower_layer_tracer_plots_AD[5], lower_layer_tracer_plots_AD[6])

#Now can load the output of the .jld2 file created to create histograms and plots of Concentration ~ normalised area.

data = load("Honours thesis/hist_conc.jld2")

hist_top = Plots.Plot{Plots.GRBackend}[]
hist_bottom = Plots.Plot{Plots.GRBackend}[]
for i ∈ step_nums
    push!(hist_top, plot(data["Histograms/step"*string(i)][1], 
                             label = false, 
                            xlabel = "Concentration", 
                            ylabel = "Normalised area",
                             xlims = (0, max_conc[1] + 0.01)))
    push!(hist_bottom, plot(data["Histograms/step"*string(i)][2],
                             label = false, 
                            xlabel = "Concentration", 
                            ylabel = "Normalised area",
                             xlims = (0, max_conc[2] + 0.01)))
end
hist_top = plot(hist_top[1], hist_top[2], hist_top[3], hist_top[4], hist_top[5], hist_top[6])
hist_bottom = plot(hist_bottom[1], hist_bottom[2], hist_bottom[3], hist_bottom[4], hist_bottom[5], hist_bottom[6])

plot(plot_top, hist_top, layout=(2, 1), size=(1200, 1200))
plot(plot_bottom, hist_bottom, layout=(2, 1), size=(1200, 1200))

#Concentration ~ area example plot
plot(data["ConcentrationData/step400"][1], data["Histograms/step400"][1].edges,
         label = false,
        xlabel = "Normalised area",
        ylabel = "Concentration",
        ylims = (0, max_conc[1] + 0.01)
        )

#Animation of #Concentration ~ area 
#This animation is already saved so is commented out
#=
ConcVsArea = @animate for i in 1:10:nsteps
    p1 = plot(data["ConcentrationData/step"*string(i)][1], data["Histograms/step"*string(i)][1].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 ylims = (0, max_conc[1] + 0.01),
                 title = "Top layer"
                )
    p2 = plot(data["ConcentrationData/step"*string(i)][2], data["Histograms/step"*string(i)][2].edges,
                 label = false,
                xlabel = "Normalised area",
                ylabel = "Concentration",
                 ylims = (0, max_conc[2] + 0.01),
                 title = "Bottom layer"
                )
    plot(p1, p2)
end

mp4(ConcVsArea, "Movies/ConcVsArea.mp4", fps=18)
=#