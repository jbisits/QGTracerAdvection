#QGflow and create a histogram of the concnetration at each time step that the plot is created at

#Passive tracer advection using a two layer QG flow (from the geophysical flows package).

using .TracerAdvDiff_QG

using GeophysicalFlows.MultiLayerQG, Plots, Distributions, StatsBase, LinearAlgebra

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

initial_data = reshape(AD_prob.vars.c[:, :, 1], :) #Get this here to play with down below

#Plot of initial condition in the upper layer.
IC_upper = heatmap(x, y, v_AD.c[:, :, 1]',
            title = "Upper layer initial tracer concentration",
            xlabel = "x",
            ylabel = "y",
            color = :balance,
            aspecetratio = 1,
            colorbar = true,
            xlim = (-g_AD.Lx/2, g_AD.Lx/2),
            ylim = (-g_AD.Ly/2, g_AD.Ly/2))
IC_lower = heatmap(x, y, v_AD.c[:, :, 2]',
            title = "Lower layer initial tracer concentration",
            xlabel = "x",
            ylabel = "y",
            color = :balance,
            aspecetratio = 1,
            colorbar = true,
            xlim = (-g_AD.Lx/2, g_AD.Lx/2),
            ylim = (-g_AD.Ly/2, g_AD.Ly/2))
plot(IC_upper, IC_lower, size = (900, 400))

#Define blank arrays in which to store the plots of tracer diffusion in each layer.
lower_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]
upper_layer_tracer_plots_AD = Plots.Plot{Plots.GRBackend}[]

#Define blank arrays to save a histogram at each time step a plot is saved.
upper_concentration_hist = Plots.Plot{Plots.GRBackend}[]
lower_concentration_hist = Plots.Plot{Plots.GRBackend}[]

#Define frequency at which to save a plot.
#plot_time_AD is when to get the first plot, plot_time_inc is at what interval subsequent plots are created.
#Setting them the same gives plots at equal time increments. (Might be a better work around)
plot_time_AD, plot_time_inc = 0.2, 0.2
#Define arguments for plots.
kwargs = (
         xlabel = "x",
         ylabel = "y",
    aspectratio = 1,
          color = :balance,
       colorbar = true,
           xlim = (-g_AD.Lx/2, g_AD.Lx/2),
           ylim = (-g_AD.Ly/2, g_AD.Ly/2)
)

#Step the tracer advection problem forward and plot at the desired time step.
while cl_AD.step <= nsteps
    if cl_AD.step == 0
        tp_u = heatmap(x, y, v_AD.c[:, :, 1]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)));
        push!(upper_layer_tracer_plots_AD, tp_u)
        hist_upper = histogram(reshape(AD_prob.vars.c[:, :, 1], :), label = false, normalize = :probability,
                               title = "Normalised histogram of \ntracer concentration, t = "*string(round(cl_AD.t; digits = 2)))
        push!(upper_concentration_hist, hist_upper)
        tp_l = heatmap(x, y, v_AD.c[:, :, 2]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_layer_tracer_plots_AD, tp_l)
        hist_lower = histogram(reshape(AD_prob.vars.c[:, :, 2], :), label = false, normalize = :probability,
                               title = "Normalised histogram of \ntracer concentration, t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_concentration_hist, hist_lower)
    elseif round(Int64, cl_AD.step) == round(Int64, plot_time_AD*nsteps)
        tp_u = heatmap(x, y, v_AD.c[:, :, 1]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(upper_layer_tracer_plots_AD, tp_u)
        hist_upper = histogram(reshape(AD_prob.vars.c[:, :, 1], :), label = false, normalize = :probability,
                               title = "Normalised histogram of \ntracer concentration, t = "*string(round(cl_AD.t; digits = 2)))
        push!(upper_concentration_hist, hist_upper)
        tp_l = heatmap(x, y, v_AD.c[:, :, 2]',
                aspectratio = 1,
                c = :balance,
                xlabel = "x",
                ylabel = "y",
                colorbar = true,
                xlim = (-g_AD.Lx/2, g_AD.Lx/2),
                ylim = (-g_AD.Ly/2, g_AD.Ly/2),
                title = "C(x,y,t), t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_layer_tracer_plots_AD, tp_l)
        hist_lower = histogram(reshape(AD_prob.vars.c[:, :, 2], :), label = false, normalize = :probability,
                                       title = "Normalised histogram of \ntracer concentration, t = "*string(round(cl_AD.t; digits = 2)))
        push!(lower_concentration_hist, hist_lower)
        global plot_time_AD += plot_time_inc
    end
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

#Histograms in upper layer
hist_top = plot(upper_concentration_hist[1], upper_concentration_hist[2],
                upper_concentration_hist[3], upper_concentration_hist[4],
                upper_concentration_hist[5], upper_concentration_hist[6])

#Histograms in lower layer
hist_bottom = plot(lower_concentration_hist[1], lower_concentration_hist[2],
                   lower_concentration_hist[3], lower_concentration_hist[4],
                   lower_concentration_hist[5], lower_concentration_hist[6])

plot(plot_top, hist_top, layout=(2, 1), size=(1200, 1200))
plot(plot_bottom, hist_bottom, layout=(2, 1), size=(1200, 1200))

#Want to integrate (cumumlative sum) over the historgram to get the area of the tracer patch
#To do this need to fit the histogram as an object rather than just using the histogram from Plots

histi = fit(Histogram, initial_data, nbins = 30)
probhisti = normalize(histi, mode = :probability)
#Now cumlulative sum to each bin and will get data that you can plot.
#Use cumlulative sum from highest concentration to lowest concentration. That is why need reverse argument on the weights.
hist_data = cumsum(reverse(probhisti.weights)) 
hist_data = vcat(0, hist_data)
plot(probhisti, label = false)
plot!(probhisti.edges, reverse(hist_data), label = false, linewidth = 4, color = :red,
        xlabel = "Concentration", ylabel = "Normalised area")
#Looks to just be a swap between the x and y axes to get the correct plot.
initial_conc_area = plot(reverse(hist_data), probhisti.edges, label = "Initial concentration level over grid", 
                        xlabel = "Normalised area", ylabel = "Concentration")


final_data = reshape(AD_prob.vars.c[:, :, 1], :)
histf = fit(Histogram, final_data)
probhistf = normalize(histf, mode = :probability)
hist_dataf = cumsum(reverse(probhistf.weights))
hist_dataf = vcat(0, hist_dataf)
#This plots only the part of the domain where there is concentration so need to make zero 
#either side of this. In reality would asymptote to zero as a Gaussian but will start by just
#setting it to zero. 
plot(probhistf, label = false)
plot!(probhistf.edges, reverse(hist_dataf), label = false, linewidth = 4, color = :red,
    xlabel = "Concentration", ylabel = "Normalised area")
final_conc_area = plot(reverse(hist_dataf), probhistf.edges, label = "Initial concentration level over grid", 
                        xlabel = "Normalised area", ylabel = "Concentration")

#I think a plot of this at every time step (then animate) would show the mixing of the tracer into the domain as we want it to.
#Running this for longer the line stays horizontal!! That means that the mixing in the domain is linear (I think). On the right track!
plot!(initial_conc_area, reverse(hist_dataf), probhistf.edges, label = "Final concentration level over grid", xlabel = "Normalised area", ylabel = "Concentration")

#This now needs to be done at every timestep. Could then animate it to get a movie which I think would show things best.
