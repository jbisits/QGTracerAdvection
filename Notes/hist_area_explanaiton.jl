#How to calculate the histograms and then concentration and area data from each histograms

#=
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
=#