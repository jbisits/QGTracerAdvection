#Load packages required for the simulation

ModulePath = joinpath(pwd(), "../../Modules")
include(joinpath(ModulePath, "TracerAdvDiff_QG.jl"))
include(joinpath(ModulePath, "MeasureMixing.jl"))
include(joinpath(ModulePath, "ExperimentSetup.jl"))

#import the local modules
using .TracerAdvDiff_QG, .MeasureMixing, .ExperimentSetup
#import required packages
using GeophysicalFlows.MultiLayerQG, Plots, Distributions, StatsBase, LinearAlgebra, JLD2
using Random: seed!