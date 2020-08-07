
###setup
using Interpolations, Plots, Parameters, Distributions, DataFrames, CSV, LinearAlgebra, Statistics, Random
using Optim: optimize

dir = "C:\\Users\\Garrett\\Documents\\Grad_School\\Papers\\TSO\\Model"
cd(dir)

include("tso_structs.jl")
include("tso_readin.jl")
include("tso_utilities.jl")
include("tso_background_functions.jl")
include("tso_model.jl")
include("tso_estimate.jl")
include("tso_simulate.jl")

#package = Readin_utilities(dir)
#package_moments = Readin_moments(dir)









###########
