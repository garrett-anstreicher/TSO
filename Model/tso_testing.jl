
###setup
using Interpolations, Plots, Parameters, Distributions, DataFrames, CSV, LinearAlgebra, Statistics, Random
using Optim: optimize

dir = "C:\\Users\\Garrett\\Documents\\Grad_School\\Papers\\TSO\\Model"
cd(dir)

include("tso_readin.jl")
include("tso_structs.jl")
include("tso_utilities.jl")
include("tso_background_functions.jl")
include("tso_model.jl")
#include("tso_simulate.jl")
#include("tso_estimate.jl")

#initialize parameters and primitives for the first time
guess_init = Param_init(2, 3, 2, 2)
@elapsed Solve_model(guess_init)

###########
