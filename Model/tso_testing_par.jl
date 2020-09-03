using Distributed
#add processes
workers()
addprocs(3)


###setup
@everywhere using Interpolations, Plots, Parameters, Distributions, DataFrames, CSV, LinearAlgebra, Statistics, Random, SharedArrays
@everywhere using Optim: optimize

@everywhere dir = "C:\\Users\\Garrett\\Documents\\Grad_School\\Papers\\TSO\\Model"
@everywhere cd(dir)

@everywhere include("tso_readin.jl")
@everywhere include("tso_structs.jl")
@everywhere include("tso_utilities.jl")
@everywhere include("tso_background_functions.jl")
@everywhere include("tso_model_par.jl")
@everywhere include("tso_model.jl")
#@everywhere include("tso_simulate_par.jl")
#include("tso_estimate.jl")


#initialize parameters and primitives for the first time
@everywhere guess_init = Param_init(2, 3, 2, 2)
#Solve_model(guess_init; nsim=100000)
#@elapsed data_simul = Solve_model(guess_init; nsim=10000)


#=
@everywhere function testpar(guess::Array{Any,1})
    #obtain vanilla value functions
    println("Solving vanilla value functions ...")
    prim, prim_grp, param, res = Initialize(guess_init)
    Backward_induct(prim, prim_grp, param, res) #backward induction protocol

    println("Solving parallelized value functions ...")
    prim_par, prim_grp_par, param_par, res_par = Initialize_par(guess_init)
    Backward_induct_par(prim_par, prim_grp_par, param_par, res_par) #backward induction protocol

    println("Checking for equivalence ...")
    diff_coll = res.v_coll .- res_par.v_coll
    println("Maximum difference for coll fvunc: ", maximum(diff_coll))

    diff_a = res.v_work_a .- res_par.v_work_a
    println("Maximum difference for v_work_a fvunc: ", maximum(diff_a))

    diff_b = res.v_work_b .- res_par.v_work_b
    println("Maximum difference for v_work_b fvunc: ", maximum(diff_b))

    diff_d = res.v_work_d .- res_par.v_work_d
    println("Maximum difference for v_work_d fvunc: ", maximum(diff_d))
    println("All done!")
end

testpar(guess_init)
=#
###########
