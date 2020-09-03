#struct for model primitives and state spaces
@with_kw struct Primitives
    β::Float64 = 0.95 #discount rate
    Π::Float64 = 0.95^4 #college-to-work discount rate (4 years)
    J::Int64 = 2 #number of occupations (home=0, non-teaching, teaching)
    T::Int64 = 10 #number of periods
    poly_exp::Int64 = 3 #degree of polynomial for experience

    ######state spaces and dimensions######
    #demographics and college major
    cq_grid::Array{Int64,1} = collect(0:1:1) #college quality (continuous)
    gender_grid::Array{Int64,1} = [0, 1] #gender state space. 1=female
    race_grid::Array{Int64,1} = [0, 1] #race (asian/white, black/hispanic)
    ng::Int64 = length(gender_grid)
    ncq::Int64 = length(cq_grid)
    nrace::Int64 = length(race_grid)

    #unobserved ability and preferences
    θ_grid::Array{Int64,1} = [1,2] #unobserved ability
    ν_grid::Array{Int64,1} = [1,2] #low, high taste for each occupation
    nθ::Int64 = length(θ_grid)
    nν::Int64 = length(ν_grid)

    #college major, teacher MA degree and teaching license
    m_grid::Array{Int64,1} = [0,1] #0 = not teaching. 1 = teaching
    MA_grid::Array{Int64,1} = [0,1] #0 = no master, 1 = master
    l_grid::Array{Int64,1} = [0,1] #dummy
    nm::Int64 = length(m_grid)
    nMA::Int64 = length(MA_grid)
    nl::Int64 = length(l_grid)

    #teacher quality
    ξ_grid::Array{Float64,1} = [0.0, -1.0, 0.0, 1.0] #levels of teaching quality
    nξ::Int64 = length(ξ_grid) #additioanl state for teaching quality being unknown

    #previous job and teaching job offer
    dt_grid::Array{Int64,1} = collect(0:1:J) #number of possible previous occupations
    𝒥_grid::Array{Int64,1} = [0,1] #dummy for offer in current period
    ndt::Int64 = length(dt_grid)
    n𝒥::Int64 = length(𝒥_grid)
end

@with_kw struct Primitives_collect
    T::Int64 = 10 #working periods

    #demographics
    X_grid::Array{Array{Int64,1}} = [] #all together
    nX::Int64 = length(X_grid)

    #heterogeneity
    χ_grid::Array{Array{Int64,1}} = [] #all together
    nχ::Int64 = length(χ_grid)

    #experience
    e_grid::Array{Array{Int64,1}} = [] #experience levels; a vector of vectors
    ne::Int64 = length(e_grid)
end

#struct for model parameters to be estimated. Going to need something nicer than what I usually write . . .
@with_kw mutable struct Params
    guess = zeros(20) #parameter guess
    γ_jw::Vector{Vector{Float64}} = guess[1]
    δ_j::Vector{Float64} = guess[2]
    σ_η::Vector{Float64} = guess[3]
    γ_0::Vector{Float64} = guess[4]
    δ_0::Float64 = guess[5]
    σ_ς::Float64 = guess[6]
    σ_ξ::Float64 = guess[7]
    α::Float64 = guess[8]
    λ::Float64 = guess[9]
    ν::Vector{Float64} = guess[10]
    γ_ju::Vector{Vector{Float64}} = guess[11]
    κ_j::Vector{Float64} = guess[12]
    γ_MA::Vector{Float64} = guess[13]
    δ_MA::Float64 = guess[14]
    γ_l::Vector{Float64} = guess[15]
    δ_l::Float64 = guess[16]
    γ_m::Vector{Vector{Float64}} = guess[17]
    δ_m::Vector{Float64} = guess[18]
    ρ_m::Vector{Float64} = guess[19]
    μ::Vector{Float64} = guess[20]
end

#struct for holding model value functions
#phases c and e are when teaching offers arise and people work. no choices to be made.
mutable struct Results
    v_coll::SharedArray{Float64, 3} #college value function
    v_work_a::SharedArray{Float64, 10} #working value function, phase 1
    v_work_b::SharedArray{Float64, 10} #working value function, phase 2
    v_work_d::SharedArray{Float64, 11} #working value function, phase 4
end

mutable struct Results_iter
    v_coll::Array{Float64, 1} #college value function
    v_work_a::Array{Float64, 8} #working value function, phase 1
    v_work_b::Array{Float64, 8} #working value function, phase 2
    v_work_d::Array{Float64, 9} #working value function, phase 4
end

#non-parallel results struct
mutable struct Results_npar
    v_coll::Array{Float64, 3} #college value function
    v_work_a::Array{Float64, 10} #working value function, phase 1
    v_work_b::Array{Float64, 10} #working value function, phase 2
    v_work_d::Array{Float64, 11} #working value function, phase 4
end
