####File that contains functions governing wages, job offers, preferences, and unobserved heterogeneity

#function that governs whether people receive teaching offers
function μ(prim::Primitives, param::Params, Ω::Array{Any,1}, l::Int64)
    #arguments: license, last-period teaching choice, and total experience
    @unpack J, poly_exp = prim
    @unpack μ = param
    e, dt = Ω[6], Ω[8]
    e_teach = e[J] #teaching experience
    teach_last = 1 * (dt == J) #dummy for having teached in pervious period

    #unpack stuff and get stuff ready
    γ_dt, γ_l, γ_exp = μ[1], μ[2], μ[3:(3+poly_exp-1)]
    index = 3 + poly_exp
    γ_l_int, γ_exp_int = μ[index], μ[(index+1):(index + 1 + poly_exp - 1)]
    exp_vec = []
    for i = 1:poly_exp
        push!(exp_vec, e_teach^i)
    end

    #construct probability
    prob = 0
    prob += γ_dt*teach_last + γ_l*l + sum(γ_exp.*exp_vec)
    if teach_last == 1 #taught last period
        prob += γ_l_int * l + sum(γ_exp_int.*exp_vec) #add interacted license/experience coefficients
    end

    #bound
    if prob>1.0
        prob = 1.0
    elseif prob<0.0
        prob = 0.0
    end

    prob #return deliverable
end

#occupation utility given parameters and current state
function util(prim::Primitives, param::Params, Ω::Array{Any,1}, j::Int64, J::Int64)

    ###take care of some unpleasant setup
    @unpack α, κ_j, γ_ju, ν, λ, σ_η = param
    ν_index = Ω[2][j+1] #location of taste state for occupation j -- in first go, will be 1 or 2
    X = Ω[1] #demographics
    γ, κ, σ_η = γ_ju[j], κ_j[j], σ_η[j]

    #begin construction
    util = α* exp(w(prim, param, Ω, j)) #start with wage
    util += X[2] * γ[1] #gender utility component
    util += X[3] * γ[2] #race component
    util -= κ * (j!=Ω[8]) #add job switching cost
    util += ν[ν_index] #add unobserved taste for occupation j

    #add teaching VA
    if j == J #in teaching
        util += λ * va(prim, param, Ω) #add preference for teaching output
    end

    util += α * (σ_η^2)/2 #account for expected wage shock
    util #return
end

#function that governs wage process
function w(prim::Primitives, param::Params, Ω::Array{Any,1}, j::Int64)
    @unpack J, nm, nm, poly_exp = prim
    @unpack γ_jw, δ_j = param
    X, θ, m, MA, l, exp = Ω[1], Ω[2][1], Ω[3], Ω[4], Ω[5], Ω[6][j] #relevant state variables
    γ, δ = γ_jw[j], δ_j[j]
    wage = θ * δ #initialize with ability component (note: automatically 0 if in teaching)

    #polynomial of experience
    exp_vec = []
    for i = 1:poly_exp
        push!(exp_vec, exp^i)
    end

    if j == J #in teaching
        #prepare parameters
        γ_int, γ_exp, γ_ma = γ[1], γ[2:(2+poly_exp)-1], γ[2+poly_exp]
        wage += (γ_int + γ_ma*MA + sum(γ_exp.*exp_vec)) #add terms
    elseif j!=J #not in teaching
        γ_int, γ_exp, γ_cq = γ[1], γ[2:(2+poly_exp)-1], γ[(2+poly_exp)] #intercept, experience, cq
        index = 2+ poly_exp + 1 #for brevity
        γ_s, γ_r, γ_m = γ[index], γ[index+1], γ[index+2] #sex, race, major

        wage += (γ_int + sum(γ_exp.*exp_vec)) #add terms
        wage += (γ_s * X[2] + γ_r * X[3] + γ_m * m) #add sex, race, teaching major effects
    end
    wage
end

#function that governs teaching-va
function va(prim::Primitives, param::Params, Ω::Array{Any,1})
    @unpack γ_0, δ_0 = param
    @unpack J, nm, poly_exp = prim
    X, θ, m, MA, l, exp = Ω[1], Ω[2][1], Ω[3], Ω[4], Ω[5], Ω[6][J] #relevant state variables
    γ, δ = γ_0, δ_0

    #setup
    γ_int, γ_exp, γ_cq = γ[1], γ[2:(2+poly_exp-1)], γ[2+poly_exp]
    index = 2 + poly_exp + 1
    γ_g, γ_r, γ_m, γ_MA, γ_l = γ[index], γ[index+1], γ[index+2], γ[index+3], γ[index+4]

    #polynomial of experience
    exp_vec = []
    for i = 1:poly_exp
        push!(exp_vec, exp^i)
    end

    #individuals make decisions based off their expected VA, so if they don't know thier state yet we don't have to add it.
    va = Ω[7] #idiosyncratic term. Note that value is set to zero if agent is in state corresponding to not knowing
    va += δ_0 * θ #add ability effect
    va += γ_int + sum(γ_exp .* exp_vec) + γ_cq * X[1] + γ_g*X[2] + γ_r*X[3] + γ_m*m + γ_MA*MA + γ_l*l
    va #return
end

#function that governs how costly it is to obtain a teaching masters, ignoring utliity shocks
function cost_MA(prim::Primitives, param::Params, Ω::Array{Any,1})
    @unpack γ_MA, δ_MA = param
    @unpack J = prim
    X, θ, m = Ω[1], Ω[2][1], Ω[3] #relevant state variables
    γ, δ = γ_MA, δ_MA

    #setup
    γ_int, γ_cq, γ_g, γ_r, γ_m = γ[1], γ[2], γ[3], γ[4], γ[5]

    #construct
    cost = δ * θ #initialize with ability effect
    cost += γ_int + γ_cq*X[1] + γ_g * X[2] + γ_r*X[3] + γ_m * m #add other stuff
    return exp(cost)
end


#function that governs how costly it is to obtain a license, ignoring utliity shocks
function cost_license(prim::Primitives, param::Params, Ω::Array{Any,1}; override::Int64 = 0)
    @unpack γ_l, δ_l = param
    @unpack J = prim
    X, θ, m, MA = Ω[1], Ω[2][1], Ω[3], Ω[4] #relevant state variables
    γ, δ = γ_l, δ_l

    if override == 1 #switch MA state if agetn is evaluating continuatino value from getting an MA
        MA = 1
    end

    #setup
    γ_int, γ_cq, γ_g, γ_r, γ_m, γ_MA = γ[1], γ[2], γ[3], γ[4], γ[5], γ[6]

    #construct
    cost = δ * θ #initialize with ability effect
    cost += (γ_int + γ_cq * X[1] + γ_g*X[2] + γ_r*X[3] + γ_m*m + γ_MA*MA)
    return exp(cost)
end

#function that governs how costly it is to obtain certain college majors
function util_major(prim::Primitives, param::Params, Ω::Array{Any,1}, i_m::Int64)
    @unpack γ_m, δ_m, ρ_m, ν = param
    @unpack J = prim  #not quite so simple with teh ν

    X, θ = Ω[1], Ω[2][1] #relevant state variables
    ν_index = Ω[2][i_m + 1] #index of taste state for occupation j -- in first go, will be 1 or 2
    ν_m = ν[ν_index] #unobserved preference for the individaul, based on their state
    γ, δ, ρ = γ_m[i_m], δ_m[i_m], ρ_m[i_m]


    γ_int, γ_cq, γ_g, γ_r = γ[1], γ[2], γ[3], γ[4]

    #construct
    util = δ*θ + ρ*ν_m #initialize with unobserved effects
    util += γ_int + γ_cq*X[1] + γ_g * X[2] + γ_r * X[3] #intercpet, gender, race\
    util
end
##########
