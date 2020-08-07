###equal mass discretization of normal
function discretize_normal(μ_a::Float64, σ_a::Float64, n::Int64)
    dist = Normal(μ_a, σ_a) #set up log normal distribution
    draws= zeros(n) #preallocate market luck states

    for i = 1:n #begin loop to fill discretized vector
        quant = (2*i-1)/(2*n) #desired quantile
        draws[i] = quantile(dist, quant) #fill element of state vector
    end
    draws #return deliverable
end

#function for getting interpolated index of a given Float with an arbitrary grid, assuming linear interpolation
function get_index(val::Float64, grid::Array{Float64,1})
    n = length(grid)
    index = 0 #preallocation
    if val<=grid[1] #LEQ smallest element
        index = 1
    elseif val>=grid[n] #GEQ biggest element
        index = n
    else
        index_upper = findfirst(z->z>val, grid)
        index_lower = index_upper - 1
        val_upper, val_lower = grid[index_upper], grid[index_lower] #values
        index = index_lower + (val - val_lower)  / (val_upper - val_lower) #weighted average
    end
    index #return
end

###function to reshape flat vector of model Parameters into more manageable groups
function Reshape_param(prim::Primitives, guess::Vector{Float64})
    @unpack J, nν = prim #get number of occupations cq, race, etc
    poly_exp = 3 #degree of experience polynomial
    index = 1 #keep track of where we are in the vector

    #######first up: wage process parameters ########
    #non-teaching
    γ_jw, δ_j, σ_η = Any[], Any[], Any[]
    for i = 1:J-1 #loop over occupations
        nparam_wage = (1 + poly_exp + 1 + 1 + 1 + 1)-1 #intercept, exp poly, cq, sex, race, major
        push!(γ_jw, guess[index:(index + nparam_wage)]) #wage parameters for occupation j
        push!(δ_j, guess[index + nparam_wage + 1] #ability effects in ocupation j
        push!(σ_η, guess[index + nparam_wage + 2) #spread of wage shocks in occupation j
        index+=(nparam_wage + 3) #onto the next set of parameters!
    end

    #teaching
    nparam_wage_teach = (1 + poly_exp + 1)-1 #intercept, teacher experience, MA
    push!(γ_jw, guess[index:(index + nwage_params_teach)]) #teacher wage parameters
    push!(δ_j, 0.0) #no return to ability for teachers
    push!(σ_η, guess[index + nwage_params_teach + 1])
    index+=(nwage_params_teach+2) #onto the next set of parameters!

    #######Next up: teacher VA production#####
    nparam_va = (1 + poly_exp + 1 + 1 + 1 + 1 + 1 + 1)-1 #int, exp, cq, sex, race, major, MA, license
    γ_0 = guess[index:(index+nparam_va)] #VA parameters
    δ_0 = guess[index + nparam_va + 1] #ability effects
    σ_ς = guess[index + nparam_va + 2] #spread of VA shocks
    σ_ξ = guess[index + nparam_va + 3] #spread of teaching ability
    index += nparam_va + 4 #update index

    ######Preferences#######
    α = guess[index] #wage preference
    λ = guess[index + 1] #preference for teaching ability
    index += 2

    ###flow utility from work
    γ_ju, κ_j, ν_j = Any[], Any[], Any[]
    for i = 1:J #loop over occupations, including home work
        nparam_util = (1 + 1)-1  #sex, race
        push!(γ_ju, guess[index:(index + nparam_util)]) #demographic preference parameters for occupation J
        push!(κ_j, guess[index + nparam_util + 1) #switching cost
        push!(ν_j, guess[(index + nparam_util + 2):(index + nparam_util + 2 + nν-1)) #unobserved heterogeneity states (can reduce later via Gauss quadrature)
        index+=(nparam_util + 2 + nv) #onto the next set of parameters!
    end

    ###costs of teaching masters
    nparam_ma = (1 + 1 + 1 + 1 + 1)-1 #intercpet, cq, sex, race, teaching major
    γ_ma = guess[index:(index + nparam_ma)]
    δ_ma = guess[index + nparam_ma + 1] #ability effect
    index += (nparam_ma + 2)

    ###cost of certification
    nparam_l = (1 + 1 + 1 + 1 + 1 + 1)-1 #intercpet, cq, sex, race, teaching major, ma
    γ_l = guess[index:(index + nparam_l)]
    δ_l = guess[index + nparam_l + 1] #ability effect
    index += (nparam_ma + 2)

    ###first-period major choice
    γ_m, δ_m, ρ_m = Any[], Any[], Any[]
    for i = 1:nm #loop over major choices
        nparam_m = 1 + 1 + 1 + 1 + 1 #intercpet, cq, sex, race
        push!(γ_m, guess[index:(index + nparam_m)]) #demopgrahic major effect
        push!(δ_m, guess[index + nparam_m + 1]) #ability effects
        push!(ρ_m, guess[index + nparam_m + 2]) #unobserved taste comopnent
        index += (nparam_m + 3)
    end

    ###offer probabilities: license and experience polynomial interacted with last-period choice
    num = (poly_exp + 1)*2
    μ = guess[index:num]

    ###later: stuff governing correlation between unobserved ability and unobserved tastes

    #group together
    params = [γ_jw, δ_j, σ_η, γ_0, δ_0, σ_ς, σ_ξ, α, λ, γ_ju, κ_j, ν_j, γ_ma, δ_ma, γ_l, δ_l, γ_m, δ_m, μ]
    params
end


###function that defines and flattens initial guess of model parameters
function Param_init(prim::Primitives)
    @unpack J, poly_exp = prim

    ####initialize###

    ###wage parameters
    γ_jw =[ones(1 + poly_exp + 1 + 1 + 1 + 1) for j = 1:J-1, ones(1 + poly_exp + 1)]
    δ_j = [1.0 for j = 1:J-1, 0.0]
    σ_η = [1.0 for j = 1:J]

    ###teacher production
    γ_0 = ones(1 + poly_exp + 1 + 1 + 1 + 1 + 1 + 1) #intercept, experience, cq, sex, race, major, masters, license
    δ_0 = 1.0
    σ_ς = 1.0
    σ_ξ = 1.0

    ###preferences
    α = 1.0
    λ = 1.0

    #flow utility
    γ_ju = ones(2)
    κ_j = [1.0 for j = 1:J]
    ν_j = [1.0 for j = 1:J]

    #cost of masters
    γ_ma = ones(1 + 1 + 1 + 1 + 1) #intercept, cq, sex, race, major
    δ_ma = 1.0

    #cost of licensure
    γ_l = ones(1 + 1 + 1 + 1 +1 + 1) #intercept, cq, sex, race, major, masters
    δ_l = 1.0

    #first-period major choice
    γ_m = [ones(1 + 1 + 1 + 1) for m = 1:m] #intercpet, cq, sex, race
    δ_m = [1.0 for m = 1:nm]
    ρ_m = [1.0 for m = 1:nm]

    #teaching offer parameters
    μ = ones(9)

    ###later: stuff governing correlation between unobserved ability and unobserved tastes


    ####flatten####
    params_flat = []

    #wage effects
    for j = 1:J
        params_flat = vcat(params_flat, γ_jw[j])
        params_flat = vcat(params_flat) [δ_j[j], σ_η[j]])
    end

    #teacher production and some preferneces
    params_flat = vcat(params_flat, γ_0)
    params_flat = vcat(params_flat, [δ_0, σ_ς, σ_ξ, α, λ])

    #occupation preferences
    for j = 1:J
        params_flat = vcat(params_flat, γ_ju[j])
        params_flat = vcat(params_flat) [k_j[j])
        params_flat = vcat(params_flat) [ν_j[j])
    end

    #preferences for masters, licenses, and majors
    params_flat = vcat(params_flat, γ_ma)
    params_flat = vcat(params_flat, δ_ma)
    params_flat = vcat(params_flat, γ_l)
    params_flat = vcat(params_flat, δ_l)
    for m = 1:nm
        params_flat = vcat(params_flat, γ_m[m])
        params_flat = vcat(params_flat, δ_m[m])
        params_flat = vcat(params_flat, ρ_m[m])
    end
    params_flat = vcat(params_flat, μ)
    params_flat #return
end

######################
