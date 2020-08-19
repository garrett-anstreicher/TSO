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
function Reshape_param(prim::Primitives, guess::Array{Any,1})
    @unpack J, nν, poly_exp, nm = prim #get number of occupations cq, race, etc
    index = 1 #keep track of where we are in the vector

    #######first up: wage process parameters ########
    #non-teaching
    γ_jw, δ_j, σ_η = Any[], Any[], Any[]
    for i = 1:J-1 #loop over occupations
        nparam_wage = (1 + poly_exp + 1 + 1 + 1 + 1)-1 #intercept, exp poly, cq, sex, race, major
        push!(γ_jw, guess[index:(index + nparam_wage)]) #wage parameters for occupation j
        push!(δ_j, guess[index + nparam_wage + 1]) #ability effects in ocupation j
        push!(σ_η, guess[index + nparam_wage + 2]) #spread of wage shocks in occupation j
        index+=(nparam_wage + 3) #onto the next set of parameters!
    end


    #teaching
    nparam_wage_teach = (1 + poly_exp + 1)-1 #intercept, teacher experience, MA
    push!(γ_jw, guess[index:(index + nparam_wage_teach)]) #teacher wage parameters
    push!(δ_j, guess[index + nparam_wage_teach + 1]) #no return to ability for teachers
    push!(σ_η, guess[index + nparam_wage_teach + 2])
    index+=(nparam_wage_teach+3) #onto the next set of parameters!

    #######Next up: teacher VA production#####
    nparam_va = (1 + poly_exp + 1 + 1 + 1 + 1 + 1 + 1)-1 #int, exp, cq, sex, race, major, MA, license
    γ_0 = guess[index:(index+nparam_va)] #VA parameters
    δ_0 = guess[index + nparam_va + 1] #ability effects
    σ_ς = guess[index + nparam_va + 2] #spread of VA shocks
    σ_ξ = guess[index + nparam_va + 3] #spread of teaching ability
    index += (nparam_va + 4) #update index


    ######Preferences#######
    α = guess[index] #wage preference
    λ = guess[index + 1] #preference for teaching ability
    ν = guess[(index + 2):(index + 2 + nν - 1)]
    index += (2 + nν)

    ###flow utility from work
    γ_ju, κ_j = Any[], Any[]
    for i = 1:J #loop over occupations, including home work
        nparam_util = (1 + 1)-1  #sex, race
        push!(γ_ju, guess[index:(index + nparam_util)]) #demographic preference parameters for occupation J
        push!(κ_j, guess[index + nparam_util + 1]) #switching cost
        index+=(nparam_util + 2) #onto the next set of parameters!
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
    index += (nparam_l + 2)

    ###first-period major choice
    γ_m, δ_m, ρ_m = Any[], Any[], Any[]
    for i = 1:nm #loop over major choices
        nparam_m = (1 + 1 + 1 + 1)-1  #intercpet, cq, sex, race
        push!(γ_m, guess[index:(index + nparam_m)]) #demopgrahic major effect
        push!(δ_m, guess[index + nparam_m + 1]) #ability effects
        push!(ρ_m, guess[index + nparam_m + 2]) #unobserved taste comopnent
        index += (nparam_m + 3)
    end

    ###offer probabilities: license and experience polynomial interacted with last-period choice
    num = (poly_exp + 1)*2
    μ = guess[index:(index + num)]

    ###later: stuff governing correlation between unobserved ability and unobserved tastes

    #group together
    params = [γ_jw, δ_j, σ_η, γ_0, δ_0, σ_ς, σ_ξ, α, λ, ν, γ_ju, κ_j, γ_ma, δ_ma, γ_l, δ_l, γ_m, δ_m, ρ_m, μ]
    params
end


###function that defines and flattens initial guess of model parameters.
#Accepts as arguments number of occupations and legnth of experience polynomial
function Param_init(J::Int64, poly_exp::Int64, nm::Int64, nν::Int64)
    ###wage parameters
    γ_jw =[ones(1 + poly_exp + 1 + 1 + 1 + 1) for j = 1:J-1]./1000
    push!(γ_jw, ones(1 + poly_exp + 1)./100) #add on teacher parameters
    δ_j = [0.0 for j = 1:J-1]./100
    push!(δ_j, 0.0) #add teacher ability effect (0)
    σ_η = [1.0 for j = 1:J]

    ###teacher production
    γ_0 = ones(1 + poly_exp + 1 + 1 + 1 + 1 + 1 + 1)./100 #intercept, experience, cq, sex, race, major, masters, license
    δ_0 = 0.0
    σ_ς = 1.0
    σ_ξ = 1.0

    ###preferences
    α = 0.01
    λ = 0.01
    ν = [0.0 for j = 1:nν] #unobserved taste levels for occupations

    #flow utility
    γ_ju = [zeros(2) for j = 1:J]
    κ_j = [1.0 for j = 1:J]

    #cost of masters
    γ_ma = zeros(1 + 1 + 1 + 1 + 1) #intercept, cq, sex, race, major
    δ_ma = 1.0

    #cost of licensure
    γ_l = zeros(1 + 1 + 1 + 1 +1 + 1) #intercept, cq, sex, race, major, masters
    δ_l = 1.0

    #first-period major choice
    γ_m = [zeros(1 + 1 + 1 + 1) for m = 1:nm] #intercpet, cq, sex, race
    δ_m = [1.0 for m = 1:nm]
    ρ_m = [0.0 for m = 1:nm]

    #teaching offer parameters
    μ = ones(9)

    ###later: stuff governing correlation between unobserved ability and unobserved tastes

    ####flatten####
    params_flat = []

    #wage effects
    for j = 1:J
        params_flat = vcat(params_flat, γ_jw[j])
        params_flat = vcat(params_flat, [δ_j[j], σ_η[j]])
    end

    #teacher production and some preferneces
    params_flat = vcat(params_flat, γ_0)
    params_flat = vcat(params_flat, [δ_0, σ_ς, σ_ξ, α, λ])
    params_flat = vcat(params_flat, ν)

    #occupation preferences
    for j = 1:J
        params_flat = vcat(params_flat, γ_ju[j])
        params_flat = vcat(params_flat, κ_j[j])
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
