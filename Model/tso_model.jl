##function to initialze model primitives and parameters
function Initialize(guess::Array{Any,1})
    prim = Primitives() #initialize model primitives
    @unpack nEJ, nξ, nm, nMA, nl, ndt, n𝒥, T, J = prim
    @unpack cq_grid, gender_grid, race_grid, θ_grid, ν_grid = prim

    #collect some state variables for ease
    X_grid = [] #demographics
    for c in cq_grid, g in gender_grid, r in race_grid
        push!(X_grid, [c, g, r])
    end

    χ_grid = [] #unobserved heterogeneity
    for θ in θ_grid
        for ν1 in ν_grid, ν2 in ν_grid #need as many nests as their are occupations
            push!(χ_grid, [θ, ν1, ν2])
        end
    end

    e_grid = [] #experience
    for e1 = 0:1:T #loop for constructing e_grid
        for e2 = 0:1:T #for future reference: need as many nests as you have jobs (don't know how to automate better)
            temp = [e1, e2]
            if e1 + e2<=T #nix too much experience
                push!(e_grid, temp)
            end
        end
    end

    prim_grp = Primitives_collect(T=T, X_grid = X_grid, χ_grid = χ_grid, e_grid = e_grid)
    @unpack nX, nχ, ne = prim_grp

    #initialize model parametesr
    p = Reshape_param(prim, guess)
    param = Params(guess = p)

    #initialize model value functions
    v_coll = zeros(nX, nχ, nm) #add states for major options
    v_work_a = zeros(nX, nχ, nm, nMA, nl, ne, nξ, nEJ, ndt, T, 2) #add states for time and choice of major
    v_work_b = zeros(nX, nχ, nm, nMA, nl, ne, nξ, nEJ, ndt, T, 2) #same, but license
    v_work_d = zeros(nX, nχ, nm, nMA, nl, ne, nξ, nEJ, ndt, n𝒥, T, J+1) #same, but now ocupation (and home work
    res = Results(v_coll, v_work_a, v_work_b, v_work_d)
    prim, prim_grp, param, res #return all the stuff
end

##function that initializes and solvse model
function Solve_model(guess::Array{Any,1})
    prim, prim_grp, param, res = Initialize(guess) #initialize important stuff
    Backward_induct(prim, prim_grp, param, res) #backward induction protocol
    #prim, param, res #return important stuff
end

#backward induction protocol
function Backward_induct(prim::Primitives, prim_grp::Primitives_collect, param::Params, res::Results)
    @unpack T = prim
    #for i = 1:T #loop over time periods
    for i = 1:T #loop over time periods
        t = T - i + 1 #now backwards
        Bellman_d(prim, prim_grp, param, res, t) #solve phase-4 choices and compute value functions
        Bellman_b(prim, prim_grp, param, res, t) #solve phase-4 choices and compute value functions
        Bellman_a(prim, prim_grp, param, res, t) #solve phase-4 choices and compute value functions
    end
    Bellman_coll(prim, prim_grp, param, res) #run period-0 Bellman
end

#College-period Bellman
function Bellman_coll(prim::Primitives, prim_grp::Primitives_collect, param::Params, res::Results)
    @unpack Π, T, nm = prim #unpack state space sizes
    @unpack X_grid, χ_grid, nX, nχ = prim_grp #grids
    @unpack v_work_a = res
    γ_eul = Base.MathConstants.eulergamma  #Euler's constant

    #begin main loop over state space
    for i_x = 1:nX, i_χ = 1:nχ
        X, χ = X_grid[i_x], χ_grid[i_χ] #initialize
        Ω = [X, χ] #collect state space

        for i_m = 1:nm #loop over major selection. We doing anything with priors of teaching ability?
            val = γ_eul + Π * log(sum(exp.(v_work_a[i_x, i_χ, i_m, 1, 1, 1, 1, 1, 1, 1, :]))) #expected value
            val += util_major(prim, param, Ω, i_m) #add on non-pecuniary utility for major
            res.v_coll[i_x, i_χ, i_m] = val #update
        end
    end
end

#phase-1 bellman function
function Bellman_a(prim::Primitives, prim_grp::Primitives_collect, param::Params, res::Results, t::Int64)
    @unpack nm, nMA, nl, nξ, nEJ, ndt = prim #unpack state space sizes
    @unpack m_grid, MA_grid, l_grid, EJ_grid, ξ_grid, dt_grid = prim #grids
    @unpack X_grid, χ_grid, nX, nχ, ne, e_grid = prim_grp
    @unpack v_work_b = res
    γ_eul = Base.MathConstants.eulergamma  #Euler's constant

    #begin main loop over state space
    for i_x = 1:nX, i_χ = 1:nχ, i_m = 1:nm, i_MA = 1:nMA, i_l = 1:nl, i_e = 1:ne, i_ξ = 1:nξ, i_EJ = 1:nEJ, i_d = 1:ndt

        #initialize state space and collect
        X, χ, m, MA, l, e, ξ, EJ, dt = X_grid[i_x], χ_grid[i_χ], m_grid[i_m], MA_grid[i_MA],
        l_grid[i_l], e_grid[i_e], ξ_grid[i_ξ], EJ_grid[i_EJ], dt_grid[i_d]
        Ω = [X, χ, m, MA, l, e, ξ, EJ, dt] #collect state space

        #check for inadmissible experience state for speed
        if sum(e) >= t #total experience must be less than t
            continue #skip
        end

        #nwo check whether we currently have a masters
        if MA == 1 #no choice to be made
            val = γ_eul + log(sum(exp.(v_work_b[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, :]))) #expected value

            if l == 1 #already have a license
                val = v_work_b[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, 1]
            end

            res.v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, :] .= val
        elseif MA == 0 #now we're making a choice
            for i_s = 1:2 #1 = no masters. 2 = masters
                val = γ_eul + log(sum(exp.(v_work_b[i_x, i_χ, i_m, i_s, i_l, i_e, i_ξ, i_EJ, i_d, t, :]))) #expected value

                if l == 1
                    val = v_work_b[i_x, i_χ, i_m, i_s, i_l, i_e, i_ξ, i_EJ, i_d, t, 1]
                end

                cost = cost_MA(prim, param, Ω) * (i_s-1) #cost of obtaining masters
                res.v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, i_s] = val - cost
            end
        end
    end
end

#phase-2 bellman function. Looks a lot like phase 1
function Bellman_b(prim::Primitives, prim_grp::Primitives_collect, param::Params, res::Results, t::Int64)
    @unpack J, nm, nMA, nl, nξ, nEJ, ndt = prim #unpack state space sizes
    @unpack m_grid, MA_grid, l_grid, EJ_grid, ξ_grid, dt_grid = prim #grids
    @unpack X_grid, χ_grid, nX, nχ, ne, e_grid = prim_grp
    @unpack v_work_d = res
    γ_eul = Base.MathConstants.eulergamma  #Euler's constant

    #begin main loop over state space
    for i_x = 1:nX, i_χ = 1:nχ, i_m = 1:nm, i_MA = 1:nMA, i_l = 1:nl, i_e = 1:ne, i_ξ = 1:nξ, i_EJ = 1:nEJ, i_d = 1:ndt

        #initialize state space and collect
        X, χ, m, MA, l, e, ξ, EJ, dt = X_grid[i_x], χ_grid[i_χ], m_grid[i_m], MA_grid[i_MA],
        l_grid[i_l], e_grid[i_e], ξ_grid[i_ξ], EJ_grid[i_EJ], dt_grid[i_d]
        Ω = [X, χ, m, MA, l, e, ξ, EJ, dt] #collect state space

        #check for inadmissible experience state for speed
        if sum(e) >= t #total experience must be less than t
            continue #skip
        end

        #nwo check whether we currently have a license
        if l == 1 #no choice to be made
            μ_j = μ(prim, param, Ω, 1) #probability of teaching offer given parameters, state space, and 1 for license
            val = γ_eul + μ_j * log(sum(exp.(v_work_d[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, 2, t, 1:J+1]))) #expected value, offer
            val += (1-μ_j) * log(sum(exp.(v_work_d[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, 1, t, 1:J]))) #expected value, no offer
            res.v_work_b[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, :] .= val
        elseif l == 0 #now we're making a choice
            for i_s = 1:2 #1 = no license. 2 = license
                μ_j = μ(prim, param, Ω, i_s - 1) #probability of teaching offer given parameters, state space, and 1 for license
                val = γ_eul + μ_j * log(sum(exp.(v_work_d[i_x, i_χ, i_m, i_MA, i_s, i_e, i_ξ, i_EJ, i_d, 2, t, 1:J+1]))) #expected value, offer
                val += (1-μ_j) * log(sum(exp.(v_work_d[i_x, i_χ, i_m, i_MA, i_s, i_e, i_ξ, i_EJ, i_d, 1, t, 1:J]))) #expected value, no offer
                cost = cost_license(prim, param, Ω) * (i_s-1) #cost of obtaining masters
                res.v_work_b[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, i_s] = val - cost
            end
        end
    end
end

#phase-4 bellman function
function Bellman_d(prim::Primitives, prim_grp::Primitives_collect, param::Params, res::Results, t::Int64)
    @unpack β, J, T, nm, nMA, nl, nξ, nEJ, ndt, n𝒥 = prim #unpack state space sizes
    @unpack m_grid, MA_grid, l_grid, EJ_grid, ξ_grid, dt_grid, 𝒥_grid = prim #grids
    @unpack X_grid, χ_grid, nX, nχ, ne, e_grid = prim_grp
    @unpack v_work_a = res
    γ_eul = Base.MathConstants.eulergamma  #Euler's constant

    #first: check if in terminal period
    if t == T #in terminal period
        #loop over state space
        for i_x = 1:nX, i_χ = 1:nχ, i_m = 1:nm, i_MA = 1:nMA, i_l = 1:nl, i_e = 1:ne, i_ξ = 1:nξ, i_EJ = 1:nEJ, i_d = 1:ndt, i_𝒥 = 1:n𝒥

            #initialize state space and collect
            X, χ, m, MA, l, e, ξ, EJ, dt, 𝒥 = X_grid[i_x], χ_grid[i_χ], m_grid[i_m], MA_grid[i_MA],
            l_grid[i_l], e_grid[i_e], ξ_grid[i_ξ], EJ_grid[i_EJ], dt_grid[i_d], 𝒥_grid[i_𝒥]
            Ω = [X, χ, m, MA, l, e, ξ, EJ, dt, 𝒥] #collect state space

            #check for inadmissible experience state for speed
            if sum(e) >= t #total experience must be less than t
                continue #skip
            end

            #since availability of teaching offer doesn't impact flow utility from other professions, we can ignore
            #the 𝒥 state and just loop over occupation choice
            for i_j = 1:J #loop over occupation choices
                val = util(prim, param, Ω, i_j, J)
                res.v_work_d[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, i_𝒥, t, i_j+1] = val #update
            end
            res.v_work_d[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, i_𝒥, t, 1] = 0.0 #update home work option\
        end
    elseif t!=T #not in terminal period

        #loop over state space
        for i_x = 1:nX, i_χ = 1:nχ, i_m = 1:nm, i_MA = 1:nMA, i_l = 1:nl, i_e = 1:ne, i_ξ = 1:nξ, i_EJ = 1:nEJ, i_d = 1:ndt, i_𝒥 = 1:n𝒥
            #initialize state space and collect
            X, χ, m, MA, l, e, ξ, EJ, dt, 𝒥 = X_grid[i_x], χ_grid[i_χ], m_grid[i_m], MA_grid[i_MA],
            l_grid[i_l], e_grid[i_e], ξ_grid[i_ξ], EJ_grid[i_EJ], dt_grid[i_d], 𝒥_grid[i_𝒥]
            Ω = [X, χ, m, MA, l, e, ξ, EJ, dt, 𝒥] #collect state space

            #check for inadmissible experience state for speed
            if sum(e) >= t #total experience must be less than t
                continue #skip
            end

            #since availability of teaching offer doesn't impact flow utility from other professions, we can ignore
            #the 𝒥 state and just loop over occupation choice
            for i_j = 1:J #loop over occupation choices
                val = util(prim, param, Ω, i_j, J) #flow utility from occupation choice
                e_next = e #next-period experience given job choice
                e_next[i_j] += 1 #add one to occupation-specific experience.
                i_e_next = findfirst(x->x==e_next, e_grid)

                #construct continuation value
                val_cont = β * (γ_eul + log(sum(exp.(v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e_next, i_ξ, i_EJ, i_j+1, t+1, :])))) #add continuation value
                if EJ == 0 && i_j == J #first time teaching
                    val_cont = 0
                    for i_ξ_next = 1:nξ #loop over potential teacher quality realizations
                        val_cont += (β/nξ) * (γ_eul + log(sum(exp.(v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e_next, i_ξ_next, 2, i_j+1, t+1, :])))) #add continuation value
                    end
                end

                if MA == 1 #already have a masters
                    val_cont = β * v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e_next, i_ξ, i_EJ, i_j+1, t+1, 1]
                    if EJ == 0 && i_j == J #first time teaching
                        val_cont = 0
                        for i_ξ_next = 1:nξ #loop over potential teacher quality realizations
                            val_cont +=  (β/nξ) * v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e_next, i_ξ_next, i_EJ, i_j+1, t+1, 1] #add continuation value
                        end
                    end
                end
                val += val_cont
                res.v_work_d[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, i_𝒥, t, i_j+1] = val #update
            end

            val_nwork = β * (γ_eul + log(sum(exp.(v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, 1, t+1, :]))))
            if MA == 1
                val_nwork = β * v_work_a[i_x, i_χ, i_m, i_MA, i_l, i_e_next, i_ξ, i_EJ, 1, t+1, 1]
            end
            res.v_work_d[i_x, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, i_𝒥, t, 1] = val_nwork #update home work option
        end
    end
end
########
