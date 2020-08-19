#########Simulation routine
function Simulate(prim::Primitives, prim_grp::Primitives_collect, param::Params, res::Results; nsim::Int64=100000)
    @unpack β, J, T, nm, nMA, nl, nξ, nEJ, ndt, n𝒥 = prim #unpack state space sizes
    @unpack m_grid, MA_grid, l_grid, EJ_grid, ξ_grid, dt_grid, 𝒥_grid = prim #grids
    @unpack X_grid, χ_grid, nX, nχ, ne, e_grid = prim_grp
    @unpack v_work_a, v_work_b, v_work_d, v_coll = res
    @unpack σ_η, σ_ς = param
    data_simul = Any[] #preallocate simulated data
    dist_gumbel = Gumbel(0,1)



    for i = 1:nsim #number of simulations
        i_X, i_χ, i_ξ = rand(1:nX), rand(1:nχ), rand(1:nξ) #to start: randomly draw demographics and unobserved heterogeneity
        i_m = findmax(v_coll[i_X, i_χ, :] .+ rand(dist_gumbel, nm))[2] #major choice

        #things to record each period: demographics, major, MA, license, experience, teacher quality/VA, last/current occupation, time, wage
        cq, g, r, θ, m, ξ = X_grid[i_X][1], X_grid[i_X][2], X_grid[i_X][3], χ_grid[i_χ][1], m_grid[i_m], ξ_grid[i_ξ]
        MA, l = 0, 0 #preallocate
        e = zeros(J) #initialize no experience
        EJ, dt, j = 0, 0, 0 #preallocate

        #initialize other state indices that will be updated over time
        i_MA, i_l, i_e, i_EJ, i_d = 1, 1, 1, 1, 1

        for t = 1:T #loop over time periods
            #phase A choice
            if i_MA == 1
                MA_choice = findmax(v_work_a[i_X, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, :] .+ rand(dist_gumbel,2))[2]
                if MA_choice == 2 #individual gets an MA
                    i_MA, MA = 2, 1 #permanently update states
                end
            end

            #phase B choice
            if i_l == 1 #no license
                l_choice = findmax(v_work_b[i_X, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, t, :] .+ rand(dist_gumbel,2))[2]
                if MA_choice == 2 #individual gets a license
                    i_l, l = 2, 1 #permanently update states
                end
            end

            #initialzie other states and collect for ease
            X, χ, EJ, dt = X_grid[i_X], χ_grid[i_χ], EJ_grid[i_EJ], dt_grid[i_d]
            Ω = [X, χ, m, MA, l, e, ξ, EJ, dt] #collect state space

            #Phase C: draw of offer
            offer_prob = μ(prim, param, Ω, l)
            draw = rand() #random draw
            𝒥 = 0
            if draw<offer_prob #individual receives an offer
                𝒥 = 1 #update
            end

            #phase D choice
            j = findmax(v_work_d[i_X, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, 𝒥 + 1, t, :] + rand(dist_gumbel, J+1))[2] - 1
            if 𝒥 == 0 #no teaching offer
                j = findmax(v_work_d[i_X, i_χ, i_m, i_MA, i_l, i_e, i_ξ, i_EJ, i_d, 𝒥 + 1, t, 1:J] + rand(dist_gumbel, J))[2] - 1 #limited menu
            end

            if j>0 #update experience
                e[j]+=1
                i_e = findfirst(x->x==e, e_grid)
            end

            if j == J #update EJ
                EJ, i_EJ = 1, 2
            end

            #Phase E: realized wages and VA (if teaching)
            wage, value = 0,0
            if j>0 #wage
                wage = exp(w(prim, param, Ω, j) + rand(Normal(0, σ_η[j])))
            end

            if j == J #va
                value = va(prim, param, Ω) + rand(Normal(0, σ_ς))
            end

            line = [i t cq g r θ m MA l] #first pargt
            line = hcat(line, e') #add expoerience
            line = hcat(line, [ξ dt j wage value]) #add other stuff
            dt, i_d = j, j+1 #upate dt at end
            push!(data_simul, line) #write new line of data
        end
    end
    data_simul = vcat(data_simul...) #all together
    data_simul #return
end








####
