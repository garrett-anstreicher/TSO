####struct of model primitives
@with_kw struct Primitives
    delta::Float64 = 0.95 #discount factor

    #age grid
    a_grid::Array{Int64,1} = collect(45:1:64) #grid
    na::Int64 = length(a_grid) #number of elements

    #husband FE grid
    h_grid::Array{Float64,1} = collect(-1.8:0.2:1.8)
    nh::Int64 = length(h_grid)

    #experience grid
    x_grid::Array{Int64,1} = collect(0:1:33) #grid
    nx::Int64 = length(x_grid) #length

    #schooling grid
    s_grid::Array{Int64,1} = collect(5:1:17) #grid
    ns::Int64 = length(s_grid) #length

    #lfp status grid
    p_grid::Array{Int64,1} = collect(0:1:1)
    np::Int64 = length(p_grid)

    ##parameters of husband wage process; depends on wife age and schooling
    beta_0::Float64 = 6.340498 #constant
    beta_1::Float64 = 0.061502 # #school
    beta_2::Float64 = 0.109508 #wife age
    beta_3::Float64 = -0.001 #wife age ^2
    beta_4::Float64 = 0.984249 #fe
end

#struct of model estimands
mutable struct Estimands
    alpha_1::Float64 #utility term
    alpha_2::Float64 #utility term
    alpha_3::Float64 #hassle cost
    gamma_0::Float64 #wage term
    gamma_1::Float64 #wage term
    gamma_2::Float64 #wage term
    gamma_3::Float64 #wage term
    sig_eps::Float64 #wage shock
    sig_eta::Float64 #wage error
end

#struct of model results
mutable struct Results
    e_maxes::Array{Float64,5} #e_max values
    cutoffs::Array{Float64,5} #cutoff rules
end

#function to initialize model primitives and feed in initial guess of estimands
function Initialize(guess::Array{Float64,1})
    prim = Primitives() #initialize model primitives
    @unpack na, nh, nx, ns, np = prim
    e_maxes = zeros(na, nh, nx, ns, np)
    cutoffs = zeros(na, nh, nx, ns, np)
    res = Results(e_maxes, cutoffs) #initialize results struct

    #initialize estimand struct
    estim = Estimands(guess[1], guess[2], guess[3], guess[4], guess[5],
    guess[6], guess[7], guess[8], guess[9])

    prim, res, estim #return stuff
end

#function to quickly spit out husband income given observed states and husband FE
function Hus_wage(prim::Primitives, a::Int64, s::Int64, h::Float64)
    @unpack beta_0, beta_1, beta_2, beta_3, beta_4 = prim #unpack
    h_wage = exp(beta_0 + beta_1*s + beta_2*a + beta_3*a^2 + beta_4*h) #deterministic husband wage
    h_wage #return
end

#progressive income taxation function
function Prog_tax(inc::Float64)
    earnings = 0.9*inc
    if inc>50000
        earnings = 45000 + 0.8*(inc-50000)
    end
    earnings #return progressive after-tax earnings
end

#function to solve non-terminal e_max and cutoff rules
function Bellman(prim::Primitives, res::Results, estim::Estimands, a::Int64, cfact::Int64)
    @unpack na, np, nh, nx, ns, p_grid, a_grid, h_grid, x_grid, s_grid, delta = prim #unpack primitives
    @unpack alpha_1, alpha_2, alpha_3, gamma_0, gamma_1, gamma_2, gamma_3, sig_eps = estim #unpack estimands
    dist = Normal() #initialize normal distribution
    a_index = a-44 #get index of age array that current value corresponds to
    root = sqrt(sig_eps) #standard error of wage shock

    for i1 = 1:nh, i2 = 1:nx, i3 = 1:ns, i4 = 1:np
        h, x, s, p = h_grid[i1], x_grid[i2], s_grid[i3], p_grid[i4]#state space
        y = Hus_wage(prim, a, s, h) #fetch husband income
        temp = exp(gamma_0 + gamma_1*s + gamma_2*x + gamma_3*(x^2)) #for exposition
        n_emax_0, n_emax_1 = 0, 0 #preallocate

        if a_index!=na #fill in next-period utilities if not in the terminal period
            n_emax_0 = delta*res.e_maxes[a_index+1, i1, i2, i3, 1] #no change to x
            n_emax_1 = delta*res.e_maxes[a_index+1, i1, min(i2+1, nx), i3, 2]  #keep at max level
        end
        emax_diff = n_emax_0 - n_emax_1
        cutoff, emax_1, emax_2 = 0, 0, 0 #preallocation

        multiplier = -1 #thing to store that keeps track o fhow we add in alpha_3 into cutoffs
        if p==1
            multiplier = 1
        end
        cutoff = max(-Inf,log(max((alpha_1 + alpha_2*y + multiplier*alpha_3 + emax_diff),0))) - log(temp)

        #begin counterfactual edits to cutoff rule and continuation values
        if cfact == 2 #first counterfactual
            y = 0.9*y #update husband wage
            cutoff = max(-Inf,log(max((alpha_1 + alpha_2*y + multiplier*alpha_3 + emax_diff),0))) + log(1.0/0.9) - log(temp)
            temp = 0.9*temp #update this for use in e_max term
        elseif cfact == 3 #this one takes a bit more work
            y_new = 0 #preallocation
            if y>50000 #husband makes more than 50k alone
                 y_new = 45000 + 0.8*(y - 50000) #this one's not so bad
                 cutoff = max(-Inf,log(max((alpha_1 + alpha_2*y + multiplier*alpha_3 + emax_diff),0))) + log(1.0/0.8) - log(temp)
                 temp = 0.8*temp #update this for use in e_max term
            elseif y<=50000 #husband makes less than 50k; a bit trickier
                y_new = 0.9 * y #this part's easier, at least

                #final finesse: updating cutoff rule. perhaps the subtlest part. do by minimizing differene between alterate utilities
                thresh = y_new + alpha_1 + alpha_2*y + multiplier*alpha_3 + emax_diff
                obj(shock) = abs(Prog_tax(y_new + temp*exp(shock)) - thresh)
                cutoff = optimize(obj, -2.0, 2.0).minimizer

                temp = 0.9*temp #update expected temp value
                if temp+y>50000
                    temp = 0.9*(50000-y) + 0.8*(temp - 50000 + y)
                end
            end
            y = y_new #update value for husband income to use in emax part
        end
        res.cutoffs[a_index, i1, i2, i3, i4] = cutoff #update cutoff vector

        if p==0 #update emax
            emax_1 = (y + alpha_3 + n_emax_1)*(1-cdf(dist, cutoff/root)) + temp*exp(0.5*sig_eps)*(1-cdf(dist, (cutoff-sig_eps)/root)) #first part of emax
            emax_2 = (y*(1+alpha_2)+alpha_1 + n_emax_0)*cdf(dist, cutoff/root) #second part of emax
        elseif p == 1 #the same
            emax_1 = (y + n_emax_1)*(1-cdf(dist, cutoff/root)) + temp*exp(0.5*sig_eps)*(1-cdf(dist, (cutoff-sig_eps)/root)) #first part of emax
            emax_2 = (y*(1+alpha_2)+alpha_1 + alpha_3 + n_emax_0)*cdf(dist, cutoff/root) #second part of emax
        end
        res.e_maxes[a_index, i1, i2, i3, i4] = emax_1 + emax_2 #update emax vector
    end
end

#function that solves the e_max values and cutoff rules given primitives and estimated parameters via backward induction
function Backward_induct(prim::Primitives, res::Results, estim::Estimands; cfact::Int64 = 0)
    @unpack na = prim

    for i = 1:na #begin backward induction loop
        a = 65 - i #now looping from 64 to 45
        Bellman(prim, res, estim, a, cfact) #update cutoff rules and emaxes
    end
end

#function to compute likelihood given cutoff rules, estimands and data
function Likelihood(prim::Primitives, res::Results, estim::Estimands, data::Array{Float64,2})
    @unpack alpha_1, alpha_2, gamma_0, gamma_1, gamma_2, gamma_3, sig_eps, sig_eta = estim #unpack estimands
    @unpack h_grid = prim
    likelihood = 0.0 #start off likelihood value
    dist = Normal() #initialize standard normal distribution
    root_eps, root_eta = sqrt(sig_eps), sqrt(sig_eta)
    root_u = sqrt(sig_eps + sig_eta)
    rho = root_eps/root_u #various terms

    for i = 1:15060 #loop through rows of data
        prob = 0 #preallocation of likelihood to add
        row = data[i,:] #isolate row
        a_index = Int64(row[2]-44) #age index
        x_index = Int64(row[4]+1) #experience index
        s_index = Int64(row[6]-4) #school index
        h_index = findall(x->x==row[10], h_grid)[1] #husband FE index
        p_index = Int64(row[7]+1) #lfp0

        if row[2] != 45 #use last row in this case
            p_index = Int64(data[i-1, 3]+1) #p index; use lfp from pervious row
        end
        cutoff = res.cutoffs[a_index, h_index, x_index, s_index, p_index] #get cutoff value

        if row[3] == 0 #update lfp; not working
            prob = cdf(dist, cutoff/root_eps)
        elseif row[3] == 1 #update lfp; working
            #step 1: compute total error
            x, s = row[4], row[6]
            w_predict = gamma_0 + gamma_1*s + gamma_2*x + gamma_3*x^2 #preeicted wage
            u = log(row[5]) - w_predict

            frac = (cutoff-(rho^2)*u)/(root_eps*sqrt(1-rho^2)) #assembly of probability
            prob_1 = 1-cdf(dist, frac) #assembly of probability
            prob_2 = (1/root_u)*pdf(dist, (u/root_u)) #assembly of probability
            prob = prob_1*prob_2 #formation of probability
        end
        likelihood+=log(prob) #update likelihood value, logged to avoid shooting to zero
    end
    likelihood #return likelihood value
end

##objective function that we wish to minimize in estimation
function Objective_function(guess::Array{Float64,1})
    prim, res, estim = Initialize(guess) #initialize vectors used in computation
    likelihood = 0

    try
        Backward_induct(prim, res, estim) #backward induction
        likelihood = Likelihood(prim, res, estim, data) #likelihood
    catch err #if guesses produce domain errors, just make the likelihood really bad
        if isa(err, DomainError) #in case something goes wrong
            likelihood = -Inf #punishment!!
        end
    end
    -likelihood #return negative likelihood to make sure we're minimizing instead
end

#function to simulate data used to assess model fit
function Simulate(prim::Primitives, res::Results, estim::Estimands, data::Array{Float64,2}; cfact::Int64 = 0)
    @unpack h_grid, s_grid, x_grid = prim #get variance of wage shock
    @unpack sig_eps, sig_eta, gamma_0, gamma_1, gamma_2, gamma_3 = estim
    dist = Normal(0, sqrt(sig_eps)) #distribution of wage shock to draw things from
    dist_eta = Normal(0, sqrt(sig_eta))
    data_simul = Any[] #preallocate

    #begin looping through rows of data
    for i = 1:15060 #loop through rows of data
        row = data[i,:] #row of data
        if row[2] == 45 #new family to simulate!
            for j = 1:20   #simulate family 20 times
                fam_data = zeros(20,11) #simulated data array to fill -- 20 years for each couple
                fam_data[1,4] = row[4]  #starting x
                fam_data[:,2] = collect(45:1:64) #age

                #simulation begins
                for k = 1:20 #populate rows of simulated data

                    #fill in some of the data
                    fam_data[k,1] = row[1] #id
                    fam_data[k,6] = row[6] #edu
                    fam_data[k,7] = row[7] #starting lfp
                    fam_data[k,9] = row[9] #husband FE
                    fam_data[k,10] = row[10] #rounded husband FE
                    fam_data[k, 8] = Hus_wage(prim, Int64(fam_data[k,2]), Int64(fam_data[k,6]), fam_data[k,10])

                    #get indices of state space
                    a_index = Int64(fam_data[k,2]-44) #age index
                    x_index = Int64(fam_data[k,4]+1) #experience index
                    s_index = Int64(fam_data[k,6]-4) #school index
                    h_index = findall(x->x==fam_data[k,10], h_grid)[1] #husband FE index
                    s, x = fam_data[k,6], fam_data[k,4] #for later use

                    #last period lfp index
                    p_index = 0
                    if k==1 #use lfp0 in this case
                        p_index = fam_data[k,7]+1
                    elseif k>1
                        p_index = fam_data[k-1, 3]+1 #use lfp from pervious row
                    end
                    p_index = Int64(p_index)

                    cutoff = res.cutoffs[a_index, h_index, x_index, s_index, p_index] #get cutoff value
                    shock = rand(dist) #draw a shock

                    if shock>=cutoff #agent chooses to work
                        fam_data[k,3] = 1 #fill in lfp
                        error = rand(dist_eta) #draw some measurement error
                        fam_data[k,5] =  exp(gamma_0 + gamma_1*s + gamma_2*x + gamma_3*(x^2) + shock + error) #fill in wage
                        fam_data[k,11] =  exp(gamma_0 + gamma_1*s + gamma_2*x + gamma_3*(x^2) + shock) #fill in true wage

                        if k!=20 #if not terminal period, write in next period x value
                            if fam_data[k,4]==33 #max value
                                fam_data[k+1,4] = fam_data[k,4]
                            elseif fam_data[k,4]!=33 #not max value
                                fam_data[k+1,4] = fam_data[k,4]+1
                            end
                        end

                    elseif shock<cutoff #agent does not
                        fam_data[k,3] = 0 #fill in lfp
                        fam_data[k,5] =  -9 #fill in wage
                        fam_data[k,11] =  -9 #fill in true wage (also -9 if not working!)

                        if k!=20 #if not terminal period, write in next period x value
                            fam_data[k+1,4] = fam_data[k,4]
                        end
                    end
                end
                push!(data_simul, fam_data) #add to simulated ata
            end
        elseif row[2]!=45 #skip if not a new family
            continue
        end
    end
    data_simul = vcat(data_simul...) #stick all together into one array
    data_simul #return simulated data
end
#########################################################################
