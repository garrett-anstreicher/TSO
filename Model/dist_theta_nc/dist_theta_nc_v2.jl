using Parameters, Distributions

#parametesr that Chris estimated
@with_kw struct Primitives
    θ::Vector{Float64} = [-2.440, -1.825, -1.342, -0.853, -0.388, 0.010, 0.343, 0.898] #values of θ
    ϵ_0::Float64 = 0.05

    #probabilities of θ by demographic group
    prob_θ_wm::Vector{Float64} = [0.016, 0.085, 0.154, 0.214, 0.216, 0.134, 0.117, 0.064]
    prob_θ_ww::Vector{Float64} = [0.025, 0.117, 0.197, 0.255, 0.190, 0.101, 0.077, 0.038]
    prob_θ_mm::Vector{Float64} = [0.117, 0.202, 0.190, 0.208, 0.086, 0.101, 0.049, 0.046]
    prob_θ_mw::Vector{Float64} = [0.167, 0.240, 0.210, 0.183, 0.091, 0.060, 0.035, 0.015]
    prob_θ::Vector{Vector{Float64}} = [prob_θ_wm, prob_θ_ww, prob_θ_mm, prob_θ_mw] #all together

    #parameters for SAT Math/Verbal (note: does not vary by demographics)
    param_m::Vector{Float64} = [0.0, 1.0, 0.644, 0.450, 0.382, 0.022, -0.003, 0.005]
    param_v::Vector{Float64} = [-0.261, 1.030, 0.592, 0.211, 0.347, 0.151, 0.023, 0.005]

    #parameters for CQ (do depend on demographics)
    param_cq_wm::Vector{Float64} = [0.669, 0.610, 0.496, 0.448, -0.075, 0.063, -0.002, 0.016] #white men
    param_cq_ww::Vector{Float64} = [1.447, 0.544, 0.693, 0.343, -0.129, 0.086, -0.061, 0.011] #white women
    param_cq_mm::Vector{Float64} = [0.663, 0.703, 0.579, 0.488, -0.117, 0.017, 0.022, 0.021] #minority men
    param_cq_mw::Vector{Float64} = [0.660, 0.645, 0.637, 0.497, -0.091, 0.038, 0.006, 0.015] #minority women
    param_cq::Vector{Vector{Float64}} = [param_cq_wm, param_cq_ww, param_cq_mm, param_cq_mw] #all together
end


#main function for SAT math. Inputs normalized test score/quality, score we're using, and agent type. outputs theta probabilities
function post_dist(prim::Primitives, scores::Vector{Float64}, type::Int64, dem::Int64)
    @unpack θ, prob_θ, param_m, param_v, param_cq, ϵ_0 = prim
    nθ = length(θ)
    p_θs = prob_θ[dem] #probabilities of θ depending on demographic group
    prob_θ_post = copy(p_θs) #starting off equal to unconditional probabilities

    #all parameters we might use
    params_all = [copy(param_m), copy(param_v), copy(param_cq[dem])]

    for s = 1:3 #loop over scores
        if scores[s] == 999.0 #value indicating "missingness," i.e. skip
            continue
        end

        b0 = params_all[s][1] #constant
        b1 = params_all[s][2] #ability slope
        σ = params_all[s][3] #sigma
        x_poly = params_all[s][4:8] #density polynomial

        for i = 1:nθ #loop over possible θs provided we have a valid score to use
            x = (scores[s] - (b0 + b1 * θ[i]))/σ #compute measurement error, normalized by σ
            x_vec = [x^0, x^1, x^2, x^3, x^4] #polynomial

            p_score_θ = (((sum(x_vec.*x_poly))^2) * exp((-x^2)/2) + ϵ_0 * pdf(Normal(), x)) / σ
            prob_θ_post[i] *= p_score_θ
        end
    end


    prob_θ_post = prob_θ_post ./ sum(prob_θ_post)

    #return deliverable
    return prob_θ_post
end


#testing behavior with college quality
prim = Primitives()
temp = post_dist(prim, [999.0, 999.0, 0.0], 3, 1)
temp = post_dist(prim, [0.0, 0.0, 1.0], 3, 1)


#higher CQ score results in higher probabilities of large θs
temp = post_dist(prim, 1.0, 3, 1)

#but lowr probabilities of high θs for minority women (for example) due to different measuremnet error distrbution
temp = post_dist(prim, 1.0, 3, 4)

#converse also works.
temp = post_dist(prim, -1.0, 3, 1)
temp = post_dist(prim, -1.0, 3, 4)


####
