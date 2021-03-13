using Parameters, Distributions

#parametesr that Chris estimated
@with_kw struct Primitives
    θ::Vector{Float64} = [-2.046, -1.395, -0.903, -0.398, 0.079, 0.447, 0.750, 1.320] #values of θ
    ϵ_0::Float64 = 0.05

    #probabilities of θ by demographic group
    prob_θ_wm::Vector{Float64} = [0.016, 0.088, 0.163, 0.220, 0.220, 0.103, 0.124, 0.065]
    prob_θ_ww::Vector{Float64} = [0.024, 0.123, 0.204, 0.263, 0.184, 0.082, 0.081, 0.038]
    prob_θ_mm::Vector{Float64} = [0.113, 0.211, 0.198, 0.211, 0.072, 0.105, 0.042, 0.048]
    prob_θ_mw::Vector{Float64} = [0.162, 0.254, 0.215, 0.178, 0.089, 0.052, 0.035, 0.015]
    prob_θ::Vector{Vector{Float64}} = [prob_θ_wm, prob_θ_ww, prob_θ_mm, prob_θ_mw] #all together

    #parameters for SAT Math/Verbal (note: does not vary by demographics)
    param_m::Vector{Float64} = [0.0, 1.0, 0.611, 0.616, 0.165, -0.040, -0.004, 0.006]
    param_v::Vector{Float64} = [-0.661, 1.029, 0.611, 0.616, 0.165, -0.040, -0.004, 0.006]

    #parameters for CQ (do depend on demographics)
    param_cq_wm::Vector{Float64} = [1.313, 0.626, 0.693, 0.351, -0.128, 0.067, -0.063, 0.013] #white men
    param_cq_ww::Vector{Float64} = [1.234, 0.544, 0.694, 0.341, -0.128, 0.083, -0.062, 0.011] #white women
    param_cq_mm::Vector{Float64} = [0.365, 0.706, 0.583, 0.488, -0.116, 0.022, 0.022, 0.020] #minority men
    param_cq_mw::Vector{Float64} = [0.403, 0.643, 0.628, 0.492, -0.092, 0.037, 0.004, 0.016] #minority women
    param_cq::Vector{Vector{Float64}} = [param_cq_wm, param_cq_ww, param_cq_mm, param_cq_mw] #all together
end


#main function for SAT math. Inputs normalized test score/quality, score we're using, and agent type. outputs theta probabilities
function post_dist(prim::Primitives, score::Float64, type::Int64, dem::Int64)
    @unpack θ, prob_θ, param_m, param_v, param_cq, ϵ_0 = prim
    nθ = length(θ)
    prob_θ_post = zeros(nθ)

    #figure out which kind of paramters we're using. If type=1, doing SAT math. 2, SAT verbal. 3, College quality
    p_θs = prob_θ[dem] #probabilities of θ depending on demographic group
    params = copy(param_m) #don't depend on demographics
    if type==2
        params = copy(param_v) #don't depend on demographics
    elseif type == 3
        params = copy(param_cq[dem]) #do depend on demographics.
    end

    #unpack a bit to make things more familiar
    b0 = params[1] #constant
    b1 = params[2] #ability slope
    σ = params[3] #sigma
    x_poly = params[4:8] #density polynomial

    #solve for probabilities of θ given demographics and scorwe
    p_score = pdf(Normal(), score) #pdf of score, assuming standard normal.
    for i = 1:nθ #loop over possible θs
        p_θ = p_θs[i] #probability of specific θ given demographic group
        x = (score - (b0 + b1 * θ[i]))/σ #compute measurement error, normalized by σ
        x_vec = [x^0, x^1, x^2, x^3, x^4] #polynomial

        #density of measurement error
        p_score_θ = (((sum(x_vec.*x_poly))^2) * exp((-x^2)/2) + ϵ_0 * pdf(Normal(), x))

        #application of bayes rule
        prob_θ_post[i] = (p_score_θ * p_θ) / p_score
    end

    #return deliverable
    return prob_θ_post
end


#testing behavior with college quality
prim = Primitives()
temp = post_dist(prim, 0.0, 3, 1)

#higher CQ score results in higher probabilities of large θs
temp = post_dist(prim, 1.0, 3, 1)

#but lowr probabilities of high θs for minority women (for example) due to different measuremnet error distrbution
temp = post_dist(prim, 1.0, 3, 4)

#converse also works.
temp = post_dist(prim, -1.0, 3, 1)
temp = post_dist(prim, -1.0, 3, 4)


####
