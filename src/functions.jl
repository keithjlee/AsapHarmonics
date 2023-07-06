function spherical_gaussian(func::Matrix{Float64}, θ::Float64, ϕ::Float64; δ = 1)
    val = 0.0

    px = cos(ϕ)sin(θ)
    py = sin(ϕ)sin(θ)
    pz = cos(θ)

    for i = 1:size(func)[1] #for each force
        val += func[i,4] * exp(-δ*((px-func[i,1])^2 + (py-func[i,2])^2 + (pz-func[i,3])^2))
    end
    return val
end

function generate_feature_vector(forcefunction::Matrix{Float64}; dims = 16)

    coefficients = sph_transform(forcefunction)

    C = Vector{Vector{Float64}}()

    for l = 0:dims-1
        push!(C, [coefficients[sph_mode(l,m)] for m = -l:l])
    end

    [norm(sum(C[i] .* Y[i])) for i = 1:dims]
end