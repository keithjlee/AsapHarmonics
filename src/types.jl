"""
    NodeForces(node, model, C, forces; δ = 20, normalize_forces = false,
               store_function = false, resolution = 90)

Per-node force signature: outward unit directions (columns of
`forcepositions`) and signed magnitudes of all forces acting at `node` —
incident member axial forces, applied `NodeLoad`s, and the support reaction.

The sampled spherical signature `forcefunction` (for visualization and
numeric cross-validation) is only computed when `store_function = true`, on
the `sph_points(resolution + 1)` grid; it is otherwise left empty. Feature
vectors do not need it — they are computed in closed form from
`forcepositions`/`forcemagnitudes`.
"""
mutable struct NodeForces
    node::Node{Float64}
    forcepositions::Matrix{Float64}
    forcemagnitudes::Vector{Float64}
    forcefunction::Matrix{Float64}

    function NodeForces(node::Node, model::Model, C::SparseMatrixCSC{Int64, Int64}, forces::Vector{Float64}; δ = 20, normalize_forces = false, store_function = false, resolution = 90)
        i = node.index

        #indices of elements connected to node i
        i_connected = findall(.!iszero.(C[:, i]))

        #-1 if node i is start point, 1 if node i is end point
        start_end = vec(C[:, i][i_connected])

        #external forces
        # external_loads = getproperty.(model.loads[node.loadIDs], :value)
        external_loads = [load.value for load in model.loads if load isa NodeLoad && load.node.index == i]


        #reaction if applicable
        rxn3 = collect(reaction(model.results, node)[1:3])
        if !iszero(norm(rxn3))
            push!(external_loads, rxn3)
        end

        if normalize_forces
            external_loads = normalize.(external_loads)
        end

        #force values
        node_forces = forces[i_connected]
        
        if normalize_forces
            node_forces ./= abs.(node_forces)
        end

        if !isempty(external_loads)
            node_forces = [node_forces; norm.(external_loads)]
        end

        #force positions
        force_positions = [-normalize(e.nodeEnd.position - e.nodeStart.position) .* factor for (e, factor) in zip(model.elements[i_connected], start_end)]
        if !isempty(external_loads)
            force_positions = [force_positions; normalize.(external_loads)]
        end


        #sampled spherical signature (visualization / numeric validation only)
        func = store_function ?
            sampled_force_function(force_positions, node_forces; delta = δ, nlat = resolution + 1) :
            Matrix{Float64}(undef, 0, 0)

        new(node, hcat(force_positions...), node_forces, func)
        
    end
end

"""
    HarmonicAnalysis(model; delta = 20, dims = 16, normalize_forces = false,
                     store_functions = false, resolution = 90)

Spherical-harmonic shape-descriptor analysis of every node of a solved 3D
`Model`: builds a [`NodeForces`](@ref) signature per node and computes its
rotation-invariant feature vector (degrees `l = 0:dims-1`) in closed form via
[`spherical_feature_vector`](@ref).

`delta` is the sharpness of the Gaussian force bumps (each bump approaches a
Dirac spike as `delta → ∞`). With `store_functions = true` the sampled
spherical signatures are also retained in `forcefunctions` for visualization.
"""
mutable struct HarmonicAnalysis
    model::Model{Float64}
    nodeforces::Vector{NodeForces}
    forcefunctions::Vector{Matrix{Float64}}
    featurevectors::Vector{Vector{Float64}}

    function HarmonicAnalysis(model::Model; delta = 20, dims = 16, normalize_forces = false, store_functions = false, resolution = 90)
        C = connectivity(model)
        forces = [axial_force(model.results, el) for el in model.elements]

        nodeforces = [NodeForces(node, model, C, forces; δ = delta, normalize_forces = normalize_forces, store_function = store_functions, resolution = resolution) for node in model.nodes]

        forcefunctions = getproperty.(nodeforces, :forcefunction)

        featurevectors = [spherical_feature_vector(nf.forcepositions, nf.forcemagnitudes; delta = delta, dims = dims) for nf in nodeforces]

        new(model, nodeforces, forcefunctions, featurevectors)
    end
end

"""
    NodeForces2d(node, model, C, forces; δ = 0.1, n = 90,
                 normalize_forces = false, store_function = false)

Per-node planar force signature: unit directions and signed magnitudes of all
forces acting at `node` in the XY plane (member axial forces via
`local_frame`, applied loads, reactions — compression/tension encoded by
sign and direction flip). The `n`-point sampled circular signature
`forcefunction` is only computed when `store_function = true`.
"""
mutable struct NodeForces2d
    node::Node{Float64}
    forcepositions::Vector{Vector{Float64}}
    forcemagnitudes::Vector{Float64}
    forcefunction::Vector{Float64}

    function NodeForces2d(node::Node, model::Model, C::SparseMatrixCSC{Int64, Int64}, forces::Vector{Float64}; δ = 0.1, n = 90, normalize_forces = false, store_function = false)
        i = node.index

        #indices of elements connected to node i
        i_connected = findall(.!iszero.(C[:, i]))

        #-1 if node i is start point, 1 if node i is end point
        start_end = collect(C[:, i][i_connected])

        #external forces
        external_loads = Vector{Float64}()
        external_load_positions = Vector{Vector{Float64}}()

        for load in model.loads
            if load isa NodeLoad && load.node.index == i
                val = collect(load.value[1:2])

                #position of load vector (default top)
                position_vector = normalize(val)
                if val[2] < 0 #this is a compressive force acting downwards
                    push!(external_load_positions, position_vector .* -1)
                    push!(external_loads, -norm(val))
                else #else it is a tensile force acting upwards
                    push!(external_load_positions, normalize(val))
                    push!(external_loads, norm(val))
                end
            end
        end

        #reaction if applicable
        rxn6 = reaction(model.results, model.nodes[i])
        if !iszero(norm(rxn6))

            rxn = collect(rxn6[1:2])

            #position of reaction vector (default bottom
            position_vector = normalize(rxn)
            if rxn[2] > 0 #this is a compressive force acting upwards
                push!(external_loads, -norm(rxn))
                push!(external_load_positions, position_vector .* -1)
            else #else it is a tensile force acting downwards
                push!(external_loads, norm(rxn))
                push!(external_load_positions, position_vector)
            end
        end

        #force values
        node_forces = forces[i_connected]
        if !isempty(external_loads)
            node_forces = [node_forces; external_loads]
        end

        if normalize_forces
            node_forces .= 1.0
        end

        #force positions
        force_positions = [-local_frame(e)[1, 1:2] .* factor for (e, factor) in zip(model.elements[i_connected], start_end)]
        if !isempty(external_loads)
            force_positions = [force_positions; external_load_positions]
        end

        #sampled circular signature (visualization / numeric validation only)
        force_function = store_function ? circular_gaussian(force_positions, node_forces, δ, n) : Float64[]

        #return
        new(node, force_positions, node_forces, force_function)
    end
end

"""
    HarmonicAnalysis2d(model; delta = 0.1, n = 90, dims = 16,
                       normalize_forces = false, store_functions = false)

Fourier shape-descriptor analysis of every node of a solved **planar (XY)**
`Model`: builds a [`NodeForces2d`](@ref) signature per node and computes its
rotation-invariant feature vector (frequencies `k = 0:dims-1`) in closed form
via [`circular_feature_vector`](@ref).

`delta` is the angular width σ of the Gaussian force bumps. With
`store_functions = true` the `n`-point sampled circular signatures are also
retained in `forcefunctions` for visualization (note: a raw `rfft` of those
samples is ≈ `n×` the closed-form feature values).
"""
mutable struct HarmonicAnalysis2d
    model::Model{Float64}
    nodeforces::Vector{NodeForces2d}
    forcefunctions::Vector{Vector{Float64}}
    featurevectors::Vector{Vector{Float64}}

    function HarmonicAnalysis2d(model::Model; delta = 0.1, n = 90, dims = 16, normalize_forces = false, store_functions = false)
        C = connectivity(model)
        forces = [axial_force(model.results, el) for el in model.elements]

        nodeforces = [NodeForces2d(node, model, C, forces; δ = delta, n = n, normalize_forces = normalize_forces, store_function = store_functions) for node in model.nodes]

        forcefunctions = getproperty.(nodeforces, :forcefunction)

        featurevectors = [circular_feature_vector(nf.forcepositions, nf.forcemagnitudes; sigma = delta, dims = dims) for nf in nodeforces]

        new(model, nodeforces, forcefunctions, featurevectors)
    end

end