mutable struct NodeForce
    node::TrussNode
    forcepositions::Matrix{Float64}
    forcemagnitudes::Vector{Float64}
    forcefunction::Matrix{Float64}

    function NodeForce(node::TrussNode, model::TrussModel, C::SparseMatrixCSC{Int64, Int64}, forces::Vector{Float64}; δ = 20)
        i = node.nodeID

        #indices of elements connected to node i
        i_connected = findall(.!iszero.(C[:, i]))

        #-1 if node i is start point, 1 if node i is end point
        start_end = vec(C[:, i][i_connected])

        #external forces
        external_loads = getproperty.(model.loads[node.loadIDs], :value)

        #reaction if applicable
        if !iszero(norm(node.reaction))
            push!(external_loads, node.reaction)
        end

        #force values
        node_forces = forces[i_connected]
        if !isempty(external_loads)
            node_forces = [node_forces; norm.(external_loads)]
        end

        #force positions
        force_positions = [-normalize(e.nodeEnd.position - e.nodeStart.position) .* factor for (e, factor) in zip(model.elements[i_connected], start_end)]
        if !isempty(external_loads)
            force_positions = [force_positions; normalize.(external_loads)]
        end


        #make spherical gaussian function
        sph = [getindex.(force_positions, 1) getindex.(force_positions, 2) getindex.(force_positions, 3) node_forces]
        func = [spherical_gaussian(sph, t, p; δ = δ) for t in thetarange, p in phirange]

        new(node, hcat(force_positions...), node_forces, func)
        
    end

    function NodeForce(i::Int64, model::TrussModel, C::SparseMatrixCSC{Int64, Int64}, forces::Vector{Float64}; δ = 20)
        
        node = model.nodes[i]

        #indices of elements connected to node i
        i_connected = findall(.!iszero.(C[:, i]))

        #-1 if node i is start point, 1 if node i is end point
        start_end = vec(C[:, i][i_connected])

        #external forces
        external_loads = getproperty.(model.loads[node.loadIDs], :value)

        #reaction if applicable
        if !iszero(norm(node.reaction))
            push!(external_loads, node.reaction)
        end

        #force values
        node_forces = forces[i_connected]
        if !isempty(external_loads)
            node_forces = [node_forces; norm.(external_loads)]
        end

        #force positions
        force_positions = [-normalize(e.nodeEnd.position - e.nodeStart.position) .* factor for (e, factor) in zip(model.elements[i_connected], start_end)]
        if !isempty(external_loads)
            force_positions = [force_positions; normalize.(external_loads)]
        end

        #make spherical gaussian function
        sph = [getindex.(force_positions, 1) getindex.(force_positions, 2) getindex.(force_positions, 3) node_forces]
        func = [spherical_gaussian(sph, t, p; δ = δ) for t in thetarange, p in phirange]

        new(node, hcat(force_positions...), node_forces, func)
        
    end
end

mutable struct HarmonicAnalysis
    model::TrussModel
    nodeforces::Vector{NodeForce}
    forcefunctions::Vector{Matrix{Float64}}
    featurevectors::Vector{Vector{Float64}}

    function HarmonicAnalysis(model::TrussModel; delta = 20, dims = 16)
        C = connectivity(model)
        forces = axialforce(model.elements)

        nodeforces = [NodeForce(node, model, C, forces; δ = delta) for node in model.nodes]

        forcefunctions = getproperty.(nodeforces, :forcefunction)

        featurevectors = generate_feature_vector.(forcefunctions; dims = dims)

        new(model, nodeforces, forcefunctions, featurevectors)
    end
end