mutable struct NodeForces
    node::Node{Float64}
    forcepositions::Matrix{Float64}
    forcemagnitudes::Vector{Float64}
    forcefunction::Matrix{Float64}

    function NodeForces(node::Node, model::Model, C::SparseMatrixCSC{Int64, Int64}, forces::Vector{Float64}; δ = 20, normalize_forces = false)
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


        #make spherical gaussian function
        sph = [getindex.(force_positions, 1) getindex.(force_positions, 2) getindex.(force_positions, 3) node_forces]
        func = [spherical_gaussian(sph, t, p; δ = δ) for t in thetarange, p in phirange]

        new(node, hcat(force_positions...), node_forces, func)
        
    end
end

mutable struct HarmonicAnalysis
    model::Model{Float64}
    nodeforces::Vector{NodeForces}
    forcefunctions::Vector{Matrix{Float64}}
    featurevectors::Vector{Vector{Float64}}

    function HarmonicAnalysis(model::Model; delta = 20, dims = 16, normalize_forces = false)
        C = connectivity(model)
        forces = [axial_force(model.results, el) for el in model.elements]

        nodeforces = [NodeForces(node, model, C, forces; δ = delta, normalize_forces = normalize_forces) for node in model.nodes]

        forcefunctions = getproperty.(nodeforces, :forcefunction)

        featurevectors = spherical_feature_vector.(forcefunctions; dims = dims)

        new(model, nodeforces, forcefunctions, featurevectors)
    end
end

mutable struct NodeForces2d
    node::Node{Float64}
    forcepositions::Vector{Vector{Float64}}
    forcemagnitudes::Vector{Float64}
    forcefunction::Vector{Float64}

    function NodeForces2d(node::Node, model::Model, C::SparseMatrixCSC{Int64, Int64}, forces::Vector{Float64}; δ = 0.1, n = 90, normalize_forces = false)
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

        #force functions
        force_function = circular_gaussian(force_positions, node_forces, δ, n)
            
        #return
        new(node, force_positions, node_forces, force_function)
    end
end

mutable struct HarmonicAnalysis2d
    model::Model{Float64}
    nodeforces::Vector{NodeForces2d}
    forcefunctions::Vector{Vector{Float64}}
    featurevectors::Vector{Vector{Float64}}

    function HarmonicAnalysis2d(model::Model; delta = 0.1, n = 90, dims = 16, normalize_forces = false)
        C = connectivity(model)
        forces = [axial_force(model.results, el) for el in model.elements]

        nodeforces = [NodeForces2d(node, model, C, forces; δ = delta, n = n, normalize_forces = normalize_forces) for node in model.nodes]

        forcefunctions = getproperty.(nodeforces, :forcefunction)

        featurevectors = [circular_feature_vector(ff; dims = dims) for ff in forcefunctions]

        new(model, nodeforces, forcefunctions, featurevectors)
    end

end