# AsapHarmonics v1.0 smoke suite (first test suite for this package):
# harmonic force signatures build from a solved v1.0 model and produce
# finite feature vectors of the expected shape.

using AsapHarmonics
using Asap
using Test
using LinearAlgebra

function small_truss()
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]
    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([false, false, false], rot))
    n3 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot))
    n4 = Node([2.0, 3.0, 4.0], vcat([true, true, true], rot))
    els = AbstractElement{Float64}[
        TrussElement(n1, n3, sec), TrussElement(n2, n3, sec),
        TrussElement(n1, n4, sec), TrussElement(n2, n4, sec),
        TrussElement(n3, n4, sec)]
    loads = AbstractLoad{Float64}[
        NodeForce(n3, [0.0, -50.0, 0.0]), NodeForce(n4, [10.0, -20.0, 5.0])]
    model = Model([n1, n2, n3, n4], els, loads)
    solve!(model)
    return model
end

@testset "AsapHarmonics v1.0" begin
    model = small_truss()

    ha = HarmonicAnalysis(model; dims = 8)
    @test length(ha.nodeforces) == length(model.nodes)
    @test all(fv -> length(fv) == 8 && all(isfinite, fv), ha.featurevectors)
    @test all(nf -> all(isfinite, nf.forcefunction), ha.nodeforces)

    # loaded node's signature includes its external load; supported node its reaction
    loaded = ha.nodeforces[3]
    @test length(loaded.forcemagnitudes) >= 3          # connected members + load

    # the 2D analysis expects a PLANAR (XY) truss
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]
    p1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    p2 = Node([4.0, 0.0, 0.0], vcat([false, false, false], rot))
    p3 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot))
    pels = AbstractElement{Float64}[TrussElement(p1, p3, sec), TrussElement(p2, p3, sec)]
    pmodel = Model([p1, p2, p3], pels,
        AbstractLoad{Float64}[NodeForce(p3, [0.0, -50.0, 0.0])])
    solve!(pmodel)

    ha2 = HarmonicAnalysis2d(pmodel; dims = 8, n = 45)
    @test length(ha2.nodeforces) == length(pmodel.nodes)
    @test all(fv -> all(isfinite, fv), ha2.featurevectors)
end
