module TestPackage


include("integrate.jl")
using .Integrate

include("shares.jl")


export AbstractIntegrator, FixedIntegrator, MCIntegrator, SobolIntegrator, GaussHermite, SparseGrid, CubatureIntegrator, share, delta


end
