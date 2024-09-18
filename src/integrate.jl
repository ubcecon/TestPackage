module Integrate

import Statistics, HCubature, Sobol, FastGaussQuadrature, SparseGrids
import Base.Iterators: take, product, repeated
import LinearAlgebra: I, cholesky
import Distributions: Distribution, Normal, MvNormal, quantile, pdf, AbstractMvNormal

export AbstractIntegrator, FixedIntegrator, MCIntegrator, SobolIntegrator, GaussHermite, SparseGrid, CubatureIntegrator

abstract type AbstractIntegrator end

struct FixedIntegrator{TX, TW} <: AbstractIntegrator
  x::TX
  w::TW
end

function (∫::FixedIntegrator)(f::Function)
  sum(f(∫.x[i])*∫.w[i] for i in eachindex(∫.x))
end

MCIntegrator(dx::Distribution,ndraw) = FixedIntegrator([rand(dx) for i in 1:ndraw], fill(1/ndraw, ndraw))

function SobolIntegrator(dx::Distribution,ndraw)
  if (length(dx)) != 1
    error("SobolIntegrator only works for 1D distributions or MvNormal")
  end
  invcdf(x) = quantile(dx,x)
  ss = skip(Sobol.SobolSeq(length(dx)),ndraw)
  x = [invcdf(x[1]) for x in take(ss,ndraw)]
  w = fill(1/ndraw, ndraw)
  FixedIntegrator(x,w)
end

function SobolIntegrator(dx::AbstractMvNormal,ndraw)
  marginals = [Normal(0.0,1.0) for i in 1:length(dx)]
  invcdf(x) = quantile.(marginals,x)
  ss = skip(Sobol.SobolSeq(length(dx)),ndraw)
  L = cholesky(dx.Σ).L
  x = [L*invcdf(x) + dx.μ for x in take(ss,ndraw)]
  w = fill(1/ndraw, ndraw)
  FixedIntegrator(x,w)
end

function GaussHermite(dx::MvNormal, ndraw)
  n = Int(ceil(ndraw^(1/length(dx))))
  x, w = FastGaussQuadrature.gausshermite(n)
  L = cholesky(dx.Σ).L
  x = [(√2*L*vcat(xs...) + dx.μ) for xs in product(repeated(x, length(dx))...)]
  w = [prod(ws)/π^(length(dx)/2) for ws in product(repeated(w, length(dx))...)]
  FixedIntegrator(x,w)
end

function GaussHermite(dx::Normal, n)
  x, w = FastGaussQuadrature.gausshermite(n)
  x = (sqrt(2)*dx.σ*x .+ dx.μ)
  w = w ./ sqrt(π)
  FixedIntegrator(x,w)
end


function SparseGrid(dx::MvNormal, order)
  X, W = SparseGrids.sparsegrid(length(dx), order, FastGaussQuadrature.gausshermite, sym=true)
  L = cholesky(dx.Σ).L
  x = ([(√2*L*x + dx.μ) for x in X])
  w = W ./ π^(length(dx)/2)
  FixedIntegrator(x,w)
end

struct CubatureIntegrator{VT, F} <: AbstractIntegrator
  rtol::Float64
  atol::Float64
  lb::VT
  ub::VT
  transform::F
end

function CubatureIntegrator(dx::Union{Normal, AbstractMvNormal}; rtol=1e-6, atol=1e-6)
  x(t) = t./(1 .- t.^2)
  Dx(t) = prod((1 .+ t.^2)./(1 .- t.^2).^2)
  CubatureIntegrator(rtol, atol, -ones(length(dx)), ones(length(dx)),
    (f,t)->f(x(t))*pdf(dx,x(t))*Dx(t))
end

function (∫::CubatureIntegrator)(f::Function)
  HCubature.hcubature(t->∫.transform(f,t), ∫.lb, ∫.ub, rtol=∫.rtol, atol=∫.atol)[1]
end

end
