module TestPackage

using Distributions, Statistics
using HCubature
import Sobol
import Base.Iterators: take
using FastGaussQuadrature, LinearAlgebra
import Base.Iterators: product, repeated
using SparseGrids

export ∫mc, ∫s, ∫q, ∫sgq, ∫cuba


∫mc(f, dx; ndraw=100)=mean(f(rand(dx)) for i in 1:ndraw)

function ∫s(f,dx::AbstractMvNormal;ndraw=100)
  marginals = [Normal(dx.μ[i], sqrt(dx.Σ[i,i])) for i in 1:length(dx)]
  invcdf(x) = quantile.(marginals,x)
  ss = skip(Sobol.SobolSeq(length(dx)),ndraw)
  mean(f(invcdf(x)) for x in take(ss,ndraw))
end

function ∫s(f,dx::Normal;ndraw=100)
  invcdf(x) = quantile(dx,x)
  ss = skip(Sobol.SobolSeq(length(dx)),ndraw)
  mean(f(invcdf(x[1])) for x in take(ss,ndraw))
end

function ∫q(f, dx::MvNormal; ndraw=100)
  n = Int(ceil(ndraw^(1/length(dx))))
  x, w = gausshermite(n)
  L = cholesky(dx.Σ).L
  sum(f(√2*L*vcat(xs...) + dx.μ)*prod(ws)
      for (xs,ws) ∈ zip(product(repeated(x, length(dx))...),
                        product(repeated(w, length(dx))...))
        )/(π^(length(dx)/2))
end

function ∫sgq(f, dx::MvNormal; order=5)
  X, W = sparsegrid(length(dx), order, gausshermite, sym=true)
  L = cholesky(dx.Σ).L
  sum(f(√2*L*x + dx.μ)*w for (x,w) ∈ zip(X,W))/(π^(length(dx)/2))
end


function ∫cuba(f, dx; rtol=1e-4)
  D = length(dx)
  x(t) = t./(1 .- t.^2)
  Dx(t) = prod((1 .+ t.^2)./(1 .- t.^2).^2)
  hcubature(t->f(x(t))*pdf(dx,x(t))*Dx(t), -ones(D),ones(D), rtol=rtol)[1]
end


end
