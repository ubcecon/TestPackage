import SimpleNonlinearSolve
using LinearAlgebra

@doc raw"""
    share(δ, Σ, x, ∫)

Computes shares in random coefficient logit with mean tastes `δ`, observed characteristics `x`, unobserved taste distribution that can be integrated by ∫, and taste covariances `Σ`.

# Arguments

- `δ` vector of length `J`
- `Σ` `K` by `K` matrix
- `x` `J` by `K` array
- `∫` function that integrates over the distribution of ν

# Returns

- vector of length `J` consisting of $s_1$, ..., $s_J$
"""
function share(δ, Σ, x, ∫)
  J,K = size(x)
  (length(δ) == J) || error("length(δ)=$(length(δ)) != size(x,1)=$J")
  (K,K) === size(Σ) || error("size(x,2)=$K != size(Σ)=$(size(Σ))")
  function shareν(ν)
    s = δ + x*Σ*ν
    s .-= maximum(s)
    s .= exp.(s)
    s ./= sum(s)
    return(s)
  end
  return(∫(shareν))
end


function delta(s, Σ, x, ∫)
  p = (Σ=Σ, x=x)
  F(d, p) = share([0, d...],p.Σ, p.x,∫) - s
  d0 = zeros(length(s)-1)
  prob = SimpleNonlinearSolve.NonlinearProblem(F, d0, p)
  sol = SimpleNonlinearSolve.solve(prob,
    SimpleNonlinearSolve.SimpleNewtonRaphson(),show_trace=Val(false)) #, trace_level=TraceAll())
  return [0, sol.u...]
end
