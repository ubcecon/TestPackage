@testset "delta" begin
  using Distributions, LinearAlgebra

  J = 3
  K = 2
  x = randn(J,K)
  C = randn(K,K)
  Σ = C'*C + I
  dFν = MvNormal(zeros(K), I)
  ∫ = SparseGrid(dFν, 5)
  for i in 1:10
    δ = randn(J)*0.3
    δ[1] = 0
    s = share(δ, Σ, x, ∫)
    @test isapprox(δ, delta(s, Σ, x, ∫), atol=1e-6)
  end

end
