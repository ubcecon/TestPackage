@testset "TestPackage.jl" begin
  using Distributions, LinearAlgebra
  dx = Normal(0,1)
  Σ = ones(2,2) + I
  dX = MvNormal(zeros(2), Σ)
  N = 10000
  @test isapprox(∫mc(x->x^2, dx, ndraw=N), 1, rtol=6/sqrt(N))
  @test isapprox(∫s(x->x^2, dx, ndraw=N), 1, rtol=6*log(N)/N)
  @test isapprox(norm(∫mc(x->x*x', dX, ndraw=N) - Σ), 0, atol=12/sqrt(N))
  @test isapprox(norm(∫s(x->x*x', dX, ndraw=N) - Σ), 0, atol=6*log(N)^2/N)

end
