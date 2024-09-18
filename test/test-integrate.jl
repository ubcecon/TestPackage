@testset "TestPackage.jl" begin
  using Distributions, LinearAlgebra
  dx = Normal(0,1)
  N = 10000


  ∫mc = MCIntegrator(dx, N)
  ∫s = SobolIntegrator(dx, N)
  @test ∫mc(x->x^2) ≈ 1  rtol=6/sqrt(N)
  @test isapprox(∫s(x->x^2), 1, rtol=6*log(N)/N)
  ∫gh = GaussHermite(dx, 10)
  @test isapprox(∫gh(x->x^2),1)
  ∫c = CubatureIntegrator(dx, rtol=1e-6)
  @test isapprox(∫c(x->x[1]^2)[1],1, rtol=1e-6)

  Σ = ones(2,2) + I
  dX = MvNormal(zeros(2), Σ)
  ∫mc = MCIntegrator(dX, N)
  ∫s = SobolIntegrator(dX, N)

  @test isapprox(norm(∫mc(x->x*x') - Σ), 0, atol=12/sqrt(N))
  @test isapprox(norm(∫s(x->x*x') - Σ), 0, atol=3*log(N)^2/N)
  ∫gh = GaussHermite(dX, 100)
  @test isapprox(norm(∫gh(x->x*x') - Σ), 0, atol=1e-4)
  ∫sg = SparseGrid(dX, 5)
  @test isapprox(norm(∫sg(x->x*x') - Σ), 0, atol=1e-4)

  ∫c = CubatureIntegrator(dX, rtol=1e-6)
  @test isapprox(norm(∫c(x->x*x') - Σ), 0, atol=1e-4)


end
