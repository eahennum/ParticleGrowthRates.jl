using Test

using ThermAl: BETA_MG2SI

using ParticleGrowthRates: growth_rate_sffk

temp = 773.15
dc_mg = 5.348638269137225e-14
dc_si = 1.5663062228869834e-13
x_mg = 0.0070
x_si = 0.0035

# growth_rate_sffk(BETA_MG2SI, 5e-9, 0.0070, 0.0035, dc_mg, dc_si, temp, 0.36)

@testset "SFFK growth model" begin
    drdt = growth_rate_sffk(BETA_MG2SI, 5e-9, 0.0070, 0.0035, dc_mg, dc_si, temp, 0.36)
    @test drdt == -2.9942624463712036e-8
end