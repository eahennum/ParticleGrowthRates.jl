module ParticleGrowthRates

function growth_rate_sffk(r, x_mg, x_si, dc_mg, dc_si, df, T, γ)
    d_eff = (((2/3 - x_mg)^2 / (x_mg*dc_mg)) + ((1/3 - x_si)^2 / (x_si*dc_si)))
    A = R * T * r * d_eff/1e-5
    B = df - 2γ/r

    return B/A
end

end # module ParticleGrowthRates
