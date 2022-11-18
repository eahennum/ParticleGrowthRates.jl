module ParticleGrowthRates

using ThermAl: AbstractPhase, delta_gibbs

const __UNIV_GAS_CONST = 8.314456


# Calculatates the multicomponent effective diffusion effective diffusion coefficient
# TODO: move this function to another package / base library
function __effective_diffusion_coefficient(x_mg, x_si, dc_mg, dc_si)
    return (((2/3 - x_mg)^2 / (x_mg*dc_mg)) + ((1/3 - x_si)^2 / (x_si*dc_si))) / 1e-5
end

# Calculates the chemical driving force
function __driving_force(p::AbstractPhase, x_mg, x_si, temp)

    return  delta_gibbs(p, x_mg, x_si, temp) / 1e-5
end


"""
    growth_rate_sffk(r, x_mg, x_si, dc_mg, dc_si, γ, temp)

Calculates the growth rate of a spherical particle.
"""
function growth_rate_sffk(p::AbstractPhase, r, x_mg, x_si, dc_mg, dc_si, temp, γ)
    d = __effective_diffusion_coefficient(x_mg, x_si, dc_mg, dc_si)
    f = __driving_force(p, x_mg, x_si, temp)

    return growth_rate_sffk(r, f, d, γ, temp)
end

# Growth rate for a stoichiometric particle
function growth_rate_sffk(r, f, d, γ, t)

    return (f - 2γ / r) / (__UNIV_GAS_CONST * t * r * d) 
end

end # module ParticleGrowthRates
