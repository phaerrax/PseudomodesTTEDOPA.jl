export exchange_interaction, spin_chain

"""
    spin_chain(freqs, coups, sites::Vector{<:Index})

Return an OpSum object encoding the Hamiltonian part ``-i[H, –]`` of the GKSL equation
for a spin chain of frequencies `freqs` and coupling constants `coups`, on `sites`.
"""
function spin_chain(::SiteType"vS=1/2", freqs, coups, sites::Vector{<:Index})
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumber.(sites))
        ℓ += freqs[j] * gkslcommutator("N", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumber.(sites), 2, 1))
        jws = jwstring(; start=site1, stop=site2)
        ℓ +=
            coups[j + 1] * (
                gkslcommutator("σ+", site1, jws..., "σ-", site2) +
                gkslcommutator("σ-", site1, jws..., "σ+", site2)
            )
    end
    return ℓ
end

function spin_chain(::SiteType"vElectron", freqs, coups, sites::Vector{<:Index})
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumber.(sites))
        ℓ += freqs[j] * gkslcommutator("Ntot", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumber.(sites), 2, 1))
        jws = jwstring(; start=site1, stop=site2)
        ℓ += +coups[j + 1] * gkslcommutator("Aup†F", site1, jws..., "Aup", site2)
        ℓ += -coups[j + 1] * gkslcommutator("AupF", site1, jws..., "Aup†", site2)
        ℓ += +coups[j + 1] * gkslcommutator("Adn†", site1, jws..., "FAdn", site2)
        ℓ += -coups[j + 1] * gkslcommutator("Adn", site1, jws..., "FAdn†", site2)
    end
    return ℓ
end

"""
    exchange_interaction(s1::Index, s2::Index)

Return an OpSum object encoding the Hamiltonian part ``-i[H, –]`` of an exchange interaction
between sites `s1` and `s2` term in a GKSL equation.
"""
function exchange_interaction(::SiteType"vS=1/2", s1::Index, s2::Index)
    site1 = sitenumber(s1)
    site2 = sitenumber(s2)

    ℓ = OpSum()
    jws = jwstring(; start=site1, stop=site2)
    ℓ += (
        gkslcommutator("σ+", site1, jws..., "σ-", site2) +
        gkslcommutator("σ-", site1, jws..., "σ+", site2)
    )
    return ℓ
end

function exchange_interaction(::SiteType"vElectron", s1::Index, s2::Index)
    site1 = sitenumber(s1)
    site2 = sitenumber(s2)
    # c↑ᵢ† c↑ᵢ₊₁ + c↑ᵢ₊₁† c↑ᵢ + c↓ᵢ† c↓ᵢ₊₁ + c↓ᵢ₊₁† c↓ᵢ =
    # a↑ᵢ† Fᵢ a↑ᵢ₊₁ - a↑ᵢ Fᵢ a↑ᵢ₊₁† + a↓ᵢ† Fᵢ₊₁ a↓ᵢ₊₁ - a↓ᵢ Fᵢ₊₁ a↓ᵢ₊₁†
    ℓ = OpSum()
    jws = jwstring(; start=site1, stop=site2)
    ℓ += (
        gkslcommutator("Aup†F", site1, jws..., "Aup", site2) -
        gkslcommutator("AupF", site1, jws..., "Aup†", site2) +
        gkslcommutator("Adn†", site1, jws..., "FAdn", site2) -
        gkslcommutator("Adn", site1, jws..., "FAdn†", site2)
    )
    return ℓ
end

