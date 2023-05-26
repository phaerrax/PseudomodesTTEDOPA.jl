export spin_chain, spin_chain_adjoint, spin_chain′
export exchange_interaction, exchange_interaction_adjoint, exchange_interaction′

spin_chain(::SiteType, ::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Int}) = nothing

"""
    spin_chain(
        freqs::Vector{<:Real},
        coups::Vector{<:Real},
        sites::Vector{<:Index}
    )

Return an OpSum object encoding the Hamiltonian part ``-i[H, –]`` of the GKSL equation
for a spin chain of frequencies `freqs` and coupling constants `coups`, on `sites`.
"""
function spin_chain(freqs::Vector{<:Real}, coups::Vector{<:Real}, sites::Vector{<:Index})
    stypes = sitetypes(first(sites))
    for st in stypes
        # Check if all sites have this type (otherwise skip this tag).
        if all(i -> st in sitetypes(i), sites)
            # If the type is shared, then try calling the function with it.
            ℓ = spin_chain(st, freqs, coups, sitenumber.(sites))
            # If the result is something, return that result.
            if !isnothing(ℓ)
                return ℓ
            end
        end
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"spin_chain\" function not found for Index tags $(tags(sites[1]))"
        ),
    )
end

function spin_chain(
    ::SiteType"Fermion",
    freqs::Vector{<:Real},
    coups::Vector{<:Real},
    sitenumbers::Vector{<:Int},
)
    h = OpSum()

    for (j, site) in enumerate(sitenumbers)
        h += freqs[j], "n", site
    end

    for (j, (site1, site2)) in enumerate(partition(sitenumbers, 2, 1))
        h += coups[j + 1], "c†", site1, "c", site2
        h += coups[j + 1], "c†", site2, "c", site1
    end

    return h
end

function spin_chain(
    ::SiteType"Electron",
    freqs::Vector{<:Real},
    coups::Vector{<:Real},
    sitenumbers::Vector{<:Int},
)
    h = OpSum()

    for (j, site) in enumerate(sitenumbers)
        h += freqs[j], "n↓", site
        h += freqs[j], "n↑", site
    end

    for (j, (site1, site2)) in enumerate(partition(sitenumbers, 2, 1))
        h += coups[j + 1], "c†↓", site1, "c↓", site2
        h += coups[j + 1], "c†↓", site2, "c↓", site1

        h += coups[j + 1], "c†↑", site1, "c↑", site2
        h += coups[j + 1], "c†↑", site2, "c↑", site1
    end

    return h
end

function spin_chain(
    ::SiteType"vS=1/2",
    freqs::Vector{<:Real},
    coups::Vector{<:Real},
    sitenumbers::Vector{<:Int},
)
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumbers)
        ℓ += freqs[j] * gkslcommutator("N", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumbers, 2, 1))
        jws = jwstring(; start=site1, stop=site2)
        ℓ +=
            coups[j + 1] * (
                gkslcommutator("σ+", site1, jws..., "σ-", site2) +
                gkslcommutator("σ-", site1, jws..., "σ+", site2)
            )
    end
    return ℓ
end

function spin_chain(
    ::SiteType"vElectron",
    freqs::Vector{<:Real},
    coups::Vector{<:Real},
    sitenumbers::Vector{<:Int},
)
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumbers)
        ℓ += freqs[j] * gkslcommutator("Ntot", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumbers, 2, 1))
        jws = jwstring(; start=site1, stop=site2)
        ℓ += +coups[j + 1] * gkslcommutator("Aup†F", site1, jws..., "Aup", site2)
        ℓ += -coups[j + 1] * gkslcommutator("AupF", site1, jws..., "Aup†", site2)
        ℓ += +coups[j + 1] * gkslcommutator("Adn†", site1, jws..., "FAdn", site2)
        ℓ += -coups[j + 1] * gkslcommutator("Adn", site1, jws..., "FAdn†", site2)
    end
    return ℓ
end

############################################################################################

function spin_chain_adjoint(::SiteType, ::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Int})
    return nothing
end

"""
    spin_chain_adjoint(
        freqs::Vector{<:Real},
        coups::Vector{<:Real},
        sites::Vector{<:Index}
    )

Return an OpSum object encoding the adjoint of the Hamiltonian part ``-i[H, –]`` of the
GKSL equation for a spin chain of frequencies `freqs` and coupling constants `coups`, on
`sites`.
"""
function spin_chain_adjoint(
    freqs::Vector{<:Real}, coups::Vector{<:Real}, sites::Vector{<:Index}
)
    stypes = sitetypes(first(sites))
    for st in stypes
        # Check if all sites have this type (otherwise skip this tag).
        if all(i -> st in sitetypes(i), sites)
            # If the type is shared, then try calling the function with it.
            ℓ = spin_chain_adjoint(st, freqs, coups, sitenumber.(sites))
            # If the result is something, return that result.
            if !isnothing(ℓ)
                return ℓ
            end
        end
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"spin_chain_adjoint\" function not found for " *
            "Index tags $(tags(sites[1]))",
        ),
    )
end

# In the GKSL equation, the adjoint of -i[H,-] is simply i[H,-].
function spin_chain_adjoint(
    st::SiteType"vS=1/2",
    freqs::Vector{<:Real},
    coups::Vector{<:Real},
    sitenumbers::Vector{<:Int},
)
    return -spin_chain(st, freqs, coups, sitenumbers)
end

function spin_chain_adjoint(
    st::SiteType"vElectron",
    freqs::Vector{<:Real},
    coups::Vector{<:Real},
    sitenumbers::Vector{<:Int},
)
    return -spin_chain(st, freqs, coups, sitenumbers)
end

const spin_chain′ = spin_chain_adjoint

############################################################################################

exchange_interaction(::SiteType, ::Int, ::Int) = nothing

"""
    exchange_interaction(s1::Index, s2::Index)

Return an OpSum object encoding the Hamiltonian part ``-i[H, –]`` of an exchange interaction
between sites `s1` and `s2` term in a GKSL equation.
"""
function exchange_interaction(s1::Index, s2::Index)
    stypes1 = sitetypes(s1)
    stypes2 = sitetypes(s2)
    for st in stypes1
        # Check if all sites have this type (otherwise skip this tag).
        if st in stypes2
            # If the type is shared, then try calling the function with it.
            ℓ = exchange_interaction(st, sitenumber(s1), sitenumber(s2))
            # If the result is something, return that result.
            if !isnothing(ℓ)
                return ℓ
            end
        end
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"exchange_interaction\" function not found for " *
            "Index tags $(tags(s1))",
        ),
    )
end

function exchange_interaction(::SiteType"Fermion", site1::Int, site2::Int)
    h = OpSum()

    h += "c†", site1, "c", site2
    h += "c†", site2, "c", site1

    return h
end

function exchange_interaction(::SiteType"Electron", site1::Int, site2::Int)
    h = OpSum()

    h += "c†↑", site1, "c↑", site2
    h += "c†↑", site2, "c↑", site1

    h += "c†↓", site1, "c↓", site2
    h += "c†↓", site2, "c↓", site1

    return h
end

function exchange_interaction(::SiteType"vS=1/2", site1::Int, site2::Int)
    ℓ = OpSum()
    jws = jwstring(; start=site1, stop=site2)
    ℓ += (
        gkslcommutator("σ+", site1, jws..., "σ-", site2) +
        gkslcommutator("σ-", site1, jws..., "σ+", site2)
    )
    return ℓ
end

function exchange_interaction(::SiteType"vElectron", site1::Int, site2::Int)
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

############################################################################################

exchange_interaction_adjoint(::SiteType, ::Int, ::Int) = nothing

"""
    exchange_interaction_adjoint(s1::Index, s2::Index)

Return an OpSum object encoding the adjoint of the Hamiltonian part ``-i[H, –]`` of an
exchange interaction between sites `s1` and `s2` term in a GKSL equation.
"""
function exchange_interaction_adjoint(s1::Index, s2::Index)
    stypes1 = sitetypes(s1)
    stypes2 = sitetypes(s2)
    for st in stypes1
        # Check if all sites have this type (otherwise skip this tag).
        if st in stypes2
            # If the type is shared, then try calling the function with it.
            ℓ = exchange_interaction_adjoint(st, sitenumber(s1), sitenumber(s2))
            # If the result is something, return that result.
            if !isnothing(ℓ)
                return ℓ
            end
        end
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"exchange_interaction_adjoint\" function not found for " *
            "Index tags $(tags(s1))",
        ),
    )
end

# In the GKSL equation, the adjoint of -i[H,-] is simply i[H,-].
function exchange_interaction_adjoint(st::SiteType"vS=1/2", site1::Int, site2::Int)
    return -exchange_interaction(st, site1, site2)
end

function exchange_interaction_adjoint(st::SiteType"vElectron", site1::Int, site2::Int)
    return -exchange_interaction(st, site1, site2)
end

const exchange_interaction′ = exchange_interaction_adjoint
