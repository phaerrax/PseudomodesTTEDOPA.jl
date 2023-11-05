export exchange_interaction, exchange_interaction_adjoint, exchange_interaction′

exchange_interaction(::SiteType, ::SiteType, ::Int, ::Int; kwargs...) = nothing

"""
    exchange_interaction(s1::Index, s2::Index; kwargs...)

Return an OpSum object encoding the Hamiltonian part ``-i[H, –]`` of an exchange interaction
between sites `s1` and `s2` term in a GKSL equation.
"""
function exchange_interaction(s1::Index, s2::Index; kwargs...)
    for (st1, st2) in zip(sitetypes(s1), sitetypes(s2))
        ℓ = exchange_interaction(st1, st2, sitenumber(s1), sitenumber(s2); kwargs...)
        # If the result is something, return that result.
        if !isnothing(ℓ)
            return ℓ
        end
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"exchange_interaction\" function not found for " *
            "Index tags $(tags(s1)) and $(tags(s2))",
        ),
    )
end

function exchange_interaction(
    ::SiteType"Fermion", ::SiteType"Fermion", site1::Int, site2::Int; coupling_constant=1.0
)
    h = OpSum()

    h += coupling_constant, "c†", site1, "c", site2
    h += conj(coupling_constant), "c†", site2, "c", site1

    return h
end

function exchange_interaction(
    ::SiteType"vFermion",
    ::SiteType"vFermion",
    site1::Int,
    site2::Int,
    coupling_constant=1.0,
)
    ℓ = OpSum()
    jws = jwstring(; start=site1, stop=site2)
    ℓ += (
        coupling_constant * gkslcommutator("A†", site1, jws..., "A", site2) +
        conj(coupling_constant) * gkslcommutator("A", site1, jws..., "A†", site2)
    )
    return ℓ
end

function exchange_interaction(
    ::SiteType"vS=1/2", ::SiteType"vS=1/2", site1::Int, site2::Int, coupling_constant=1.0
)
    ℓ = OpSum()
    ℓ += (
        coupling_constant * gkslcommutator("σ+", site1, "σ-", site2) +
        conj(coupling_constant) * gkslcommutator("σ-", site1, "σ+", site2)
    )
    return ℓ
end

function exchange_interaction(
    ::SiteType"Electron",
    ::SiteType"Electron",
    site1::Int,
    site2::Int;
    coupling_constant_up=1.0,
    coupling_constant_dn=1.0,
)
    h = OpSum()

    h += coupling_constant_up, "c†↑", site1, "c↑", site2
    h += conj(coupling_constant_up), "c†↑", site2, "c↑", site1

    h += coupling_constant_dn, "c†↓", site1, "c↓", site2
    h += conj(coupling_constant_dn), "c†↓", site2, "c↓", site1

    return h
end

function exchange_interaction(
    ::SiteType"vElectron",
    ::SiteType"vElectron",
    site1::Int,
    site2::Int;
    coupling_constant_up=1.0,
    coupling_constant_dn=1.0,
)
    # c↑ᵢ† c↑ᵢ₊₁ + c↑ᵢ₊₁† c↑ᵢ + c↓ᵢ† c↓ᵢ₊₁ + c↓ᵢ₊₁† c↓ᵢ =
    # a↑ᵢ† Fᵢ a↑ᵢ₊₁ - a↑ᵢ Fᵢ a↑ᵢ₊₁† + a↓ᵢ† Fᵢ₊₁ a↓ᵢ₊₁ - a↓ᵢ Fᵢ₊₁ a↓ᵢ₊₁†
    ℓ = OpSum()
    jws = jwstring(; start=site1, stop=site2)
    ℓ += (
        coupling_constant_up * gkslcommutator("Aup†F", site1, jws..., "Aup", site2) -
        conj(coupling_constant_up) * gkslcommutator("AupF", site1, jws..., "Aup†", site2) +
        coupling_constant_dn * gkslcommutator("Adn†", site1, jws..., "FAdn", site2) -
        conj(coupling_constant_dn) * gkslcommutator("Adn", site1, jws..., "FAdn†", site2)
    )
    return ℓ
end

############################################################################################

exchange_interaction_adjoint(::SiteType, ::SiteType, ::Int, ::Int; kwargs...) = nothing

"""
    exchange_interaction_adjoint(s1::Index, s2::Index; kwargs...)

Return an OpSum object encoding the adjoint of the Hamiltonian part ``-i[H, –]`` of an
exchange interaction between sites `s1` and `s2` term in a GKSL equation.
"""
function exchange_interaction_adjoint(s1::Index, s2::Index; kwargs...)
    for (st1, st2) in zip(sitetypes(s1), sitetypes(s2))
        ℓ = exchange_interaction_adjoint(
            st1, st2, sitenumber(s1), sitenumber(s2); kwargs...
        )
        # If the result is something, return that result.
        if !isnothing(ℓ)
            return ℓ
        end
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"exchange_interaction_adjoint\" function not found for " *
            "Index tags $(tags(s1)) and $(tags(s2))",
        ),
    )
end

# In the GKSL equation, the adjoint of -i[H,-] is simply i[H,-].

function exchange_interaction_adjoint(
    st1::SiteType"vFermion", st2::SiteType"vFermion", site1::Int, site2::Int; kwargs...
)
    return -exchange_interaction(st1, st2, site1, site2; kwargs...)
end

function exchange_interaction_adjoint(
    st1::SiteType"vS=1/2", st2::SiteType"vS=1/2", site1::Int, site2::Int; kwargs...
)
    return -exchange_interaction(st1, st2, site1, site2; kwargs...)
end

function exchange_interaction_adjoint(
    st1::SiteType"vElectron", st2::SiteType"vElectron", site1::Int, site2::Int; kwargs...
)
    return -exchange_interaction(st1, st2, site1, site2; kwargs...)
end

const exchange_interaction′ = exchange_interaction_adjoint
