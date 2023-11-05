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

function exchange_interaction(
    ::SiteType"Electron",
    ::SiteType"Fermion",
    electron_site::Int,
    fermion_site::Int;
    coupling_constant_up=1.0,
    coupling_constant_dn=1.0,
)
    h = OpSum()

    jws = jwstring(; start=electron_site, stop=fermion_site)

    h += coupling_constant_up, "a†↑ * F↓", electron_site, jws..., "a", fermion_site
    h += conj(coupling_constant_up), "a↑ * F↓", electron_site, jws..., "a†", fermion_site
    h += coupling_constant_dn, "a†↓", electron_site, jws..., "a", fermion_site
    h += conj(coupling_constant_dn), "a↓", electron_site, jws..., "a†", fermion_site

    return h
end

function exchange_interaction(
    ::SiteType"Fermion",
    ::SiteType"Electron",
    fermion_site::Int;
    electron_site::Int,
    coupling_constant_up=1.0,
    coupling_constant_dn=1.0,
)
    h = OpSum()

    jws = jwstring(; start=electron_site, stop=fermion_site)

    h += coupling_constant_up, "a", fermion_site, jws..., "a†↑", electron_site
    h += conj(coupling_constant_up), "a†", fermion_site, jws..., "a↑", electron_site
    h += coupling_constant_dn, "a", fermion_site, jws..., "F↑ * a†↓", electron_site
    h += conj(coupling_constant_dn), "a†", fermion_site, jws..., "F↑ * a↓", electron_site

    return h
end

function exchange_interaction(
    ::SiteType"vElectron",
    ::SiteType"vFermion",
    electron_site::Int,
    fermion_site::Int,
    coupling_constant_up=1.0,
    coupling_constant_dn=1.0,
)
    jws = jwstring(; start=electron_site, stop=fermion_site)
    return (
        coupling_constant_up * gkslcommutator("Aup†F", dot_site, jws..., "a", other_site) -
        conj(coupling_constant_up) *
        gkslcommutator("AupF", dot_site, jws..., "a†", other_site) +
        coupling_constant_dn * gkslcommutator("Adn†", dot_site, jws..., "a", other_site) -
        conj(coupling_constant_dn) *
        gkslcommutator("Adn", dot_site, jws..., "a†", other_site)
    )
end

function exchange_interaction(
    ::SiteType"vFermion",
    ::SiteType"vElectron",
    fermion_site::Int,
    electron_site::Int,
    coupling_constant_up=1.0,
    coupling_constant_dn=1.0,
)
    jws = jwstring(; start=electron_site, stop=fermion_site)
    return (
        conj(coupling_constant_up) *
        gkslcommutator("a", fermion_site, jws..., "Aup†", electron_site) -
        coupling_constant_up *
        gkslcommutator("a†", fermion_site, jws..., "Aup", electron_site) +
        conj(coupling_constant_dn) *
        gkslcommutator("a", fermion_site, jws..., "FAdn†", electron_site) -
        coupling_constant_dn *
        gkslcommutator("a†", fermion_site, jws..., "FAdn", electron_site)
    )
    # TODO: creare "FdnAdn" e simili
    # coupling_constant_up, "a", fermion_site, jws..., "a†↑", electron_site
    # conj(coupling_constant_up), "a†", fermion_site, jws..., "a↑", electron_site
    # coupling_constant_dn, "a", fermion_site, jws..., "F↑ * a†↓", electron_site
    # conj(coupling_constant_dn), "a†", fermion_site, jws..., "F↑ * a↓", electron_site
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

function exchange_interaction_adjoint(
    st1::SiteType"vElectron", st2::SiteType"vFermion", site1::Int, site2::Int; kwargs...
)
    return -exchange_interaction(st1, st2, site1, site2; kwargs...)
end

function exchange_interaction_adjoint(
    st1::SiteType"vFermion", st2::SiteType"vElectron", site1::Int, site2::Int; kwargs...
)
    return -exchange_interaction(st1, st2, site1, site2; kwargs...)
end

const exchange_interaction′ = exchange_interaction_adjoint
