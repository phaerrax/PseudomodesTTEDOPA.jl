export Closure, closure
export freq, innercoup, outercoup, damp, freqs, innercoups, outercoups, damps
export closure_op, closure_op_adjoint, closure_op′
export filled_closure_op, filled_closure_op_adjoint, filled_closure_op′

"""
A Closure object stores the parameters that define a Markovian closure, providing different
ways to input them.
"""
struct Closure
    frequency
    innercoupling
    outercoupling
    damping
    @doc"""
        Closure(ω::Vector{<:Real}, γ::Vector{<:Real}, g::Vector{<:Real}, ζ::Vector{<:Complex})

    Closure is an aggregate type which stores the parameters needed to describe the set of
    pseudomodes that make up a Markovian closure.
    """
    function Closure(
        ω::Vector{<:Real}, γ::Vector{<:Real}, g::Vector{<:Real}, ζ::Vector{<:Complex}
    )
        if (length(ω) - 1 != length(g) || length(ω) != length(γ) || length(ω) != length(ζ))
            error("Lengths of input parameters do not match.")
        end
        return new(ω, g, ζ, γ)
    end
end

"""
    closure(
        Ω::Real, K::Real,
        α::Matrix{<:Real}, β::Matrix{<:Real}, w::Matrix{<:Real}
    )

Construct a Closure object with asymptotic frequency `Ω` and coupling `K`, with the
universal parameter set given by `α`, `β` and `w`, each of which is assumed to be a
two-column matrix with the real parts of the parameters in the first column and the
imaginary parts in the second one.
"""
function closure(Ω::Real, K::Real, α::Matrix{<:Real}, β::Matrix{<:Real}, w::Matrix{<:Real})
    return closure(
        Ω, K, α[:, 1] .+ im .* α[:, 2], β[:, 1] .+ im .* β[:, 2], w[:, 1] .+ im .* w[:, 2]
    )
end

"""
    closure(
        Ω::Real, K::Real,
        α::Vector{<:Complex}, β::Vector{<:Complex}, w::Vector{<:Complex}
    )

Construct a Closure object with asymptotic frequency `Ω` and coupling `K`, with the
universal parameter set given by `α`, `β` and `w`, which are assumed to be vectors of
complex numbers.
"""
function closure(
    Ω::Real, K::Real, α::Vector{<:Complex}, β::Vector{<:Complex}, w::Vector{<:Complex}
)
    frequency = @. Ω - 2K * imag(α)
    damping = @. -4K * real(α)
    innercoupling = @. -2K * imag(β)
    outercoupling = @. K * w
    return Closure(frequency, damping, innercoupling, outercoupling)
end

"""
    length(mc::Closure)

Return the length of a Markovian closure, i.e. the number of pseudomodes.
"""
Base.length(mc::Closure) = Base.length(mc.frequency)

"""
    freqs(mc::Closure)

Return the list of frequencies of the pseudomodes in the given Markovian closure.
"""
freqs(mc::Closure) = mc.frequency

"""
    innercoups(mc::Closure)

Return the list of inner couplings of the pseudomodes in the Markovian closure `mc`,
i.e. the coupling constant between the pseudomodes.
"""
innercoups(mc::Closure) = mc.innercoupling

"""
    outercoups(mc::Closure)

Return the list of outer couplings of the pseudomodes in the Markovian closure `mc`,
i.e. the coupling constants with the last point of the TEDOPA chain.
"""
outercoups(mc::Closure) = mc.outercoupling

"""
    damps(mc::Closure)

Return the list of damping coefficients of the pseudomodes in the Markovian closure `mc`.
"""
damps(mc::Closure) = mc.damping

"""
    freq(mc::Closure, j::Int)

Return the frequency of the `j`-th pseudomode of the Markovian closure `mc`.
"""
freq(mc::Closure, j::Int) = mc.frequency[j]

"""
    innercoup(mc::Closure, j::Int)

Return the coupling coefficient between the `j`-th and the `j+1`-th pseudomodes of the
Markovian closure `mc`.
"""
innercoup(mc::Closure, j::Int) = mc.innercoupling[j]

"""
    outercoup(mc::Closure, j::Int)

Return the coupling coefficient of the `j`-th pseudomode in the Markovian closure `mc`,
with the last point of the TEDOPA chain.
"""
outercoup(mc::Closure, j::Int) = mc.outercoupling[j]

"""
    damp(mc::Closure, j::Int)

Return the damping coefficient of the `j`-th pseudomode in the Markovian closure `mc`.
"""
damp(mc::Closure, j::Int) = mc.damping[j]

############################################################################################

closure_op(::SiteType, ::Closure, ::Vector{<:Int}, ::Int) = nothing

"""
    closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`.
This closure replaces a chain starting from an empty state.
"""
function closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)
    @assert length(mc) == length(sites)
    stypes = sitetypes(first(sites))
    for st in stypes
        # Check if all sites have this type (otherwise skip this tag).
        if all(i -> st in sitetypes(i), sites)
            # If the type is shared, then try calling the function with it.
            ℓ = closure_op(st, mc, sitenumber.(sites), chain_edge_site)
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
            "Overload of \"closure_op\" function not found for Index tags $(tags(sites[1]))"
        ),
    )
end

function closure_op(
    st::SiteType"vFermion", mc::Closure, sitenumbers::Vector{<:Int}, chain_edge_site::Int
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("A†", chain_edge_site, jws..., "A", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("A", chain_edge_site, jws..., "A†", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "A⋅ * ⋅A†"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function closure_op(
    st::SiteType"vS=1/2", mc::Closure, sitenumbers::Vector{<:Int}, chain_edge_site::Int
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, "σ-", site) +
            conj(outercoup(mc, j)) * gkslcommutator("σ-", chain_edge_site, "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        ℓ += damp(mc, j), "σ-⋅ * ⋅σ+", site
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function closure_op(
    st::SiteType"vElectron", mc::Closure, sitenumbers::Vector{<:Int}, chain_edge_site::Int
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        # c↑ᵢ† c↑ᵢ₊ₙ = a↑ᵢ† Fᵢ Fᵢ₊₁ ⋯ Fᵢ₊ₙ₋₁ a↑ᵢ₊ₙ
        # c↑ᵢ₊ₙ† c↑ᵢ = -a↑ᵢ Fᵢ Fᵢ₊₁ ⋯ Fᵢ₊ₙ₋₁ a↑ᵢ₊ₙ†
        # c↓ᵢ† c↓ᵢ₊ₙ = a↓ᵢ† Fᵢ₊₁ Fᵢ₊₂ ⋯ Fᵢ₊ₙ a↓ᵢ₊ₙ
        # c↓ᵢ₊ₙ† c↓ᵢ = -a↓ᵢ Fᵢ₊₁ Fᵢ₊₂ ⋯ Fᵢ₊ₙ a↓ᵢ₊ₙ†

        jws = jwstring(; start=chain_edge_site, stop=site)
        # ζⱼ c↑₀† c↑ⱼ (0 = chain edge, j = pseudomode)
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Aup†F", chain_edge_site, jws..., "Aup", site)
        # conj(ζⱼ) c↑ⱼ† c↑₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("AupF", chain_edge_site, jws..., "Aup†", site)
        # ζⱼ c↓₀† c↓ⱼ
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Adn†", chain_edge_site, jws..., "FAdn", site)
        # conj(ζⱼ) c↓ⱼ† c↓₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("Adn", chain_edge_site, jws..., "FAdn†", site)
    end

    # Dissipative part
    for (j, site) in enumerate(sitenumbers)
        # c↑ₖ ρ c↑ₖ† = F₁ ⋯ Fₖ₋₁ a↑ₖ ρ a↑ₖ† Fₖ₋₁ ⋯ F₁
        # Remember that:
        # • Fⱼ = (1 - 2 N↑ₖ) (1 - 2 N↓ₖ);
        # • Fⱼ and aₛ,ₖ commute only on different sites;
        # • {a↓ₖ, a↓ₖ†} = {a↑ₖ, a↑ₖ†} = 1;
        # • Fₖ anticommutes with a↓ₖ, a↓ₖ†, a↑ₖ and a↑ₖ†.
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Aup⋅ * ⋅Aup†"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # c↓ₖ ρ c↓ₖ† = F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ ρ a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "FAdn⋅ * ⋅Adn†F"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)

        # -½ (c↑ₖ† c↑ₖ ρ + ρ c↑ₖ† c↑ₖ) = -½ (a↑ₖ† a↑ₖ ρ + ρ a↑ₖ† a↑ₖ)
        ℓ += -0.5damp(mc, j), "Nup⋅", site
        ℓ += -0.5damp(mc, j), "⋅Nup", site
        # -½ (c↓ₖ† c↓ₖ ρ + ρ c↓ₖ† c↓ₖ) = -½ (a↓ₖ† a↓ₖ ρ + ρ a↓ₖ† a↓ₖ)
        ℓ += -0.5damp(mc, j), "Ndn⋅", site
        ℓ += -0.5damp(mc, j), "⋅Ndn", site
    end
    return ℓ
end

############################################################################################

closure_op_adjoint(::SiteType, ::Closure, ::Vector{<:Int}, ::Int, ::Int) = nothing

"""
    closure_op_adjoint(
        mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int, gradefactor::Int
    )

Return an OpSum object encoding the adjoint of the Markovian closure operators with
parameters given by `mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on
site `chain_edge_site`. The integer `gradefactor` is the parity (1: even, -1: odd) of the
operator subject to the time evolution.
This closure replaces a chain starting from an empty state.
"""
function closure_op_adjoint(
    mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int, gradefactor::Int
)
    @assert length(mc) == length(sites)
    stypes = sitetypes(first(sites))
    for st in stypes
        # Check if all sites have this type (otherwise skip this tag).
        if all(i -> st in sitetypes(i), sites)
            # If the type is shared, then try calling the function with it.
            ℓ = closure_op_adjoint(st, mc, sitenumber.(sites), chain_edge_site, gradefactor)
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
            "Overload of \"closure_op\" function not found for Index tags $(tags(sites[1]))"
        ),
    )
end

function closure_op_adjoint(
    st::SiteType"vFermion",
    mc::Closure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int,
)
    ℓ = spin_chain_adjoint(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += -(
            outercoup(mc, j) * gkslcommutator("A†", chain_edge_site, jws..., "A", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("A", chain_edge_site, jws..., "A†", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "A†⋅ * ⋅A"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function closure_op_adjoint(
    st::SiteType"vS=1/2",
    mc::Closure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int,
)
    ℓ = spin_chain_adjoint(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        ℓ += -(
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, "σ-", site) +
            conj(outercoup(mc, j)) * gkslcommutator("σ-", chain_edge_site, "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a† ρ a
        ℓ += gradefactor * damp(mc, j), "σ+⋅ * ⋅σ-", site
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function closure_op_adjoint(
    st::SiteType"vElectron",
    mc::Closure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int,
)
    ℓ = spin_chain_adjoint(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        # c↑ᵢ† c↑ᵢ₊ₙ = a↑ᵢ† Fᵢ Fᵢ₊₁ ⋯ Fᵢ₊ₙ₋₁ a↑ᵢ₊ₙ
        # c↑ᵢ₊ₙ† c↑ᵢ = -a↑ᵢ Fᵢ Fᵢ₊₁ ⋯ Fᵢ₊ₙ₋₁ a↑ᵢ₊ₙ†
        # c↓ᵢ† c↓ᵢ₊ₙ = a↓ᵢ† Fᵢ₊₁ Fᵢ₊₂ ⋯ Fᵢ₊ₙ a↓ᵢ₊ₙ
        # c↓ᵢ₊ₙ† c↓ᵢ = -a↓ᵢ Fᵢ₊₁ Fᵢ₊₂ ⋯ Fᵢ₊ₙ a↓ᵢ₊ₙ†

        jws = jwstring(; start=chain_edge_site, stop=site)
        # ζⱼ c↑₀† c↑ⱼ (0 = chain edge, j = pseudomode)
        ℓ +=
            -outercoup(mc, j) *
            gkslcommutator("Aup†F", chain_edge_site, jws..., "Aup", site)
        # conj(ζⱼ) c↑ⱼ† c↑₀
        ℓ +=
            conj(outercoup(mc, j)) *
            gkslcommutator("AupF", chain_edge_site, jws..., "Aup†", site)
        # ζⱼ c↓₀† c↓ⱼ
        ℓ +=
            -outercoup(mc, j) *
            gkslcommutator("Adn†", chain_edge_site, jws..., "FAdn", site)
        # conj(ζⱼ) c↓ⱼ† c↓₀
        ℓ +=
            conj(outercoup(mc, j)) *
            gkslcommutator("Adn", chain_edge_site, jws..., "FAdn†", site)
    end

    # Dissipative part
    for (j, site) in enumerate(sitenumbers)
        # Remember that:
        # • Fⱼ = (1 - 2 N↑ₖ) (1 - 2 N↓ₖ);
        # • Fⱼ and aₛ,ₖ commute only on different sites;
        # • {a↓ₖ, a↓ₖ†} = {a↑ₖ, a↑ₖ†} = 1;
        # • Fₖ anticommutes with a↓ₖ, a↓ₖ†, a↑ₖ and a↑ₖ†.

        # c↑ₖ† X c↑ₖ = a↑ₖ† Fₖ₋₁ ⋯ F₁ X F₁ ⋯ Fₖ₋₁ a↑ₖ
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Aup†⋅ * ⋅Aup"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )
        # c↓ₖ† X c↓ₖ = a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁ X F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Adn†F⋅ * ⋅FAdn"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )

        # -½ (c↑ₖ† c↑ₖ ρ + ρ c↑ₖ† c↑ₖ) = -½ (a↑ₖ† a↑ₖ ρ + ρ a↑ₖ† a↑ₖ)
        ℓ += -0.5damp(mc, j), "Nup⋅", site
        ℓ += -0.5damp(mc, j), "⋅Nup", site
        # -½ (c↓ₖ† c↓ₖ ρ + ρ c↓ₖ† c↓ₖ) = -½ (a↓ₖ† a↓ₖ ρ + ρ a↓ₖ† a↓ₖ)
        ℓ += -0.5damp(mc, j), "Ndn⋅", site
        ℓ += -0.5damp(mc, j), "⋅Ndn", site
    end
    return ℓ
end

const closure_op′ = closure_op_adjoint

############################################################################################

filled_closure_op(::SiteType, ::Closure, ::Vector{<:Int}, ::Int) = nothing

"""
    filled_closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`.
This closure replaces a chain starting from a completely filled state.
"""
function filled_closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)
    @assert length(mc) == length(sites)
    stypes = sitetypes(first(sites))
    for st in stypes
        # Check if all sites have this type (otherwise skip this tag).
        if all(i -> st in sitetypes(i), sites)
            # If the type is shared, then try calling the function with it.
            ℓ = filled_closure_op(st, mc, sitenumber.(sites), chain_edge_site)
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
            "Overload of \"filled_closure_op\" function not found for " *
            "Index tags $(tags(sites[1]))",
        ),
    )
end

function filled_closure_op(
    st::SiteType"vFermion", mc::Closure, sitenumbers::Vector{<:Int}, chain_edge_site::Int
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("A†", chain_edge_site, jws..., "A", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("A", chain_edge_site, jws..., "A†", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "A†⋅ * ⋅A"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # -0.5 (a a† ρ + ρ a a†)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_closure_op(
    st::SiteType"vS=1/2", mc::Closure, sitenumbers::Vector{<:Int}, chain_edge_site::Int
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, "σ-", site) +
            conj(outercoup(mc, j)) * gkslcommutator("σ-", chain_edge_site, "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a† ρ a
        ℓ += damp(mc, j), "σ+⋅ * ⋅σ-", site
        # -0.5 (a a† ρ + ρ a a†)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_closure_op(
    st::SiteType"vElectron", mc::Closure, sitenumbers::Vector{<:Int}, chain_edge_site::Int
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        # ζⱼ c↑₀† c↑ⱼ (0 = chain edge, j = pseudomode)
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Aup†F", chain_edge_site, jws..., "Aup", site)
        # conj(ζⱼ) c↑ⱼ† c↑₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("AupF", chain_edge_site, jws..., "Aup†", site)
        # ζⱼ c↓₀† c↓ⱼ
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Adn†", chain_edge_site, jws..., "FAdn", site)
        # conj(ζⱼ) c↓ⱼ† c↓₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("Adn", chain_edge_site, jws..., "FAdn†", site)
    end

    # Dissipative part
    for (j, site) in enumerate(sitenumbers)
        # c↑ₖ† ρ c↑ₖ = a↑ₖ† Fₖ₋₁ ⋯ F₁ ρ F₁ ⋯ Fₖ₋₁ a↑ₖ
        # Remember that:
        # • Fⱼ = (1 - 2 N↑ₖ) (1 - 2 N↓ₖ);
        # • Fⱼ and aₛ,ₖ commute only on different sites;
        # • {a↓ₖ, a↓ₖ†} = {a↑ₖ, a↑ₖ†} = 1;
        # • Fₖ anticommutes with a↓ₖ, a↓ₖ†, a↑ₖ and a↑ₖ†.
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Aup†⋅ * ⋅Aup"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # c↓ₖ† ρ c↓ₖ = a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁ ρ F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Adn†F⋅ * ⋅FAdn"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)

        # c↑ₖ c↑ₖ† = F₁ ⋯ Fₖ₋₁ a↑ₖ a↑ₖ† Fₖ₋₁ ⋯ F₁ =
        #          = F₁² ⋯ Fₖ₋₁² a↑ₖ a↑ₖ† =
        #          = a↑ₖ a↑ₖ† =
        #          = 1 - N↑ₖ
        ℓ += -damp(mc, j), "Id", site
        ℓ += 0.5damp(mc, j), "Nup⋅", site
        ℓ += 0.5damp(mc, j), "⋅Nup", site
        # c↓ₖ c↓ₖ† = F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁ =
        #          = F₁² ⋯ Fₖ₋₁² Fₖa↓ₖ a↓ₖ†Fₖ =
        #          = Fₖa↓ₖ a↓ₖ†Fₖ =
        #          = Fₖ² a↓ₖ a↓ₖ† =
        #          = 1 - N↓ₖ
        ℓ += -damp(mc, j), "Id", site
        ℓ += 0.5damp(mc, j), "Ndn⋅", site
        ℓ += 0.5damp(mc, j), "⋅Ndn", site
    end
    return ℓ
end

############################################################################################

filled_closure_op_adjoint(::SiteType, ::Closure, ::Vector{<:Int}, ::Int, ::Int) = nothing

"""
    filled_closure_op_adjoint(
        mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int, gradefactor::Int
    )

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`. The integer `gradefactor` is the parity (1: even, -1: odd) of the
operator subject to the time evolution.
This closure replaces a chain starting from a completely filled state.
"""
function filled_closure_op_adjoint(
    mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int, gradefactor::Int
)
    @assert length(mc) == length(sites)
    stypes = sitetypes(first(sites))
    for st in stypes
        # Check if all sites have this type (otherwise skip this tag).
        if all(i -> st in sitetypes(i), sites)
            # If the type is shared, then try calling the function with it.
            ℓ = filled_closure_op_adjoint(
                st, mc, sitenumber.(sites), chain_edge_site, gradefactor
            )
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
            "Overload of \"filled_closure_op\" function not found for " *
            "Index tags $(tags(sites[1]))",
        ),
    )
end

function filled_closure_op_adjoint(
    st::SiteType"vFermion",
    mc::Closure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int,
)
    ℓ = adjoint_spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += -(
            outercoup(mc, j) * gkslcommutator("A†", chain_edge_site, jws..., "A", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("A", chain_edge_site, jws..., "A†", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "A⋅ * ⋅A†"]
        ℓ += (gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # -0.5 (a a† ρ + ρ a a†) = -0.5 (ρ - a† a ρ + ρ a† a)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_closure_op_adjoint(
    st::SiteType"vS=1/2",
    mc::Closure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int,
)
    ℓ = adjoint_spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        ℓ += -(
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, "σ-", site) +
            conj(outercoup(mc, j)) * gkslcommutator("σ-", chain_edge_site, "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        ℓ += damp(mc, j), "σ-⋅ * ⋅σ+", site
        # -0.5 (a a† ρ + ρ a a†) = -0.5 (2ρ - a† a ρ + ρ a† a)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_closure_op_adjoint(
    st::SiteType"vElectron",
    mc::Closure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int,
)
    ℓ = adjoint_spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        # ζⱼ c↑₀† c↑ⱼ (0 = chain edge, j = pseudomode)
        ℓ +=
            -outercoup(mc, j) *
            gkslcommutator("Aup†F", chain_edge_site, jws..., "Aup", site)
        # conj(ζⱼ) c↑ⱼ† c↑₀
        ℓ +=
            conj(outercoup(mc, j)) *
            gkslcommutator("AupF", chain_edge_site, jws..., "Aup†", site)
        # ζⱼ c↓₀† c↓ⱼ
        ℓ +=
            -outercoup(mc, j) *
            gkslcommutator("Adn†", chain_edge_site, jws..., "FAdn", site)
        # conj(ζⱼ) c↓ⱼ† c↓₀
        ℓ +=
            conj(outercoup(mc, j)) *
            gkslcommutator("Adn", chain_edge_site, jws..., "FAdn†", site)
    end

    # Dissipative part
    for (j, site) in enumerate(sitenumbers)
        # Remember that:
        # • Fⱼ = (1 - 2 N↑ₖ) (1 - 2 N↓ₖ);
        # • Fⱼ and aₛ,ₖ commute only on different sites;
        # • {a↓ₖ, a↓ₖ†} = {a↑ₖ, a↑ₖ†} = 1;
        # • Fₖ anticommutes with a↓ₖ, a↓ₖ†, a↑ₖ and a↑ₖ†.
        #
        # c↑ₖ X c↑ₖ† = F₁ ⋯ Fₖ₋₁ a↑ₖ X a↑ₖ† Fₖ₋₁ ⋯ F₁
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Aup⋅ * ⋅Aup†"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )
        # c↓ₖ X c↓ₖ† = F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ X a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "FAdn⋅ * ⋅Adn†F"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )

        # c↑ₖ c↑ₖ† = F₁ ⋯ Fₖ₋₁ a↑ₖ a↑ₖ† Fₖ₋₁ ⋯ F₁ =
        #          = F₁² ⋯ Fₖ₋₁² a↑ₖ a↑ₖ† =
        #          = a↑ₖ a↑ₖ† =
        #          = 1 - N↑ₖ
        ℓ += -damp(mc, j), "Id", site
        ℓ += 0.5damp(mc, j), "Nup⋅", site
        ℓ += 0.5damp(mc, j), "⋅Nup", site
        # c↓ₖ c↓ₖ† = F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁ =
        #          = F₁² ⋯ Fₖ₋₁² Fₖa↓ₖ a↓ₖ†Fₖ =
        #          = Fₖa↓ₖ a↓ₖ†Fₖ =
        #          = Fₖ² a↓ₖ a↓ₖ† =
        #          = 1 - N↓ₖ
        ℓ += -damp(mc, j), "Id", site
        ℓ += 0.5damp(mc, j), "Ndn⋅", site
        ℓ += 0.5damp(mc, j), "⋅Ndn", site
    end
    return ℓ
end

const filled_closure_op′ = filled_closure_op_adjoint
