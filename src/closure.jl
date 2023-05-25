export Closure, closure
export freq,
       innercoup,
       outercoup,
       damp,
       freqs,
       innercoups,
       outercoups,
       damps
export closure_op, filled_closure_op

"""
    Closure(ω::Vector{<:Real}, γ::Vector{<:Real}, g::Vector{<:Real}, ζ::Vector{<:Complex})

Closure is an aggregate type which stores the parameters needed to describe the set of
pseudomodes that make up a Markovian closure.
"""
struct Closure
    frequency
    innercoupling
    outercoupling
    damping
    function Closure(ω::Vector{<:Real}, γ::Vector{<:Real}, g::Vector{<:Real}, ζ::Vector{<:Complex})
        if (
            length(ω) - 1 != length(g) ||
            length(ω) != length(γ) ||
            length(ω) != length(ζ)
        )
            error("Lengths of input parameters do not match.")
        end
        return new(ω, g, ζ, γ)
    end
end

"""
    function closure(
        Ω::Real, K::Real,
        α::Matrix{<:Real}, β::Matrix{<:Real}, w::Matrix{<:Real}
    )

Construct a Closure object with asymptotic frequency `Ω` and coupling `K`, with the
universal parameter set given by `α`, `β` and `w`, which are assumed to be two-column
with the real parts of the parameters in the first one and the imaginary parts in the
second.
"""
function closure(
        Ω::Real, K::Real,
        α::Matrix{<:Real}, β::Matrix{<:Real}, w::Matrix{<:Real}
)
    return closure(
        Ω, K, α[:, 1] .+ im .* α[:, 2], β[:, 1] .+ im .* β[:, 2], w[:, 1] .+ im .* w[:, 2]
    )
end

"""
    function closure(
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

freqs(mc::Closure) = mc.frequency
innercoups(mc::Closure) = mc.innercoupling
outercoups(mc::Closure) = mc.outercoupling
damps(mc::Closure) = mc.damping

freq(mc::Closure, j::Int) = mc.frequency[j]
innercoup(mc::Closure, j::Int) = mc.innercoupling[j]
outercoup(mc::Closure, j::Int) = mc.outercoupling[j]
damp(mc::Closure, j::Int) = mc.damping[j]

"""
    closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`.
This closure replaces a chain starting from an empty state.
"""
function closure_op(
    ::SiteType"vS=1/2", mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int
)
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumber.(sites))
        ℓ += freq(mc, j) * gkslcommutator("N", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumber.(sites), 2, 1))
        jws = jwstring(; start=site1, stop=site2)
        ℓ +=
            innercoup(mc, j) * (
                gkslcommutator("σ-", site1, jws..., "σ+", site2) +
                gkslcommutator("σ+", site1, jws..., "σ-", site2)
            )
    end
    for (j, site) in enumerate(sitenumber.(sites))
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, jws..., "σ-", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("σ-", chain_edge_site, jws..., "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumber.(sites))
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "σ-⋅ * ⋅σ+"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function closure_op(::SiteType"vElectron", mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumber.(sites))
        ℓ += freq(mc, j) * gkslcommutator("Ntot", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumber.(sites), 2, 1))
        # Add as many "F" string operator as needed.
        jws = jwstring(; start=site1, stop=site2)
        ℓ +=
            innercoup(mc, j) * (
                gkslcommutator("Aup†F", site1, jws..., "Aup", site2) -
                gkslcommutator("AupF", site1, jws..., "Aup†", site2) +
                gkslcommutator("Adn†", site1, jws..., "FAdn", site2) -
                gkslcommutator("Adn", site1, jws..., "FAdn†", site2)
            )
    end
    for (j, site) in enumerate(sitenumber.(sites))
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
    for (j, site) in enumerate(sitenumber.(sites))
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

"""
    filled_closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`.
This closure replaces a chain starting from a completely filled state.
"""
function filled_closure_op(
    ::SiteType"vS=1/2", mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int
)
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumber.(sites))
        ℓ += freq(mc, j) * gkslcommutator("N", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumber.(sites), 2, 1))
        jws = jwstring(; start=site1, stop=site2)
        ℓ +=
            innercoup(mc, j) * (
                gkslcommutator("σ-", site1, jws..., "σ+", site2) +
                gkslcommutator("σ+", site1, jws..., "σ-", site2)
            )
    end
    for (j, site) in enumerate(sitenumber.(sites))
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, jws..., "σ-", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("σ-", chain_edge_site, jws..., "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumber.(sites))
        # a† ρ a
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "σ+⋅ * ⋅σ-"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # -0.5 (a a† ρ + ρ a a†)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_closure_op(::SiteType"vElectron", mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)
    ℓ = OpSum()
    for (j, site) in enumerate(sitenumber.(sites))
        ℓ += freq(mc, j) * gkslcommutator("Ntot", site)
    end
    for (j, (site1, site2)) in enumerate(partition(sitenumber.(sites), 2, 1))
        # Add as many "F" string operator as needed.
        jws = jwstring(; start=site1, stop=site2)
        ℓ +=
            innercoup(mc, j) * (
                gkslcommutator("Aup†F", site1, jws..., "Aup", site2) -
                gkslcommutator("AupF", site1, jws..., "Aup†", site2) +
                gkslcommutator("Adn†", site1, jws..., "FAdn", site2) -
                gkslcommutator("Adn", site1, jws..., "FAdn†", site2)
            )
    end
    for (j, site) in enumerate(sitenumber.(sites))
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
    for (j, site) in enumerate(sitenumber.(sites))
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
