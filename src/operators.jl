@doc raw"""
    twositeoperators(sites::Vector{Index{Int64}}, localcfs::Vector{<:Real}, interactioncfs::Vector{<:Real})

Return a list of ITensor operators representing the two-site-gates
decomposition of the Hamiltonian or Lindbladian of the system.

The element `list[j]` is the operator ``hⱼ,ⱼ₊₁`` (or ``ℓⱼ,ⱼ₊₁``, depending on
the notation) in the expansion
```math
H=\sum_{j=1}^{N}h_{j,j+1}.
```

# Arguments
- `sites::Vector{Index}`: the sites making up the system
- `localcfs::Vector{Index}`: coefficients of one-site operators
- `interactioncfs::Vector{Index}`: coefficients of two-site operators

Index convention: the `j`-th element refers to the ``hⱼ,ⱼ₊₁``/``ℓⱼ,ⱼ₊₁`` term.

# Usage
This function is extremely specialised to our use-case: the function
automatically chooses some operators representing local and interaction terms,
but their OpNames are currently hardcoded (see `localop` and `interactionop`).
"""
function twositeoperators(
    sites::Vector{Index{Int64}}, localcfs::Vector{<:Real}, interactioncfs::Vector{<:Real}
)
    list = ITensor[]
    # The sites at the edge of the chain need to be treated separately since their
    # local coefficients must not be multiplied by 1/2. Instead of writing these
    # particular cases out of the loop, I just multiply them beforehand by 2.
    localcfs[begin] *= 2
    localcfs[end] *= 2
    for j in 1:(length(sites) - 1)
        s1 = sites[j]
        s2 = sites[j + 1]
        h =
            0.5localcfs[j] * localop(s1) * op("Id", s2) +
            0.5localcfs[j + 1] * op("Id", s1) * localop(s2) +
            interactioncfs[j] * interactionop(s1, s2)
        push!(list, h)
    end
    return list
end

"""
    localop(s::Index)

Return the local ITensor operator associated to `s`'s SiteType.

# Usage
Currently this is just a way to simplify code, the user doesn't really have a
say in which operator is used, it is just hardcoded.
"""
function localop(s::Index)
    if SiteType("S=1/2") ∈ sitetypes(s)
        t = 0.5op("σz", s)
    elseif SiteType("Osc") ∈ sitetypes(s)
        t = op("N", s)
    elseif SiteType("vecS=1/2") ∈ sitetypes(s)
        t = 0.5im * (op("σz:Id", s) - op("Id:σz", s))
    elseif SiteType("HvS=1/2") ∈ sitetypes(s)
        t = -0.5im * (op("σz⋅", s) - op("⋅σz", s))
    elseif SiteType("vecOsc") ∈ sitetypes(s)
        t = im * (op("N:Id", s) - op("Id:N", s))
    elseif SiteType("HvOsc") ∈ sitetypes(s)
        t = -im * (op("N⋅", s) - op("⋅N", s))
    else
        throw(DomainError(s, "SiteType non riconosciuto."))
    end
    return t
end

"""
    interactionop(s1::Index, s2::Index)

Return the ITensor local operator associated to the SiteTypes of `s1` and `s2`.

# Usage
Currently this is just a way to simplify code, the user doesn't really have a
say in which operator is used, it is just hardcoded.
"""
function interactionop(s1::Index, s2::Index)
    if SiteType("S=1/2") ∈ sitetypes(s1) && SiteType("S=1/2") ∈ sitetypes(s2)
        t = -0.5 * (op("σ-", s1) * op("σ+", s2) + op("σ+", s1) * op("σ-", s2))
    elseif SiteType("Osc") ∈ sitetypes(s1) && SiteType("S=1/2") ∈ sitetypes(s2)
        t = op("X", s1) * op("σx", s2)
    elseif SiteType("S=1/2") ∈ sitetypes(s1) && SiteType("Osc") ∈ sitetypes(s2)
        t = op("σx", s1) * op("X", s2)
    elseif SiteType("Osc") ∈ sitetypes(s1) && SiteType("Osc") ∈ sitetypes(s2)
        t = (op("a+", s1) * op("a-", s2) + op("a-", s1) * op("a+", s2))
        #
    elseif SiteType("vecS=1/2") ∈ sitetypes(s1) && SiteType("vecS=1/2") ∈ sitetypes(s2)
        t =
            -0.5im * (
                op("σ-:Id", s1) * op("σ+:Id", s2) + op("σ+:Id", s1) * op("σ-:Id", s2) -
                op("Id:σ-", s1) * op("Id:σ+", s2) - op("Id:σ+", s1) * op("Id:σ-", s2)
            )
    elseif SiteType("vecOsc") ∈ sitetypes(s1) && SiteType("vecS=1/2") ∈ sitetypes(s2)
        t = im * (op("asum:Id", s1) * op("σx:Id", s2) - op("Id:asum", s1) * op("Id:σx", s2))
    elseif SiteType("vecS=1/2") ∈ sitetypes(s1) && SiteType("vecOsc") ∈ sitetypes(s2)
        t = im * (op("σx:Id", s1) * op("asum:Id", s2) - op("Id:σx", s1) * op("Id:asum", s2))
        #
    elseif SiteType("HvS=1/2") ∈ sitetypes(s1) && SiteType("HvS=1/2") ∈ sitetypes(s2)
        t =
            0.5im * (
                op("σ-⋅", s1) * op("σ+⋅", s2) + op("σ+⋅", s1) * op("σ-⋅", s2) -
                op("⋅σ-", s1) * op("⋅σ+", s2) - op("⋅σ+", s1) * op("⋅σ-", s2)
            )
    elseif SiteType("HvOsc") ∈ sitetypes(s1) && SiteType("HvS=1/2") ∈ sitetypes(s2)
        t = -im * (op("asum⋅", s1) * op("σx⋅", s2) - op("⋅asum", s1) * op("⋅σx", s2))
    elseif SiteType("HvS=1/2") ∈ sitetypes(s1) && SiteType("HvOsc") ∈ sitetypes(s2)
        t = -im * (op("σx⋅", s1) * op("asum⋅", s2) - op("⋅σx", s1) * op("⋅asum", s2))
    elseif SiteType("HvOsc") ∈ sitetypes(s1) && SiteType("HvOsc") ∈ sitetypes(s2)
        # Oscillator-oscillator interaction H = a₁†a₂ + a₂†a₁ that goes in the
        # commutator -i[H,ρ] as 
        #   -i (a₁†a₂ ρ + a₂†a₁ ρ - ρ a₁†a₂ - ρ a₂†a₁)
        t =
            -im * (
                op("a+⋅", s1) * op("a-⋅", s2) + op("a+⋅", s2) * op("a-⋅", s1) -
                op("⋅a+", s1) * op("⋅a-", s2) - op("⋅a+", s2) * op("⋅a-", s1)
            )
    else
        throw(DomainError((s1, s2), "Unrecognised SiteTypes."))
    end
    return t
end

"""
    evolve(initialstate, timesteplist, nskip, STorder, linksodd, linkseven, maxerr, maxrank; fout)

Evolves the initial state in time using the TEBD algorithm.
Return a list of the values of the functions passed within `fout` at each time
instant, together with a list of such instants.

# Arguments
- `initialstate` → MPS of the initial state
- `timesteplist` → list of the time steps at which the state is evolved
- `nskip`        → number of time steps to be skipped between each evaluation
  of the output functions
- `STorder`      → order of Suzuki-Trotter expansion
- `linksodd`     → two-site gates on odd sites
- `linkseven`    → two-site gates on even sites
- `maxerr`       → maximum truncation error of the MPS during evolution
- `maxrank`      → maximum allowed rank of the MPS during evolution
- `fout`         → list of functions which take an MPS as their only argument

The state is evolved through each of the instants specified in `timesteplist`,
but functions in `fout` are evaluated only after `nskip` steps each time. This
could save some computation time if a fine-grained output is not needed.

# Example
If `fout = [a1, a2]`, the function returns the tuple `(ts, b1, b2)` where
`b1` and `b2` are the result of applying `a1` and `a2` on the state for each
time instant in `ts`, that is `b1[j]` is equal to `a1` applied to the state at
`ts[j]`.
"""
function evolve(
    initialstate, timesteplist, nskip, STorder, linksodd, linkseven, maxerr, maxrank; fout
)
    @info "Computing Suzuki-Trotter decomposition of the evolution operator."
    τ = timesteplist[2] - timesteplist[1]
    if STorder == 1
        # Nothing to do if the order is 1.
        list = [linksodd(τ), linkseven(τ)]
        seq = repeat([1, 2], nskip)
    elseif STorder == 2
        # At 2nd order we can group the operators with t/2 (from two consecutive
        # steps) in order to save some time and numerical errors.
        # Instead of a series of 3*nskip operators between a measurement and the
        # next one, we have only 2*nskip+1 of them.
        list = [linksodd(0.5τ), linkseven(τ), linksodd(τ)]
        seq = [1; repeat([2, 3], nskip - 1); 2; 1]
    elseif STorder == 4
        c = (4 - 4^(1 / 3))^(-1)
        list = [
            linksodd(c * τ),
            linkseven(2c * τ),
            linksodd((0.5 - c) * τ),
            linkseven((1 - 4c) * τ),
        ]
        seq = repeat([1, 2, 3, 4, 3, 2, 1], nskip)
        # There could be still some margin for improvement: maybe group the first
        # and last `linksodd` when nskip > 1. It doesn't make that much of a
        # difference, though.
    else
        throw(
            DomainError(
                STorder,
                "Suzuki-Trotter expansion of order $STorder " *
                "is not supported. Only orders 1, 2 and 4 are " *
                "available.",
            ),
        )
    end

    # Applying `list` to a state means evolving it for a time equal to nskip*τ,
    # where τ is the integration time step.
    @info "Starting time evolution."
    state = initialstate
    returnvalues = [[f(state) for f in fout]]

    tout = timesteplist[1:nskip:end]
    progress = Progress(length(tout), 1, "Simulation in progress", 30)
    for _ in timesteplist[(1 + nskip):nskip:end]
        for n in seq
            state .= apply(list[n], state; cutoff=maxerr, maxdim=maxrank)
        end
        push!(returnvalues, [f(state) for f in fout])
        next!(progress)
    end

    # The list `returnvalues` contains N ≡ length(timesteplist[1:nskip:end])
    # lists as elements; in each one, the j-th element is the result of fout[j]
    # when applied to the state at a given time instant.
    # We want to reorganise the output, so that the function returns a list `tout`
    # with the instant at which the observables are measures, together with as
    # many lists as the functions in `fout`: the j-th list will be the result
    # of fout[j] applied to the state for each t in `tout`.
    outresults = [[returnvalues[k][j] for k in eachindex(tout)] for j in eachindex(fout)]
    return tout, outresults...
end

# Current operators
# -----------------

"""
    current(sites, leftsite::Int, rightsite::Int)

Return the current operator from site `leftsite` to `rightside` in `sites`.

# Usage
The sites (strictly) between the given sites are assigned the ITensor operator
with OpName "Z"; `sites(leftsite)` and `sites(rightsite)` are assigned "plus"
and "minus", and other sites are assigned "Id" operators. All these operators
thus must be correctly defined in order to use this function.
"""
function fermioncurrent(sites, leftsite::Int, rightsite::Int)
    # The function doesn't care about the SiteType on which it is applied.
    # This allows us to use it also on oscillators, when needed, for example
    # to measure the current between the first oscillator sites just beyond
    # the edge of the spin chain (it works!).
    # The oscillators need to be exactly at the `leftsite` and `rightsite`
    # positions, and not between them, otherwise the function doesn't work
    # since there is no "Z" operator defined for oscillators.
    n = rightsite - leftsite
    tags1 = repeat(["Id"], length(sites))
    tags2 = repeat(["Id"], length(sites))
    tags1[leftsite] = "plus"
    tags2[leftsite] = "minus"
    for l in (leftsite + 1):(rightsite - 1)
        tags1[l] = "Z"
        tags2[l] = "Z"
    end
    tags1[rightsite] = "minus"
    tags2[rightsite] = "plus"

    if (
        SiteType("S=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("Osc") ∈ sitetypes(first(sites))
    )
        op = im * (-1)^n * (MPO(sites, tags1) - MPO(sites, tags2))
    elseif (
        SiteType("vecS=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("HvS=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("vecOsc") ∈ sitetypes(first(sites)) ||
        SiteType("HvOsc") ∈ sitetypes(first(sites))
    )
        op =
            im *
            (-1)^n *
            (
                MPS(ComplexF64, sites, string.("vec", tags1)) -
                MPS(ComplexF64, sites, string.("vec", tags2))
            )
        # The MPS(sites, states) constructor assumes Float64 as the data type,
        # but here we need it to be complex.
    else
        throw(DomainError(s, "Unrecognised SiteType."))
    end
    return op
end

function forwardflux(sites, leftsite::Int, rightsite::Int)
    # Just one of the two terms in j_{k,k'}...
    n = rightsite - leftsite
    tags1 = repeat(["Id"], length(sites))
    tags1[leftsite] = "plus"
    for l in (leftsite + 1):(rightsite - 1)
        tags1[l] = "Z"
    end
    tags1[rightsite] = "minus"

    if (
        SiteType("S=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("Osc") ∈ sitetypes(first(sites))
    )
        op = im * (-1)^n * (MPO(sites, tags1) - MPO(sites, tags2))
    elseif (
        SiteType("vecS=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("HvS=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("vecOsc") ∈ sitetypes(first(sites)) ||
        SiteType("HvOsc") ∈ sitetypes(first(sites))
    )
        op = im * (-1)^n * MPS(ComplexF64, sites, string.("vec", tags1))
    else
        throw(DomainError(s, "Unrecognised SiteType."))
    end
    return op
end
function backwardflux(sites, leftsite::Int, rightsite::Int)
    # ...and the other one.
    n = rightsite - leftsite
    tags2 = repeat(["Id"], length(sites))
    tags2[leftsite] = "minus"
    for l in (leftsite + 1):(rightsite - 1)
        tags2[l] = "Z"
    end
    tags2[rightsite] = "plus"

    if (
        SiteType("S=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("Osc") ∈ sitetypes(first(sites))
    )
        op = im * (-1)^n * (MPO(sites, tags1) - MPO(sites, tags2))
    elseif (
        SiteType("vecS=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("HvS=1/2") ∈ sitetypes(first(sites)) ||
        SiteType("vecOsc") ∈ sitetypes(first(sites)) ||
        SiteType("HvOsc") ∈ sitetypes(first(sites))
    )
        op = im * (-1)^n * MPS(ComplexF64, sites, string.("vec", tags2))
    else
        throw(DomainError(s, "Unrecognised SiteType."))
    end
    return op
end

"""
    hamiltonian_xy(length::Integer, energy::Real, coupling::Real)

Return the Hamiltonian of an XY spin chain of given `length` and energies.
"""
function hamiltonian_xy(N, energy, coupling)
    op = OpSum()
    for j in 1:(N - 1)
        op += -0.5coupling, "S+", j, "S-", j + 1
        op += -0.5coupling, "S-", j, "S+", j + 1
    end
    for j in 1:N
        op += energy, "Sz", j
    end
    return op
end

"""
    lindbladian_xy(length::Integer, spin_energy::Real, spin_coupling::Real)

Return the Lindbladian operator ``-i[H,—]``, in vectorised form, associated to the
Hamiltonian ``H`` of an XY spin chain of given `length` and energies.
"""
function lindbladian_xy(chain_length::Integer, spin_energy::Real, spin_coupling::Real)
    N = chain_length
    ε = spin_energy
    λ = spin_coupling

    op = OpSum()
    for j in 1:(N - 1)
        op += +0.5im * λ, "σ-⋅", j, "σ+⋅", j + 1
        op += +0.5im * λ, "σ+⋅", j, "σ-⋅", j + 1
        op += -0.5im * λ, "⋅σ-", j, "⋅σ+", j + 1
        op += -0.5im * λ, "⋅σ+", j, "⋅σ-", j + 1
    end
    for j in 1:N
        op += -0.5im * ε, "σz⋅", j
        op += +0.5im * ε, "⋅σz", j
    end
    return op
end

avgn(ε, T) = T == 0 ? 0 : (ℯ^(ε / T) - 1)^(-1)

"""
    dissipator_symmetric(n, excitation_energy, spin_damping_coefficient, temperature)

Return the dissipator ``D(ρ) ∝ σˣ ρ (σˣ)† - ½ {(σˣ)† σˣ, ρ}`` on site `n`.
"""
function dissipator_symmetric(
    n::Integer, excitation_energy::Real, spin_damping_coefficient::Real, temperature::Real
)
    # TODO: move the ξ constant out of the function.
    ε = excitation_energy
    κ = spin_damping_coefficient
    T = temperature
    ξ = κ * (1 + 2avgn(ε, T))

    op = OpSum()
    op += +ξ, "σx⋅ * ⋅σx", n
    op += -ξ, "Id", n
    return op
end

"""
    dissipator_asymmetric(n, frequency, damping_coefficient, temperature)

Return the dissipator ``D(ρ) = γ⁻(σ⁻ρσ⁺ - ½ {σ⁺σ⁻, ρ}) + γ⁺(σ⁺ρσ⁻ - ½ {σ⁻σ⁺, ρ})`` on
site `n`.
"""
function dissipator_asymmetric(
    n::Integer, excitation_energy::Real, spin_damping_coefficient::Real, temperature::Real
)
    ω = excitation_energy
    γ = spin_damping_coefficient
    T = temperature

    γ⁺ = γ * avgn(ω, T)
    γ⁻ = γ * (1 + avgn(ω, T))

    op = OpSum()
    op += γ⁻, "σ-⋅ * ⋅σ+", n
    op += -0.5γ⁻, "σ+⋅ * σ-⋅", n
    op += -0.5γ⁻, "⋅σ- * ⋅σ+", n
    op += γ⁺, "σ+⋅ * ⋅σ-", n
    op += -0.5γ⁺, "σ-⋅ * σ+⋅", n
    op += -0.5γ⁺, "⋅σ+ * ⋅σ-", n
    return op
end

# Easy construction of the unitary part of the GKSL equation
# ----------------------------------------------------------
function makeopsumpairs(args...)
    isodd(length(args)) &&
        throw(error("Cannot make pairs out of an odd number of arguments."))
    return [(args[i], args[i + 1]) for i in eachindex(args)[1:2:(end - 1)]]
end

"""
    gkslcommutator(operators::Tuple{String,Int}...)

Given one or more tuples `(a1, n1), (a2, n2), …`, return an OpSum representing the operation
``x ↦ -i[A_1 A_2 …, x]`` where ``A_i`` is an operator which consists of `ai` at site `ni`
and the identity elsewhere. The string `a1` must be an existing ITensor OpName whose
variants `a1⋅` and `⋅a1` are defined (this function, however, doesn't perform any checks).

# Examples
```julia-repl
julia> gkslcommutator(("σ+", 1), ("σ-", 2))
sum(
  0.0 - 1.0im σ+⋅(1,) σ-⋅(2,)
  0.0 + 1.0im ⋅σ+(1,) ⋅σ-(2,)
)
```
"""
function gkslcommutator(operators::Tuple{String,Int}...)
    l = OpSum()
    opnames = first.(operators)
    sites = last.(operators)
    l += (-im, collect(Iterators.flatten(zip(opnames .* "⋅", sites)))...)
    l += (+im, collect(Iterators.flatten(zip("⋅" .* opnames, sites)))...)
    return l
end

"""
    gkslcommutator(args...)

Given a sequence `a1, n1, a2, n2, …` where each `ai` is a `String` and each `ni` is an
`Int`, return an OpSum representing the operation ``x ↦ -i[A_1 A_2 …, x]`` where ``A_i``
is an operator which consists of `ai` at site `ni` and the identity elsewhere. The string
`a1` must be an existing ITensor OpName whose variants `a1⋅` and `⋅a1` are defined
(this function, however, doesn't perform any checks).

# Examples
```julia-repl
julia> gkslcommutator("σ+", 1, "σ-", 2)
sum(
  0.0 - 1.0im σ+⋅(1,) σ-⋅(2,)
  0.0 + 1.0im ⋅σ+(1,) ⋅σ-(2,)
)
```
"""
gkslcommutator(args...) = gkslcommutator(makeopsumpairs(args...)...)

"""
    gkslcommutator_itensor(sites::Vector{<:Index}, operators::Tuple{String,Int}...)

Given a vector of ITensors site indices `sites` and a sequence `a1, n1, a2, n2, …` where
each `ai` is a `String` and each `ni` is an `Int`, return an ITensor with the site indices
`s[n1]`, `s[n1]'`, `s[n2]`, `s[n2]'` and so on which represents the operation
``x ↦ -i[A_1 A_2 …, x]`` where ``A_i`` is an operator consisting of `ai` at site
`ni` and the identity elsewhere. Each `ai` string must be an existing ITensors OpName whose
variants `a1⋅` and `⋅a1` are defined (this function, however, doesn't perform any checks).


# Examples

We start from a system of three 1/2-spins:
```julia
sites = siteinds("vS=1/2", 3)
```

Take an operator ``U`` which is ``Sˣ`` on the second site and the identity on the
others. The commutator ``ρ ↦ -i[U,ρ]`` is given by
```julia
gkslcommutator_itensor(sites, "Sx", 2)
```

An operator ``V`` which is ``Sy`` on the first site, the identity on the second one, and
``Sz`` on the third one. The commutator ``ρ ↦ -i[V,ρ]`` is given by
```julia
gkslcommutator_itensor(sites, "Sy", 1, "Sz", 3)
```
"""
function gkslcommutator_itensor(sites::Vector{<:Index}, operators::Tuple{String,Int}...)
    operator_names = first.(operators)
    site_numbers = last.(operators)
    if !allunique(sites)
        error("Some sites are repeated in the list. Please use unique site indices.")
        # This possibility is not allowed for now, since it would create issues in the
        # multiplication loop below. Basically, if two ITensors with the same indices
        # are multiplied together (with a simple `*`) then both indices get contracted
        # and we end up with a scalar. We should use `apply` instead, but it does not
        # work with OneITensor objects...
    end
    lmult = ITensors.OneITensor()
    rmult = ITensors.OneITensor()
    for (on, j) in zip(operator_names, site_numbers)
        lmult *= op("$on⋅", sites, j)
        rmult *= op("⋅$on", sites, j)
    end
    return -im * (lmult - rmult)
end

function gkslcommutator_itensor(sites::Vector{<:Index}, args...)
    return gkslcommutator_itensor(sites, makeopsumpairs(args...)...)
end
