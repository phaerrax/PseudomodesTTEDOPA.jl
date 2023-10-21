"""
    sitenumber(i::Index)

Return the site number of the given Index, i.e. the number in the "n=N" tag.
"""
function sitenumber(i::Index)
    st = split(replace(string.(tags.(i)), "\"" => ""), ",")
    sitetag = first(filter!(s -> occursin("n=", s), st))
    # TODO what if sitetag is empty?
    siten = replace(sitetag, "n=" => "")
    return parse(Int, siten)
end

"""
    jwstring(; start, stop, op::AbstractString="F")

Return the vector ``(op, start + 1, op, start + 2, ..., op, stop - 1)``.
"""
function jwstring(; start, stop, op::AbstractString="F")
    return collect(Iterators.flatten([(op, start + k) for k in 1:(stop - start - 1)]))
end

"""
    groupresults(obs::AbstractObserver, name::String)

If the observer returns an array of values, rearrange the results in a big matrix.
"""
function groupresults(obs::AbstractObserver, name::String)
    return mapreduce(permutedims, vcat, obs[!, name])
end

"""
    disablegrifqtech()

Disable the graphical output if the script is running on Qtech.
"""
function disablegrifqtech()
    # The script crashes if it is executed on Qtech, unless we disable
    # the graphical output.
    if gethostname() == "qtech.fisica.unimi.it" || gethostname() == "qtech2.fisica.unimi.it"
        ENV["GKSwstype"] = "100"
        @info "Running on remote server. Disabling graphical output."
    else
        delete!(ENV, "GKSwstype")
        # If the key "GKSwstype" doesn't exist then nothing happens.
    end
end

# Vectorisation utilities
# =======================
# In order to use tr(x'*y) as a tool to extract coefficient the basis must
# of course be orthonormal wrt this inner product.
# The canonical basis or the Gell-Mann one are okay.

"""
    vec(A::Matrix, basis::Vector)

Compute the vector of coefficients of the matrix `A` wrt the basis `basis`.
"""
function vec(A::Matrix, basis::Vector)
    return [tr(b' * A) for b in basis]
end
"""
    vec(L::Function, basis::Vector)

Compute the matrix of coefficients of the linear map `L` wrt the basis `basis`.
The linearity of the map is not checked, so using this function with non-linear
functions leads to undefined results.
"""
function vec(L::Function, basis::Vector)
    return [tr(bi' * L(bj)) for (bi, bj) in Base.product(basis, basis)]
end

"""
    partialtrace(sites::Vector{Index{Int64}}, v::MPS, j::Int)

Compute the partial trace, on the `j`-th site of `sites`, of the matrix
represented, as vectorised, by the MPS `v`.
The result is a `Vector` containing the coordinates of the partial trace
(i.e. the result is still in vectorised form).
"""
function partialtrace(sites::Vector{Index{Int64}}, v::MPS, j::Int)
    # Custom contraction sequence: we start from v[end] and we contract it with
    # a vecId state; we contract the result with v[end-1] and so on, until we
    # get to v[j].  # Then we do the same starting from v[1].
    x = ITensor(1.0)
    for i in length(sites):-1:(j + 1)
        x = v[i] * state(sites[i], "vecId") * x
    end
    y = ITensor(1.0)
    for i in 1:(j - 1)
        y = v[i] * state(sites[i], "vecId") * y
    end
    z = y * v[j] * x
    # Now `vector(z)` is the coordinate vector of the partial trace.
    return vector(z)
end

# Generalised Gell-Mann matrices
# ==============================
"""
    gellmannmatrix(j, k, dim)

Return the `(j,k)` generalised Gell-Mann matrix of dimension `dim`, normalised
wrt the Hilbert-Schmidt inner product ``(A,B) = tr(A†B)``.
The matrices are indexed as follows:

    * if ``j > k`` the matrix is symmetric and traceless;
    * if ``j < k`` the matrix is antisymmetric;
    * if ``j = k`` the matrix is diagonal.

In particular, ``j = k = dim`` gives a matrix proportional to the identity.
The two indices `j` and `k` determine the non-zero coefficients of the matrix.
The whole set of (different) Gell-Mann matrices that can be generated with this
function is a basis of ``Mat(ℂᵈⁱᵐ)``.
"""
function gellmannmatrix(j, k, dim)
    if j > dim || k > dim || j < 0 || k < 0
        throw(DomainError)
    end
    m = zeros(ComplexF64, dim, dim)
    if j > k
        m[j, k] = 1 / sqrt(2)
        m[k, j] = 1 / sqrt(2)
    elseif k > j
        m[j, k] = -im / sqrt(2)
        m[k, j] = im / sqrt(2)
    elseif j == k && j < dim
        for i in 1:j
            m[i, i] = 1
        end
        m[j + 1, j + 1] = -j
        m .*= sqrt(1 / (j * (j + 1)))
    else
        for i in 1:dim
            m[i, i] = 1 / sqrt(dim)
        end
    end
    return m
end

"""
    gellmannbasis(dim)

Return a list containing the "Hermitian basis" of ``Mat(ℂᵈⁱᵐ)``, i.e. composed
of the ``dim²`` generalised Gell-Mann matrices.
"""
function gellmannbasis(dim)
    return [gellmannmatrix(j, k, dim) for (j, k) in [Base.product(1:dim, 1:dim)...]]
    # We need to splat the result from `product` so that the result is a list
    # of matrices (a Vector) and not a Matrix.
end

"""
    canonicalmatrix(i, j, dim)

Return the (`i`,`j`) element of the canonical basis of ``Mat(ℂᵈⁱᵐ)``, i.e. a
`dim`×`dim` matrix whose element on the `i`-th row and `j`-th column is ``1``,
and zero elsewhere.
"""
function canonicalmatrix(i, j, dim)
    m = zeros(ComplexF64, dim, dim)
    m[i, j] = 1
    return m
end

"""
    canonicalbasis(dim)

Return a list of the matrices in the canonical basis of ``Mat(ℂᵈⁱᵐ)``. 
The list is ordered corresponding to column-based vectorisation, i.e.

    canonicalbasis(dim)[j] = canonicalmatrix((j-1)%dim + 1, (j-1)÷dim + 1, dim)

with ``j ∈ {1,…,dim²}``. With this ordering,
``vec(A)ⱼ = tr(canonicalbasis(dim)[j]' * A)``.
"""
function canonicalbasis(dim)
    return [canonicalmatrix(i, j, dim) for (i, j) in [Base.product(1:dim, 1:dim)...]]
end

# Spin-half space utilities
# =========================

# Chain eigenstate basis
# ----------------------
# How to measure how much each eigenspace of the number operator (of the
# whole spin chain) "contributes" to a given state ρ?
# We build a projector operator associated to each eigenspace.
# The projector on the m-th level eigenspace can be made by simply adding the
# orthogonal projections on each element of the canonical basis which has m
# "Up" spins and N-m "Down" spins, N being the length of the chain.
# This may result in a very big MPO. If these are needed for more than one
# simulation, calculating them once and for all before the simulations start
# may help to cut down the computation time.
"""
    chain_basis_states(n::Int, level::Int)

Return a list of strings which can be used to build the MPS of all states in
the ``ℂ²ⁿ`` basis that contain `level` "Up" spins.
"""
function chain_basis_states(n::Int, level::Int)
    return unique(permutations([
        repeat(["Up"], level)
        repeat(["Dn"], n - level)
    ]))
end

"""
    level_subspace_proj(sites::Vector{Index{Int64}}, l::Int)

Return the projector on the subspace with `l` "Up" spins.
"""
function level_subspace_proj(sites::Vector{Index{Int64}}, l::Int)
    N = length(sites)
    # Check if all sites are spin-½ sites.
    if all(x -> SiteType("S=1/2") in x, sitetypes.(sites))
        projs = [
            projector(MPS(sites, names); normalize=false) for
            names in chain_basis_states(N, l)
        ]
    elseif all(x -> SiteType("vecS=1/2") in x, sitetypes.(sites)) ||
        all(x -> SiteType("HvS=1/2") in x, sitetypes.(sites))
        projs = [MPS(sites, names) for names in chain_basis_states(N, l)]
    else
        throw(
            ArgumentError(
                "level_subspace_proj works with SiteTypes " *
                "\"S=1/2\", \"vecS=1/2\" or \"HvS=1/2\".",
            ),
        )
    end
    # Somehow return sum(projs) doesn't work… we have to sum manually.
    P = projs[1]
    for p in projs[2:end]
        P = +(P, p; cutoff=1e-10)
    end
    return P
end

# First-level chain eigenstates
# -----------------------------
# It is useful to have at hand the eigenstates of the free chain Hamiltonian
# within the first-level eigenspace of the number operator ([H,N] = 0).
# States |sₖ⟩ with an up spin on site k and a down spin on the others are
# in fact not eigenstates, but they still form a basis for the eigenspace;
# wrt the {sₖ}ₖ (k ∈ {1,…,N}) basis the free chain Hamiltonian is written as
#   ⎛ ε λ 0 0 … 0 ⎞
#   ⎜ λ ε λ 0 … 0 ⎟
#   ⎜ 0 λ ε λ … 0 ⎟
#   ⎜ 0 0 λ ε … 0 ⎟
#   ⎜ ⋮ ⋮ ⋮ ⋮ ⋱ ⋮ ⎟
#   ⎝ 0 0 0 0 … ε ⎠
# whose eigenstates are |vⱼ⟩= ∑ₖ sin(kjπ /(N+1)) |sₖ⟩, con j ∈ {1,…,N}.
# Note that they are not normalised: ‖|vⱼ⟩‖² = (N+1)/2.
"""
    single_ex_state(sites::Vector{Index{Int64}}, k::Int)

Return the MPS of a state with a single excitation on the `k`-th site.
"""
function single_ex_state(sites::Vector{Index{Int64}}, k::Int)
    N = length(sites)
    if k ∈ 1:N
        states = [i == k ? "Up" : "Dn" for i in 1:N]
    else
        throw(
            DomainError(
                k,
                "Trying to build a state with an excitation localised " *
                "at site $k, which does not belong to the chain: please " *
                "insert a value between 1 and $N.",
            ),
        )
    end
    return MPS(sites, states)
end

"""
    chain_L1_state(sites::Vector{Index{Int64}}, j::Int)

Return the `j`-th eigenstate of the chain Hamiltonian in the single-excitation
subspace.
"""
function chain_L1_state(sites::Vector{Index{Int64}}, j::Int)
    # FIXME: this seems just wrong. Why not computing directly vec(|vⱼ⟩ ⊗ ⟨vⱼ|), 
    # without going through the linear combination?
    #
    # Careful with the coefficients: this isn't |vⱼ⟩ but vec(|vⱼ⟩ ⊗ ⟨vⱼ|),
    # so if we want to build it as a linear combination of the |sₖ⟩'s we
    # need to square the coefficients.
    # Note that vectorised projectors satisfy
    #     ⟨vec(a⊗aᵀ), vec(b⊗bᵀ)⟩ = tr((a⊗aᵀ)ᵀ b⊗bᵀ) = (aᵀb)².
    # We don't expect the norm of this MPS to be 1. What has to be 1, is the
    # sum of the inner products of vec(|sₖ⟩ ⊗ ⟨sₖ|) and this vector, for
    # k ∈ {1,…,N}. We get ⟨Pₙ|vⱼ⟩ = 1 only if n=j, otherwise it is 0.
    N = length(sites)
    if !(j in 1:N)
        throw(
            DomainError(
                j,
                "Trying to build a chain eigenstate with invalid index " *
                "$j: please insert a value between 1 and $N.",
            ),
        )
    end
    states = [
        2 / (N + 1) * sin(j * k * π / (N + 1))^2 * single_ex_state(sites, k) for k in 1:N
    ]
    return sum(states)
end

# Choice of the spin chain's initial state
# ----------------------------------------
"""
    parse_init_state(sites::Vector{Index{Int64}}, state::String)

Return an MPS representing a particular state of the spin chain, given
by the string `state`:

  * "empty" -> empty state (aka the ground state of the chain Hamiltonian)
  * "1locM" -> state with a single excitation at site `M` (``M ∈ {1,…,N}``)
  * "1eigM" -> single-excitation eigenstate of the chain Hamiltonian with `M` nodes (``M ∈ {0,…,N-1}``)

The string is case-insensitive. The length `N` of the chain is computed
from `sites`. 
"""
function parse_init_state(sites::Vector{Index{Int64}}, state::String)
    state = lowercase(state)
    if state == "empty"
        v = MPS(sites, "Dn")
    elseif occursin(r"^1loc", state)
        j = parse(Int, replace(state, "1loc" => ""))
        v = single_ex_state(sites, j)
    elseif occursin(r"^1eig", state)
        j = parse(Int, replace(state, "1eig" => ""))
        # The j-th eigenstate has j-1 nodes
        v = chain_L1_state(sites, j + 1)
    else
        throw(
            DomainError(
                state,
                "Unrecognised state: please choose from \"empty\", " *
                "\"1locN\" or \"1eigN\".",
            ),
        )
    end
    return v
end

"""
    parse_spin_state(site::Index{Int64}, state::String)

Return the MPS of a single spin site representing a particular state, given
by the string `state`:

- "empty", "dn", "down" → spin-down state
- "up"                  → spin-up state
- "x+"                  → ``1/√2 ( |+⟩ + |-⟩ )`` state
"""
function parse_spin_state(site::Index{Int64}, state::String)
    state = lowercase(state)
    if state == "empty" || state == "dn" || state == "down"
        v = ITensors.state(site, "Dn")
    elseif state == "up"
        v = ITensors.state(site, "Up")
    elseif state == "x+"
        v = 1 / sqrt(2) * (ITensors.state(site, "Up") + ITensors.state(site, "Dn"))
    else
        throw(DomainError(
            state,
            "Unrecognised state: please choose from \"empty\",
            \"up\", \"down\" or \"x+\".",
        ))
    end
    return MPS([v])
end

# Oscillator space utilities
# ==========================

"""
	oscdimensions(N, basedim, decay)

Compute a decreasing sequence of length `N` where the first two elements are
equal to `basedim` and the following ones are given by
`floor(2 + basedim * ℯ^(-decay * n))`.

Useful to determine the dimensions of oscillator sites in a TEDOPA chain.
"""
function oscdimensions(length, basedim, decay)
    f(j) = 2 + basedim * ℯ^(-decay * j)
    return [basedim; basedim; (Int ∘ floor ∘ f).(3:length)]
end

"""
    parse_init_state_osc(site::Index{Int64},
                         statename::String;
                         <keyword arguments>)

Return an MPS representing a particular state of a harmonic oscillator, given
by the string `statename`:

- "thermal" → thermal equilibrium state
- "fockN"   → `N`-th eigenstate of the number operator (element of Fock basis)
- "empty"   → alias for "fock0"

The string is case-insensitive. Other parameters required to build the state
(e.g. frequency, temperature) may be supplied as keyword arguments.
"""
function parse_init_state_osc(site::Index{Int64}, statename::String; kwargs...)
    # TODO: maybe remove "init" from title? It is a generic state, after all.
    statename = lowercase(statename)
    if statename == "thermal"
        s = state(site, "ThermEq"; kwargs...)
    elseif occursin(r"^fock", statename)
        j = parse(Int, replace(statename, "fock" => ""))
        s = state(site, "$j")
    elseif statename == "empty"
        s = state(site, "0")
    else
        throw(
            DomainError(
                statename,
                "Unrecognised state name; please choose from " *
                "\"empty\", \"fockN\" or \"thermal\".",
            ),
        )
    end
    return MPS([s])
end

# (Von Neumann) Entropy
# =====================
"""
    vonneumannentropy(ψ::MPS, sites::Vector{Index{Int64}}, n::Int)

Compute the entanglement entropy of the biparition ``(1,…,n)|(n+1,…,N)`` of
the system in state described by the MPS `ψ` (defined on the sites `sites`),
using its Schmidt decomposition.
"""
function vonneumannentropy(ψ₀::MPS, sites::Vector{Index{Int64}}, n::Int)
    ψ = orthogonalize(ψ₀, n)
    # Decompose ψ[n] in singular values, treating the Link between sites n-1 and n
    # and the physical index as "row index"; the remaining index, the Link
    # between n and n+1, is the "column index".
    _, S, _ = svd(ψ[n], (linkind(ψ, n - 1), sites[n]))
    # Compute the square of the singular values (aka the Schmidt coefficients
    # of the bipartition), and from them the entropy.
    sqdiagS = [S[j, j]^2 for j in ITensors.dim(S, 1)]
    return -sum(p -> p * log(p), sqdiagS; init=0.0)
end

# Chop (from Mathematica)
# =======================
# Imitation of Mathematica's "Chop" function

"""
    chop(x::Real; tolerance=1e-10)

Truncates `x` to zero if it is less than `tolerance`.
"""
function chop(x::Real; tolerance=1e-10)
    return abs(x) > tolerance ? x : zero(x)
end

"""
    chop(x::Complex; tolerance=1e-10)

Truncates the real and/or the imaginary part of `x` to zero if they are less
than `tolerance`.
"""
function chop(x::Complex; tolerance=1e-10)
    return Complex(chop(real(x)), chop(imag(x)))
end

# Manipulating ITensor objects
# ============================

"""
    sitetypes(s::Index)

Return the ITensor tags of `s` as SiteTypes. 

This function is already defined in the ITensor library, but it is not publicly
accessible.
"""
function sitetypes(s::Index)
    ts = tags(s)
    return SiteType[SiteType(ts.data[n]) for n in 1:length(ts)]
end

# Reading the parameters of the simulations
# =========================================

"""
    isjson(filename::String)

Check if `filename` ends in ".json".

By design, filenames consisting of only ".json" return `false`.
"""
function isjson(filename::String)
    return length(filename) > 5 && filename[(end - 4):end] == ".json"
    # We ignore a file whose name is only ".json".
end

"""
    load_parameters(file_list)

Load the JSON files contained in `file_list` into dictionaries, returning a
list of dictionaries, one for each file.

If `file_list` is a filename, then the list will contain just one dictionary;
if `file_list` is a directory, every JSON file within it is loaded and a
dictionary is created for each one of them.
"""
function load_parameters(file_list)
    if isempty(file_list)
        throw(ErrorException("No parameter file provided."))
    end
    first_arg = file_list[1]
    prev_dir = pwd()
    if isdir(first_arg)
        # If the first argument is a directory, we read all the JSON files within
        # and load them as parameter files; in the end all output will be saved
        # in that same directory. We ignore the remaining elements in ARGS.
        cd(first_arg)
        files = filter(isjson, readdir())
        @info "$first_arg is a directory. Ignoring other command line arguments."
    else
        # Otherwise, all command line arguments are treated as parameter files.
        # Output files will be saved in the pwd() from which the script was
        # launched.
        files = file_list
    end
    # Load parameters into dictionaries, one for each file.
    parameter_lists = []
    for f in files
        open(f) do input
            s = read(input, String)
            # Add the filename too to the dictionary.
            push!(parameter_lists, merge(Dict("filename" => f), JSON.parse(s)))
        end
    end
    cd(prev_dir)
    return parameter_lists
end

function allequal(a)
    return all(x -> x == first(a), a)
end

# Defining the time interval of the simulation
# ============================================

"""
    construct_step_list(parameters)

Return a list of time instants at which the time evolution will be evaluated.

The values run from zero up to ``parameters["simulation_end_time"]``, with a
step size equal to ``parameters["simulation_time_step"]``.
"""
function construct_step_list(parameters)
    τ = parameters["simulation_time_step"]
    end_time = parameters["simulation_end_time"]
    return collect(range(0, end_time; step=τ))
end

# Computing expectation values of observables
# ===========================================

# Function that calculates the eigenvalues of the number operator, given a set
# of projectors on the eigenspaces, and also return their sum (which should
# be 1).
# TODO: these two functions are probably outdated. Anyway they are not
# specific to the occupation levels, so they should have a more general
# description.
function levels(projs::Vector{MPO}, state::MPS)
    lev = [real(inner(state, p * state)) for p in projs]
    return [lev; sum(lev)]
end
function levels(projs::Vector{MPS}, state::MPS)
    lev = [real(inner(p, state)) for p in projs]
    return [lev; sum(lev)]
end

# MPS and MPO utilities
# =====================

"""
    vectrace(vecρ::MPS, s::Vector{Index{Int64}})

Return the trace of the density matrix ρ which is vectorized and encoded in the
given MPS `vecρ` on sites `s`.
"""
function vectrace(vecρ::MPS, s::Vector{Index{Int64}})
    return dot(MPS("vecId", s), vecρ)
end

"""
    linkdims(m::Union{MPS, MPO})

Return a list of the bond dimensions of `m`.
"""
function linkdims(m::Union{MPS,MPO})
    return [ITensors.dim(linkind(m, j)) for j in 1:(length(m) - 1)]
end

"""
    chain(left::MPS, right::MPS)

Concatenate `left` and `right`, returning `left` ``⊗`` `right`.
"""
function chain(left::MPS, right::MPS)
    # This function is like mpnum's `chain`: it takes two MPSs ans concatenates
    # them. The code is "inspired" from `ITensors.MPS` in mps.jl:308.

    midN = length(left) # The site with the missing link between the two MPSs.
    # First of all we shift the Link tags of the given MPSs, so that the final
    # enumeration of the tags is correct.
    # Note that in each MPS the numbers in the Link tags do not follow the
    # numbering of the Sites on which it is based, they always start from 1.
    for j in eachindex(left)
        replacetags!(left[j], "l=$j", "l=$j"; tags="Link")
        replacetags!(left[j], "l=$(j-1)", "l=$(j-1)"; tags="Link")
    end
    for j in eachindex(right)
        replacetags!(right[j], "l=$j", "l=$(midN+j)"; tags="Link")
        replacetags!(right[j], "l=$(j-1)", "l=$(midN+j-1)"; tags="Link")
    end
    # "Shallow" concatenation of the MPSs (there's still a missing link).
    M = MPS([left..., right...])
    # We create a "trivial" Index of dimension 1 and add it to the two sites
    # which are not yet connected.
    # The Index has dimension 1 because this is a tensor product between the
    # states so there's no correlation between them.
    missing_link = Index(1; tags="Link,l=$midN")
    M[midN] = M[midN] * state(missing_link, 1)
    M[midN + 1] = state(dag(missing_link), 1) * M[midN + 1]

    return M
end

"""
    chain(left::MPO, right::MPO)

Concatenate `left` and `right`, returning `left` ``⊗`` `right`.
"""
function chain(left::MPO, right::MPO)
    # Like the previous `chain`, but for MPOs.
    midN = length(left)# The site with the missing link between the two MPOs.
    for j in eachindex(right)
        replacetags!(right[j], "l=$j", "l=$(midN+j)"; tags="Link")
        replacetags!(right[j], "l=$(j-1)", "l=$(midN+j-1)"; tags="Link")
    end
    M = MPO([left..., right...])
    missing_link = Index(1; tags="Link,l=$midN")
    M[midN] = M[midN] * state(missing_link, 1)
    M[midN + 1] = M[midN + 1] * state(dag(missing_link), 1)
    # The order of the Indexes in M[midN] and M[midN+1] ns not what we would get
    # if we built an MPO on the whole system as usual; namely, the two "Link"
    # Indexes are swapped.
    # This however should not matter since ITensor does not care about the order.
    return M
end

# Varargs versions

"""
    chain(a::MPS, b...)

Concatenate the given MPSs into a longer MPS, returning their tensor product.
"""
chain(a::MPS, b...) = chain(a, chain(b...))
"""
    chain(a::MPO, b...)

Concatenate the given MPOs into a longer MPO, returning their tensor product.
"""
chain(a::MPO, b...) = chain(a, chain(b...))

"""
    embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPO)

Embed `slice`, defined on a subset `range` of `sites`, into a MPO which covers
the whole `sites`.

The MPO is extended by filling the empty spots with an "Id" operator, therefore
an operator with OpName "Id" is required to be defined for the SiteTypes of the
remaining sites.

# Arguments
- `sites::Array{Index{Int64}}`: the sites of the whole system.
- `range::UnitRange{Int}`: the range spanned by `slice`.
- `slice::MPO`: the MPO to be extended.
"""
function embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPO)
    # TODO: compute automatically on which sites the MPO is defined, without
    # having to supply the range explicitly as an argument.
    if length(slice) != length(range)
        throw(DimensionMismatch("slice and range must have the same size."))
    end
    if !issubset(range, eachindex(sites))
        throw(BoundsError(range, sites))
    end

    if range[begin] == 1 && range[end] == length(sites)
        mpo = slice
    elseif range[begin] == 1
        mpo = chain(slice, MPO(sites[(range[end] + 1):end], "Id"))
    elseif range[end] == length(sites)
        mpo = chain(MPO(sites[1:(range[begin] - 1)], "Id"), slice)
    else
        mpo = chain(
            MPO(sites[1:(range[begin] - 1)], "Id"),
            slice,
            MPO(sites[(range[end] + 1):end], "Id"),
        )
    end
    return mpo
end

"""
    embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPS)

Embed `slice`, defined on a subset `range` of `sites`, into a MPS which covers
the whole `sites` (to be interpreted as a vectorised operator).

The MPS is extended by filling the empty spots with a "vecId" operator,
therefore an operator with OpName "vecId" is required to be defined for the
SiteTypes of the remaining sites.

# Arguments
- `sites::Array{Index{Int64}}`: the sites of the whole system.
- `range::UnitRange{Int}`: the range spanned by `slice`.
- `slice::MPS`: the MPS to be extended.
"""
function embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPS)
    if length(slice) != length(range)
        throw(DimensionMismatch("slice e range must have the same size."))
    end
    if !issubset(range, eachindex(sites))
        throw(BoundsError(range, sites))
    end

    if range[begin] == 1 && range[end] == length(sites)
        mpo = slice
    elseif range[begin] == 1
        mpo = chain(slice, MPS(sites[(range[end] + 1):end], "vecId"))
    elseif range[end] == length(sites)
        mpo = chain(MPS(sites[1:(range[begin] - 1)], "vecId"), slice)
    else
        mpo = chain(
            MPS(sites[1:(range[begin] - 1)], "vecId"),
            slice,
            MPS(sites[(range[end] + 1):end], "vecId"),
        )
    end
    return mpo
end

# Other utilities
# ===============

"""
    filenamett(::Dict)

Extracts the filename from the parameter list, supplied as a Dict, and
formats it in a saw that's safe to use as text in TikZ nodes, plots, etc.,
in typewriter font.
"""
function filenamett(d::Dict)
    # Get basename and remove the .json extension.
    filename = replace(basename(d["filename"]), ".json" => "")
    # Sanitise string, escaping underscores.
    filename = replace(filename, "_" => "\\_")
    # Add \texttt command
    filename = raw"\texttt{" * filename * "}"
    return filename
end

"""
    consecutivepairs(v::AbstractVector)

Return a list of Strings "(a,b)" formed by all adjacent items in `v`.
"""
function consecutivepairs(v::AbstractVector)
    return string.("(", v[1:(end - 1)], ",", v[2:end], ")")
end

"""
    try_op(on::OpName, st::SiteType; kwargs...)

Return the matrix of an ITensor operator, if it exists, trying first the
```julia
op(::OpName, ::SiteType; kwargs...)
```
syntax, and then
```julia
op!(::ITensor, ::OpName, ::SiteType, ::Index...; kwargs...)
```
if the former returns nothing.
"""
function try_op(on::OpName, st::SiteType; kwargs...)
    stname = String(ITensors.tag(st))
    opstring = String(ITensors.name(on))
    # Try calling a function of the form:
    #    op(::OpName, ::SiteType; kwargs...)
    # which returns a Julia matrix
    mat = ITensors.op(on, st; kwargs...)
    if isnothing(mat)
        # Otherwise try calling a function of the form
        #    op!(::ITensor, ::OpName, ::SiteType, ::Index...; kwargs...)
        dummy = siteind(stname)
        Op = ITensor(prime(dummy), ITensors.dag(dummy))
        r = ITensors.op!(Op, on, st, dummy; kwargs...)
        if isnothing(r)
            throw(
                  ArgumentError(
                                "Overload of \"op\" or \"op!\" functions not found for " *
                                "operator name \"$opstring\" and Index tag $(tags(dummy)).",
                               ),
                 )
        end
        mat = matrix(Op)
    end
    if isnothing(mat)
        throw(
              ArgumentError(
                            "Overload of \"op\" or \"op!\" functions not found for operator " *
                            "name \"$opstring\" and Index tag \"$stname\".",
                           ),
             )
    end
    return mat
end
