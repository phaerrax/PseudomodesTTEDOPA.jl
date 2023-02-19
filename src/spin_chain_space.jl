# Additional ITensor states and operators (in addition to the ones already
# defined for S=1/2 sites):
# - null vector
ITensors.state(::StateName"0", ::SiteType"S=1/2") = [0; 0]
# - number operator, aka |↑⟩⟨↑|
ITensors.op(::OpName"N", st::SiteType"S=1/2") = op(OpName"ProjUp"(), st)
# - identity matrices
ITensors.op(::OpName"Id", ::SiteType"S=1/2") = I₂
# - null matrix
ITensors.op(::OpName"0", ::SiteType"S=1/2") = zeros(2, 2)
# - ladder operators (aliases)
ITensors.op(::OpName"σ+", st::SiteType"S=1/2") = op(OpName("S+"), st)
ITensors.op(::OpName"σ-", st::SiteType"S=1/2") = op(OpName("S-"), st)
ITensors.op(::OpName"plus", st::SiteType"S=1/2") = op(OpName("S+"), st)
ITensors.op(::OpName"minus", st::SiteType"S=1/2") = op(OpName("S-"), st)

# Spin current operators 
# ======================

function J⁺tag(::SiteType"S=1/2", left_site::Int, i::Int)
  # TODO: still useful?
  if i == left_site
    str = "σx"
  elseif i == left_site + 1
    str = "σy"
  else
    str = "Id"
  end
  return str
end
function J⁺tag(::SiteType"vecS=1/2", left_site::Int, i::Int)
  if i == left_site
    str = "vecσx"
  elseif i == left_site + 1
    str = "vecσy"
  else
    str = "vecId"
  end
  return str
end
function J⁺tag(::SiteType"HvS=1/2", left_site::Int, i::Int)
  return J⁺tag(SiteType("vecS=1/2"), left_site, i)
end

function J⁻tag(::SiteType"S=1/2", left_site::Int, i::Int)
  # Just as `J⁺tag`, but for σʸ⊗σˣ
  if i == left_site
    str = "σy"
  elseif i == left_site + 1
    str = "σx"
  else
    str = "Id"
  end
  return str
end
function J⁻tag(::SiteType"vecS=1/2", left_site::Int, i::Int)
  if i == left_site
    str = "vecσy"
  elseif i == left_site + 1
    str = "vecσx"
  else
    str = "vecId"
  end
  return str
end
function J⁻tag(::SiteType"HvS=1/2", left_site::Int, i::Int)
  return J⁻tag(SiteType("vecS=1/2"), left_site, i)
end

function spin_current_op_list(sites::Vector{Index{Int64}})
  N = length(sites)
  # Check if all sites are spin-½ sites.
  if all(x -> SiteType("S=1/2") in x, sitetypes.(sites))
    st = SiteType("S=1/2")
    MPtype = MPO
  elseif all(x -> SiteType("vecS=1/2") in x, sitetypes.(sites))
    st = SiteType("vecS=1/2")
    MPtype = MPS
  elseif all(x -> SiteType("HvS=1/2") in x, sitetypes.(sites))
    st = SiteType("HvS=1/2")
    MPtype = MPS
  else
    throw(ArgumentError("spin_current_op_list works with SiteTypes "*
                        "\"S=1/2\", \"vecS=1/2\" or \"HvS=1/2\"."))
  end
  #
  J⁺ = [MPtype(sites, [J⁺tag(st, k, i) for i = 1:N]) for k = 1:N-1]
  J⁻ = [MPtype(sites, [J⁻tag(st, k, i) for i = 1:N]) for k = 1:N-1]
  return -0.5 .* (J⁺ .- J⁻)
end


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
  return unique(permutations([repeat(["Up"], level);
                              repeat(["Dn"], n - level)]))
end

"""
    level_subspace_proj(sites::Vector{Index{Int64}}, l::Int)

Return the projector on the subspace with `l` "Up" spins.
"""
function level_subspace_proj(sites::Vector{Index{Int64}}, l::Int)
  N = length(sites)
  # Check if all sites are spin-½ sites.
  if all(x -> SiteType("S=1/2") in x, sitetypes.(sites))
    projs = [projector(MPS(sites, names); normalize=false)
             for names in chain_basis_states(N, l)]
  elseif all(x -> SiteType("vecS=1/2") in x, sitetypes.(sites)) ||
         all(x -> SiteType("HvS=1/2") in x, sitetypes.(sites))
    projs = [MPS(sites, names)
             for names in chain_basis_states(N, l)]
  else
    throw(ArgumentError("level_subspace_proj works with SiteTypes "*
                        "\"S=1/2\", \"vecS=1/2\" or \"HvS=1/2\"."))
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
    states = [i == k ? "Up" : "Dn" for i ∈ 1:N]
  else
    throw(DomainError(k,
                      "Trying to build a state with an excitation localised "*
                      "at site $k, which does not belong to the chain: please "*
                      "insert a value between 1 and $N."))
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
  if j ∉ 1:N
    throw(DomainError(j,
                      "Trying to build a chain eigenstate with invalid index "*
                      "$j: please insert a value between 1 and $N."))
  end
  states = [2/(N+1) * sin(j*k*π / (N+1))^2 * single_ex_state(sites, k)
            for k ∈ 1:N]
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
    throw(DomainError(state,
                      "Unrecognised state: please choose from \"empty\", "*
                      "\"1locN\" or \"1eigN\"."))
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
    v = 1/sqrt(2) * (ITensors.state(site, "Up") + ITensors.state(site, "Dn"))
  else
    throw(DomainError(state,
                      "Unrecognised state: please choose from \"empty\",
                      \"up\", \"down\" or \"x+\"."))
  end
  return MPS([v])
end
