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

# Space of harmonic oscillators
# =============================

ITensors.alias(::SiteType"Osc") = SiteType"Qudit"()

"""
    ITensors.space(st::SiteType"Osc";
                   dim = 2,
                   conserve_qns = false,
                   conserve_number = false,
                   qnname_number = "Number")

Create the Hilbert space for a site of type "Osc".

Optionally specify the conserved symmetries and their quantum number labels.
"""
ITensors.space(st::SiteType"Osc"; kwargs...) = space(ITensors.alias(st); kwargs...)

# Forward "Osc" definitions to "Qudit" states.
# This way each state, operator and so on defined for the "Qudit" SiteType is already
# available for the "Osc" type.
ITensors.val(vn::ValName; st::SiteType"Osc") = val(vn, ITensors.alias(st))

function ITensors.state(sn::StateName, st::SiteType"Osc", s::Index; kwargs...)
    return state(sn, ITensors.alias(st), s; kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"Osc", dims::Int...; kwargs...)
    return op(on, ITensors.alias(st), dims...; kwargs...)
end

# ITensors.op functions all require an additional parameter, the dimension: for this reason,
# we need to compute ITensors.dim.(rs), i.e. the dimensions of the Indices of the operator,
# and append them to the function arguments.
function ITensors.op(
        on::OpName, st::SiteType"Osc", s1::Index, s_tail::Index...; kwargs...
    )
    rs = reverse((s1, s_tail...))
    dims = ITensors.dim.(rs)
    opmat = ITensors.op(on, st, dims...; kwargs...)
    return ITensors.itensor(opmat, prime.(rs)..., dag.(rs)...)
end

# Operators
# ---------
function ITensors.op(::OpName"Asum", st::SiteType"Osc", d::Int)
    return ITensors.op(OpName("A"), st, d) + ITensors.op(OpName("Adag"), st, d)
end
function ITensors.op(on::OpName"X", st::SiteType"Osc", d::Int)
    return ITensors.op(OpName("Asum"), st, d) / sqrt(2)
end
function ITensors.op(on::OpName"Y", st::SiteType"Osc", d::Int)
    return im/sqrt(2)*(ITensors.op(OpName("Adag"), st, d) - ITensors.op(OpName("A"), st, d))
end

ITensors.alias(::OpName"a-") = OpName"A"()
ITensors.alias(::OpName"a+") = OpName"Adag"()
ITensors.alias(::OpName"asum") = OpName"Asum"()

ITensors.op(on::OpName"a+", st::SiteType"Osc", d::Int) = ITensors.op(ITensors.alias(on), st, d)
ITensors.op(on::OpName"a-", st::SiteType"Osc", d::Int) = ITensors.op(ITensors.alias(on), st, d)
ITensors.op(on::OpName"asum", st::SiteType"Osc", d::Int) = ITensors.op(ITensors.alias(on), st, d)

# Projection on the eigenstates of the number operator
# ------------------------------------------------

# TODO: this is a really basic function... we could do without this.
"""
    osc_levels_proj(s::Index{Int64}, n::Int)

Return an MPS representing the `n`-th occupied level of `s`'s SiteType.
"""
function osc_levels_proj(site::Index{Int64}, level::Int)
  st = state(site, "$level")
  return MPS([st])
end

# Choice of the oscillator's initial state
# ----------------------------------------

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
    throw(DomainError(statename,
                      "Unrecognised state name; please choose from "*
                      "\"empty\", \"fockN\" or \"thermal\"."))
  end
  return MPS([s])
end
