# Space of harmonic oscillators (vectorised)
# ==========================================
"""
    ITensors.space(st::SiteType"vOsc"; dim = 2)

Create the Hilbert space for a site of type "vOsc", i.e. a vectorised
harmonic oscillator of dimension `dim` (this means that the space has dimension
`dim^2`), where the vectorisation is performed wrt the (Hermitian) Gell-Mann
basis of `Mat(ℂᵈⁱᵐ)`.
"""
function ITensors.space(::SiteType"vOsc"; dim=2)
    return dim^2
end

# Just as with the "Osc" site type, "vOsc" operator and states require that we specify the
# dimension of the space, so we need to compute ITensors.dim.(rs), i.e. the dimensions of
# the Indices of the operator, and append them to the function arguments.

function ITensors.state(sn::StateName, st::SiteType"vOsc", s::Index; kwargs...)
    d = isqrt(ITensors.dim(s))
    stvec = ITensors.state(sn, st, d; kwargs...)
    return ITensors.itensor(stvec, s)
end

function ITensors.op(on::OpName, st::SiteType"vOsc", s1::Index, s_tail::Index...; kwargs...)
    rs = reverse((s1, s_tail...))
    dims = isqrt.(ITensors.dim.(rs))
    # Use `isqrt` which takes an Int and returns an Int; the decimal part is
    # truncated, but it's not a problem since by construction dim.(rs) are
    # perfect squares, so the decimal part should be zero.
    opmat = ITensors.op(on, st, dims...; kwargs...)
    return ITensors.itensor(opmat, prime.(rs)..., dag.(rs)...)
end

# Aliases (for backwards compatibility)
ITensors.alias(::SiteType"HvOsc") = SiteType"vOsc"()
ITensors.alias(::SiteType"vecOsc") = SiteType"vOsc"()

function ITensors.space(st::SiteType"HvOsc"; kwargs...)
    return ITensors.space(ITensors.alias(st); kwargs...)
end
function ITensors.space(st::SiteType"vecOsc"; kwargs...)
    return ITensors.space(ITensors.alias(st); kwargs...)
end

ITensors.val(vn::ValName, st::SiteType"HvOsc") = ITensors.val(vn, ITensors.alias(st))
ITensors.val(vn::ValName, st::SiteType"vecOsc") = ITensors.val(vn, ITensors.alias(st))

function ITensors.state(sn::StateName, st::SiteType"vecOsc", s::Index; kwargs...)
    return ITensors.state(sn, ITensors.alias(st), s; kwargs...)
end
function ITensors.state(sn::StateName, st::SiteType"HvOsc", s::Index; kwargs...)
    return ITensors.state(sn, ITensors.alias(st), s; kwargs...)
end

function ITensors.op(
    on::OpName, st::SiteType"vecOsc", s1::Index, s_tail::Index...; kwargs...
)
    return ITensors.op(on, ITensors.alias(st), s1, s_tail...; kwargs...)
end
function ITensors.op(
    on::OpName, st::SiteType"HvOsc", s1::Index, s_tail::Index...; kwargs...
)
    return ITensors.op(on, ITensors.alias(st), s1, s_tail...; kwargs...)
end

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vOsc", d::Int)
    v = ITensors.state(sn, SiteType("Osc"))
    return PseudomodesTTEDOPA.vec(kron(v, v'), gellmannbasis(d))
end
function vop(sn::StateName, ::SiteType"vOsc", d::Int)
    sn = statenamestring(sn)
    on = sn[1] == 'v' ? sn[2:end] : sn
    return PseudomodesTTEDOPA.vec(try_op(OpName(on), SiteType("Osc"), d), gellmannbasis(d))
end

# States
# ------
function ITensors.state(::StateName{N}, ::SiteType"vOsc", d::Int) where {N}
    # Eigenstates êₙ ⊗ êₙ of the number operator, wrt the Hermitian basis.
    n = parse(Int, String(N))
    v = zeros(d)
    v[n + 1] = 1.0
    return vec(kron(v, v'), gellmannbasis(d))
end

function ITensors.state(
    ::StateName"ThermEq", st::SiteType"vOsc", d::Int; frequency::Real, temperature::Real
)
    if temperature == 0
        return ITensors.state(StateName("0"), st, d)
    else
        numop = ITensors.op(OpName("N"), SiteType("Osc"), d)
        # We don't need to define our own matrix for the number operator when
        # we can call this one instead.
        ρ_eq = exp(-frequency / temperature * numop)
        ρ_eq /= tr(ρ_eq)
        return vec(ρ_eq, gellmannbasis(d))
    end
end

# Product of X = (a+a†)/√2 and of the thermal equilibrium state Z⁻¹vec(exp(-βH)).
# It is used in the computation of the correlation function of the bath.
function ITensors.state(
    ::StateName"X⋅Therm", st::SiteType"vOsc", d::Int; frequency::Real, temperature::Real
)
    xop = ITensors.op(OpName("X"), SiteType("Osc"), d)
    if temperature == 0
        ρ_eq = zeros(Float64, d, d)
        ρ_eq[1, 1] = 1.0
    else
        numop = ITensors.op(OpName("N"), SiteType("Osc"), d)
        ρ_eq = exp(-frequency / temperature * numop)
        ρ_eq /= tr(ρ_eq)
    end
    return vec(xop * ρ_eq, gellmannbasis(d))
end

# States representing vectorised operators
# ----------------------------------------
ITensors.state(sn::StateName"vAdag", st::SiteType"vOsc", d::Int) = vop(sn, st, d)
ITensors.state(sn::StateName"vA", st::SiteType"vOsc", d::Int) = vop(sn, st, d)
ITensors.state(sn::StateName"vN", st::SiteType"vOsc", d::Int) = vop(sn, st, d)
ITensors.state(sn::StateName"vId", st::SiteType"vOsc", d::Int) = vop(sn, st, d)
ITensors.state(sn::StateName"vX", st::SiteType"vOsc", d::Int) = vop(sn, st, d)
ITensors.state(sn::StateName"vY", st::SiteType"vOsc", d::Int) = vop(sn, st, d)

# Aliases (for backwards compatibility)
function ITensors.state(::StateName"veca+", st::SiteType"vOsc", d::Int)
    return ITensors.state(StateName("vAdag"), st, d)
end
function ITensors.state(::StateName"veca-", st::SiteType"vOsc", d::Int)
    return ITensors.state(StateName("vA"), st, d)
end
function ITensors.state(::StateName"vecplus", st::SiteType"vOsc", d::Int)
    return ITensors.state(StateName("vAdag"), st, d)
end
function ITensors.state(::StateName"vecminus", st::SiteType"vOsc", d::Int)
    return ITensors.state(StateName("vA"), st, d)
end
function ITensors.state(::StateName"vecN", st::SiteType"vOsc", d::Int)
    return ITensors.state(StateName("vN"), st, d)
end
function ITensors.state(::StateName"vecId", st::SiteType"vOsc", d::Int)
    return ITensors.state(StateName("vId"), st, d)
end

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vOsc", d::Int)
    return PseudomodesTTEDOPA.vec(x -> mat * x, gellmannbasis(d))
end
function postmultiply(mat, ::SiteType"vOsc", d::Int)
    return PseudomodesTTEDOPA.vec(x -> x * mat, gellmannbasis(d))
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the Osc site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vOsc", d::Int; kwargs...)
    name = strip(String(ITensors.name(on))) # Remove extra whitespace
    if name == "Id"
        return Matrix(1.0I, d^2, d^2)
    end
    dotloc = findfirst("⋅", name)
    # This returns the position of the cdot in the operator name String.
    # It is `nothing` if no cdot is found.
    if !isnothing(dotloc)
        on1, on2 = nothing, nothing
        on1 = name[1:prevind(name, dotloc.start)]
        on2 = name[nextind(name, dotloc.start):end]
        # If the OpName `on` is written correctly, i.e. it is either "A⋅" or "⋅A" for some
        # A, then either `on1` or `on2` has to be empty (not both, not neither of them).
        if (on1 == "" && on2 == "") || (on1 != "" && on2 != "")
            throw(
                ArgumentError(
                    "Invalid operator name: $name. Operator name is not \"Id\" or of the " *
                    "form \"A⋅\" or \"⋅A\"",
                ),
            )
        end
        # name == "⋅A" -> on1 is an empty string
        # name == "A⋅" -> on2 is an empty string
        if on1 == ""
            mat = try_op(OpName(on2), SiteType("Osc"), d; kwargs...)
            return postmultiply(mat, st, d)
        elseif on2 == ""
            mat = try_op(OpName(on1), SiteType("Osc"), d; kwargs...)
            return premultiply(mat, st, d)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
end

# GKSL equation terms
# -------------------
# For reference: given an Index `i`,
#   op(i, "B * A") == replaceprime(op(i', "B") * op(i, "A"), 2, 1)
# so the composition is performed from right to left.
# Look out for the correct order of operations:
#   vec(ABx) == op(A⋅) * vec(Bx) == op(A⋅) * (op(B⋅) * vec(x))
#   vec(xAB) == op(⋅B) * vec(xA) == op(⋅B) * (op(⋅A) * vec(x))

# Separate absorption and dissipation terms in GKSL equation
function dissipator_gain(n::Int)
    dsp = OpSum()
    # a† ρ a - ½ a a† ρ - ½ ρ a a†
    dsp += "Adag⋅ * ⋅A", n
    dsp += -0.5, "A⋅ * Adag⋅", n
    dsp += -0.5, "⋅Adag * ⋅A", n
    return dsp
end

function dissipator_loss(n::Int)
    dsp = OpSum()
    # a ρ a† - ½ a† a ρ - ½ ρ a† a
    dsp += "A⋅ * ⋅Adag", n
    dsp += -0.5, "N⋅", n
    dsp += -0.5, "⋅N", n
    return dsp
end

"""
    dissipator(n::Int, frequency::Real, temperature::Real)

Return an OpSum object which represents the dissipator in the GKSL equation for given
`frequency` and `temperature`, i.e. with dissipation coefficient `1/(exp(freq/temp) - 1)`.
"""
function dissipator(n::Int, frequency::Real, temperature::Real)
    if temperature == 0
        return dissipator_loss(n)
    else
        # Use `expm1(x)` instead of `e^x-1` for better precision.
        avgn = 1 / expm1(frequency / temperature)
        return (avgn + 1) * dissipator_loss(n) + avgn * dissipator_gain(n)
    end
end

"""
    mixedlindbladplus(n1::Int, n2::Int)

Return on OpSum expression which represents a mixed dissipation operator, appearing in the
equation for two pseudomodes, on sites at positions `n1` and `n2` in the system.
"""
function mixedlindbladplus(n1::Int, n2::Int)
    x = OpSum()
    x += "A⋅", n1, "⋅Adag", n2
    x += "⋅Adag", n1, "A⋅", n2
    x += -0.5, "Adag⋅", n1, "A⋅", n2
    x += -0.5, "A⋅", n1, "Adag⋅", n2
    x += -0.5, "⋅Adag", n1, "⋅A", n2
    x += -0.5, "⋅A", n1, "⋅Adag", n2
    return x
end

"""
    mixedlindbladminus(n1::Int, n2::Int)

Return on OpSum expression which represents a mixed dissipation operator, appearing in the
equation for two pseudomodes, on sites at positions `n1` and `n2` in the system.
"""
function mixedlindbladminus(n1::Int, n2::Int)
    x = OpSum()
    x += "Adag⋅", n1, "⋅A", n2
    x += "Adag⋅", n2, "⋅A", n1
    x += -0.5, "A⋅", n1, "Adag⋅", n2
    x += -0.5, "Adag⋅", n1, "A⋅", n2
    x += -0.5, "⋅A", n1, "⋅Adag", n2
    x += -0.5, "⋅Adag", n1, "⋅A", n2
    return x
end
