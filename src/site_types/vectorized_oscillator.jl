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
function ITensors.state(::StateName"vAdag", ::SiteType"vOsc", d::Int)
    return vec(ITensors.op("Adag", SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.state(::StateName"vA", ::SiteType"vOsc", d::Int)
    return vec(ITensors.op("A", SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.state(::StateName"vN", ::SiteType"vOsc", d::Int)
    return vec(ITensors.op("N", SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.state(::StateName"vId", ::SiteType"vOsc", d::Int)
    return vec(ITensors.op("Id", SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.state(::StateName"vX", ::SiteType"vOsc", d::Int)
    return vec(ITensors.op("X", SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.state(::StateName"vY", ::SiteType"vOsc", d::Int)
    return vec(ITensors.op("Y", SiteType("Osc"), d), gellmannbasis(d))
end

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

# Operators acting on vectorised oscillators
# ------------------------------------------
function ITensors.op(::OpName"Id", ::SiteType"vOsc", d::Int)
    # The basis is orthonormal wrt the trace, therefore
    #   tr(bi' * id(bj)) == tr(bi' * bj) == delta(i,j).
    # It's the identity matrix whatever the basis, it's useless to compute its
    # vectorization.
    return Matrix(1.0I, d^2, d^2)
end
function ITensors.op(::OpName"⋅Id", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("Id"), st, d)
end
function ITensors.op(::OpName"Id⋅", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("Id"), st, d)
end

function ITensors.op(::OpName"⋅Adag", ::SiteType"vOsc", d::Int)
    return vec(x -> x * ITensors.op(OpName("Adag"), SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.op(::OpName"Adag⋅", ::SiteType"vOsc", d::Int)
    return vec(x -> ITensors.op(OpName("Adag"), SiteType("Osc"), d) * x, gellmannbasis(d))
end

function ITensors.op(::OpName"⋅A", ::SiteType"vOsc", d::Int)
    return vec(x -> x * ITensors.op(OpName("A"), SiteType("Osc"), d), gellmannbasis(d))
end
function ITensors.op(::OpName"A⋅", ::SiteType"vOsc", d::Int)
    return vec(x -> ITensors.op(OpName("A"), SiteType("Osc"), d) * x, gellmannbasis(d))
end

function ITensors.op(::OpName"Asum⋅", ::SiteType"vOsc", d::Int)
    return vec(x -> ITensors.op(OpName("Asum"), SiteType("Osc"), d) * x, gellmannbasis(d))
end
function ITensors.op(::OpName"⋅Asum", ::SiteType"vOsc", d::Int)
    return vec(x -> x * ITensors.op(OpName("Asum"), SiteType("Osc"), d), gellmannbasis(d))
end

function ITensors.op(::OpName"N⋅", ::SiteType"vOsc", d::Int)
    return vec(x -> ITensors.op(OpName("N"), SiteType("Osc"), d) * x, gellmannbasis(d))
end
function ITensors.op(::OpName"⋅N", ::SiteType"vOsc", d::Int)
    return vec(x -> x * ITensors.op(OpName("N"), SiteType("Osc"), d), gellmannbasis(d))
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
        # Use `expm(x)` instead of `e^x-1` for better precision.
        avgn = 1 / expm(frequency / temperature)
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

# Aliases (for backwards compatibility)
function ITensors.op(::OpName"⋅a+", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("⋅Adag"), st, d)
end
function ITensors.op(::OpName"a+⋅", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("Adag⋅"), st, d)
end
function ITensors.op(::OpName"⋅a-", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("⋅A"), st, d)
end
function ITensors.op(::OpName"a-⋅", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("A⋅"), st, d)
end
function ITensors.op(::OpName"⋅asum", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("⋅Asum"), st, d)
end
function ITensors.op(::OpName"asum⋅", st::SiteType"vOsc", d::Int)
    return ITensors.op(OpName("Asum⋅"), st, d)
end
