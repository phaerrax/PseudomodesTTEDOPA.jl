# Space of spin-1/2 particles (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vS=1/2"; dim = 2)

Create the Hilbert space for a site of type "vS=1/2", i.e. a vectorised
spin-1/2 particle, where the vectorisation is performed wrt the generalised
Gell-Mann basis of `Mat(ℂ²)`, composed of Hermitian traceless matrices
together with the identity matrix.
"""
ITensors.space(::SiteType"vS=1/2") = 4

# Elements of and operators on Mat(ℂ²) are expanded wrt the basis {Λᵢ}ᵢ₌₁⁴ of
# generalised Gell-Mann matrices (plus a multiple of the identity).
# An element A ∈ Mat(ℂ²) is representeb by the a vector v such that
#     vᵢ = tr(Λᵢ A),
# while a linear map L : Mat(ℂ²) → Mat(ℂ²) by the matrix ℓ such that
#     ℓᵢⱼ = tr(Λᵢ L(Λⱼ)).

# Aliases (for backwards compatibility)
ITensors.alias(::SiteType"HvS=1/2") = SiteType"vS=1/2"()
ITensors.alias(::SiteType"vecS=1/2") = SiteType"vS=1/2"()

ITensors.space(st::SiteType"HvS=1/2") = ITensors.space(ITensors.alias(st))
ITensors.space(st::SiteType"vecS=1/2") = ITensors.space(ITensors.alias(st))

ITensors.val(vn::ValName, st::SiteType"HvS=1/2") = ITensors.val(vn, ITensors.alias(st))
ITensors.val(vn::ValName, st::SiteType"vecS=1/2") = ITensors.val(vn, ITensors.alias(st))

# !
# We need to replicate the signatures of the functions below: since the states are defined
# using ITensors.state(::StateName, ::SiteType"vS=1/2"), the aliasing function must follow
# the same structure, otherwise it doesn't pick up the definitions for "vS=1/2".
# The same goes for the operators.
function ITensors.state(sn::StateName, st::SiteType"vecS=1/2"; kwargs...)
    return ITensors.state(sn, ITensors.alias(st); kwargs...)
end
function ITensors.state(sn::StateName, st::SiteType"HvS=1/2"; kwargs...)
    return ITensors.state(sn, ITensors.alias(st); kwargs...)
end

function ITensors.op(
    on::OpName, st::SiteType"vecS=1/2"; kwargs...
)
    return ITensors.op(on, ITensors.alias(st); kwargs...)
end
function ITensors.op(
    on::OpName, st::SiteType"HvS=1/2"; kwargs...
)
    return ITensors.op(on, ITensors.alias(st); kwargs...)
end

# States
# ------

# "Up" = ê₊ ⊗ ê₊' = ⎡1⎤ ⊗ [1 0]
#                   ⎣0⎦
# "Dn" = ê₋ ⊗ ê₋' = ⎡0⎤ ⊗ [0 1]
#                   ⎣1⎦
function ITensors.state(::StateName"Up", ::SiteType"vS=1/2")
    return vec(kron([1; 0], [1; 0]'), gellmannbasis(2))
end
function ITensors.state(::StateName"Dn", ::SiteType"vS=1/2")
    return vec(kron([0; 1], [0; 1]'), gellmannbasis(2))
end

# States representing vectorised operators
# ----------------------------------------
function ITensors.state(::StateName"vSx", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("Sx"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vSy", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("Sy"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vSz", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("Sz"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.state(::StateName"vσx", st::SiteType"vS=1/2")
    return 2 * ITensors.state(StateName("vSx"), st)
end
function ITensors.state(::StateName"vσy", st::SiteType"vS=1/2")
    return 2 * ITensors.state(StateName("vSy"), st)
end
function ITensors.state(::StateName"vσz", st::SiteType"vS=1/2")
    return 2 * ITensors.state(StateName("vSz"), st)
end

function ITensors.state(::StateName"vId", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("Id"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vN", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("N"), SiteType("S=1/2")), gellmannbasis(2))
end

# Aliases (for backwards compatibility)
function ITensors.state(::StateName"vecσx", st::SiteType"vS=1/2")
    return ITensors.state(StateName("vσx"), st)
end
function ITensors.state(::StateName"vecσy", st::SiteType"vS=1/2")
    return ITensors.state(StateName("vσy"), st)
end
function ITensors.state(::StateName"vecσz", st::SiteType"vS=1/2")
    return ITensors.state(StateName("vσz"), st)
end
function ITensors.state(::StateName"vecplus", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("S+"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vecminus", ::SiteType"vS=1/2")
    return vec(ITensors.op(OpName("S-"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vecN", st::SiteType"vS=1/2")
    return ITensors.state(StateName("vN"), st)
end

# Operators acting on vectorised spins
# ------------------------------------
# Luckily, even when they are acting on two sites at the same times, every
# operator we need to define is factorised (or a sum of factorised operators).
# This simplifies the calculations immensely: if
#   L : Mat(ℂ²) ⊗ Mat(ℂ²) → Mat(ℂ²) ⊗ Mat(ℂ²)
# can be written as L₁ ⊗ L₂ for Lᵢ : Mat(ℂ²) → Mat(ℂ²) then
#   ⟨êᵢ₁ ⊗ êᵢ₂, L(êⱼ₁ ⊗ êⱼ₂)⟩ = ⟨êᵢ₁, L₁(êⱼ₁)⟩ ⟨êᵢ₂, L₂(êⱼ₂)⟩.

function ITensors.op(::OpName"Id", ::SiteType"vS=1/2")
    return Matrix(1.0I, 4, 4)
end
ITensors.op(::OpName"⋅Id", st::SiteType"vS=1/2") = ITensors.op(OpName("Id"), st)
ITensors.op(::OpName"Id⋅", st::SiteType"vS=1/2") = ITensors.op(OpName("Id"), st)

# N == σ⁺σ⁻
function ITensors.op(::OpName"N⋅", ::SiteType"vS=1/2")
    return vec(x -> ITensors.op(OpName("N"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅N", ::SiteType"vS=1/2")
    return vec(x -> x * ITensors.op(OpName("N"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σ+⋅", ::SiteType"vS=1/2")
    return vec(x -> ITensors.op(OpName("S+"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σ+", ::SiteType"vS=1/2")
    return vec(x -> x * ITensors.op(OpName("S+"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σ-⋅", ::SiteType"vS=1/2")
    return vec(x -> ITensors.op(OpName("S-"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σ-", ::SiteType"vS=1/2")
    return vec(x -> x * ITensors.op(OpName("S-"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σx⋅", ::SiteType"vS=1/2")
    return vec(x -> ITensors.op(OpName("σx"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σx", ::SiteType"vS=1/2")
    return vec(x -> x * ITensors.op(OpName("σx"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σy⋅", ::SiteType"vS=1/2")
    return vec(x -> ITensors.op(OpName("σy"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σy", ::SiteType"vS=1/2")
    return vec(x -> x * ITensors.op(OpName("σy"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σz⋅", ::SiteType"vS=1/2")
    return vec(x -> ITensors.op(OpName("σz"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σz", ::SiteType"vS=1/2")
    return vec(x -> x * ITensors.op(OpName("σz"), SiteType("S=1/2")), gellmannbasis(2))
end
