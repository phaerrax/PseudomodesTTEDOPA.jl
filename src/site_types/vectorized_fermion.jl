# Space of spin-1/2 particles (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vFermion"; dim = 2)

Create the Hilbert space for a site of type "vFermion", i.e. a vectorised
spin-1/2 particle, where the vectorisation is performed wrt the generalised
Gell-Mann basis of `Mat(ℂ²)`, composed of Hermitian traceless matrices
together with the identity matrix.
"""
ITensors.space(::SiteType"vFermion") = 4

# Elements of and operators on Mat(ℂ²) are expanded wrt the basis {Λᵢ}ᵢ₌₁⁴ of
# generalised Gell-Mann matrices (plus a multiple of the identity).
# An element A ∈ Mat(ℂ²) is representeb by the a vector v such that
#     vᵢ = tr(Λᵢ A),
# while a linear map L : Mat(ℂ²) → Mat(ℂ²) by the matrix ℓ such that
#     ℓᵢⱼ = tr(Λᵢ L(Λⱼ)).

# States
# ------

# "Up" = ê₊ ⊗ ê₊' = ⎡1⎤ ⊗ [1 0]
#                   ⎣0⎦
# "Dn" = ê₋ ⊗ ê₋' = ⎡0⎤ ⊗ [0 1]
#                   ⎣1⎦
# We don't use them explicitly, however, but we refer to how ITensors implements
# them; this way, we are sure that the operators act correctly on these states.
function ITensors.state(sn::StateName"Up", ::SiteType"vFermion")
    v = ITensors.state(sn, SiteType("S=1/2"))
    return vec(kron(v, v'), gellmannbasis(2))
end
function ITensors.state(sn::StateName"Dn", ::SiteType"vFermion")
    v = ITensors.state(sn, SiteType("S=1/2"))
    return vec(kron(v, v'), gellmannbasis(2))
end

# States representing vectorised operators
# ----------------------------------------
function ITensors.state(::StateName"vSx", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("Sx"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vSy", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("Sy"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vSz", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("Sz"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.state(::StateName"vσx", st::SiteType"vFermion")
    return 2 * ITensors.state(StateName("vSx"), st)
end
function ITensors.state(::StateName"vσy", st::SiteType"vFermion")
    return 2 * ITensors.state(StateName("vSy"), st)
end
function ITensors.state(::StateName"vσz", st::SiteType"vFermion")
    return 2 * ITensors.state(StateName("vSz"), st)
end

function ITensors.state(::StateName"vId", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("Id"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vN", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("N"), SiteType("S=1/2")), gellmannbasis(2))
end

# Aliases (for backwards compatibility)
function ITensors.state(::StateName"vecσx", st::SiteType"vFermion")
    return ITensors.state(StateName("vσx"), st)
end
function ITensors.state(::StateName"vecσy", st::SiteType"vFermion")
    return ITensors.state(StateName("vσy"), st)
end
function ITensors.state(::StateName"vecσz", st::SiteType"vFermion")
    return ITensors.state(StateName("vσz"), st)
end
function ITensors.state(::StateName"vecplus", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("S+"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vecminus", ::SiteType"vFermion")
    return vec(ITensors.op(OpName("S-"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vecN", st::SiteType"vFermion")
    return ITensors.state(StateName("vN"), st)
end
function ITensors.state(::StateName"vecId", st::SiteType"vFermion")
    return ITensors.state(StateName("vId"), st)
end

# Operators acting on vectorised spins
# ------------------------------------
# Luckily, even when they are acting on two sites at the same times, every
# operator we need to define is factorised (or a sum of factorised operators).
# This simplifies the calculations immensely: if
#   L : Mat(ℂ²) ⊗ Mat(ℂ²) → Mat(ℂ²) ⊗ Mat(ℂ²)
# can be written as L₁ ⊗ L₂ for Lᵢ : Mat(ℂ²) → Mat(ℂ²) then
#   ⟨êᵢ₁ ⊗ êᵢ₂, L(êⱼ₁ ⊗ êⱼ₂)⟩ = ⟨êᵢ₁, L₁(êⱼ₁)⟩ ⟨êᵢ₂, L₂(êⱼ₂)⟩.

function ITensors.op(::OpName"Id", ::SiteType"vFermion")
    return Matrix(1.0I, 4, 4)
end
ITensors.op(::OpName"⋅Id", st::SiteType"vFermion") = ITensors.op(OpName("Id"), st)
ITensors.op(::OpName"Id⋅", st::SiteType"vFermion") = ITensors.op(OpName("Id"), st)

# N == σ⁺σ⁻
function ITensors.op(::OpName"N⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("N"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅N", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("N"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σ+⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("S+"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σ+", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("S+"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σ-⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("S-"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σ-", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("S-"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σx⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("σx"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σx", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("σx"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σy⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("σy"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σy", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("σy"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"σz⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("σz"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅σz", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("σz"), SiteType("S=1/2")), gellmannbasis(2))
end

function ITensors.op(s::OpName"F⋅", ::SiteType"vFermion")
    return vec(x -> ITensors.op(OpName("F"), SiteType("S=1/2")) * x, gellmannbasis(2))
end
function ITensors.op(::OpName"⋅F", ::SiteType"vFermion")
    return vec(x -> x * ITensors.op(OpName("F"), SiteType("S=1/2")), gellmannbasis(2))
end
