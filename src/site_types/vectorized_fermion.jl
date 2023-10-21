# Space of spin-1/2 particles (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vFermion")

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

function vstate(sn::AbstractString, ::SiteType"vFermion")
    v = ITensors.state(StateName(sn), SiteType("Fermion"))
    return PseudomodesTTEDOPA.vec(kron(v, v'), gellmannbasis(2))
end
function vop(on::AbstractString, ::SiteType"vFermion")
    return PseudomodesTTEDOPA.vec(
        ITensors.op(OpName(on), SiteType("Fermion")), gellmannbasis(2)
    )
end

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
function ITensors.state(::StateName"vC", ::SiteType"vFermion")
    # FIXME: add Jordan-Wigner strings
    return vec(ITensors.op(OpName("S-"), SiteType("S=1/2")), gellmannbasis(2))
end
function ITensors.state(::StateName"vCdag", ::SiteType"vFermion")
    # FIXME: add Jordan-Wigner strings
    return vec(ITensors.op(OpName("S+"), SiteType("S=1/2")), gellmannbasis(2))
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

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vFermion")
    return PseudomodesTTEDOPA.vec(x -> mat * x, gellmannbasis(2))
end
function postmultiply(mat, ::SiteType"vFermion")
    return PseudomodesTTEDOPA.vec(x -> x * mat, gellmannbasis(2))
end

function ITensors.op(::OpName"Id", ::SiteType"Fermion")
    return Matrix(1.0I, 2, 2)
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the Fermion site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vFermion"; kwargs...)
    name = strip(String(ITensors.name(on))) # Remove extra whitespace
    if name == "Id"
        return Matrix(1.0I, 4, 4)
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
            mat = try_op(OpName(on2), SiteType("Fermion"); kwargs...)
            return postmultiply(mat, st)
        elseif on2 == ""
            mat = try_op(OpName(on1), SiteType("Fermion"); kwargs...)
            return premultiply(mat, st)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
end
