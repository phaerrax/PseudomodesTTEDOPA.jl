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

function ITensors.op(on::OpName, st::SiteType"vecS=1/2"; kwargs...)
    return ITensors.op(on, ITensors.alias(st); kwargs...)
end
function ITensors.op(on::OpName, st::SiteType"HvS=1/2"; kwargs...)
    return ITensors.op(on, ITensors.alias(st); kwargs...)
end

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vS=1/2")
    v = ITensors.state(sn, SiteType("S=1/2"))
    return PseudomodesTTEDOPA.vec(kron(v, v'), gellmannbasis(2))
end
function vop(sn::StateName, ::SiteType"vS=1/2")
    sn = statenamestring(sn)
    on = sn[1] == 'v' ? sn[2:end] : sn
    return PseudomodesTTEDOPA.vec(try_op(OpName(on), SiteType("S=1/2")), gellmannbasis(2))
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"Up", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"Dn", st::SiteType"vS=1/2") = vstate(sn, st)

# States representing vectorised operators
# ----------------------------------------
ITensors.state(sn::StateName"vSx", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"vSy", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"vSz", st::SiteType"vS=1/2") = vop(sn, st)

ITensors.state(sn::StateName"vσx", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"vσy", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"vσz", st::SiteType"vS=1/2") = vop(sn, st)

ITensors.state(sn::StateName"vId", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"vN", st::SiteType"vS=1/2") = vop(sn, st)

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
function ITensors.state(::StateName"vecId", st::SiteType"vS=1/2")
    return ITensors.state(StateName("vId"), st)
end

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vS=1/2")
    return PseudomodesTTEDOPA.vec(x -> mat * x, gellmannbasis(2))
end
function postmultiply(mat, ::SiteType"vS=1/2")
    return PseudomodesTTEDOPA.vec(x -> x * mat, gellmannbasis(2))
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the S=1/2 site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vS=1/2"; kwargs...)
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
            mat = try_op(OpName(on2), SiteType("S=1/2"); kwargs...)
            return postmultiply(mat, st)
        elseif on2 == ""
            mat = try_op(OpName(on1), SiteType("S=1/2"); kwargs...)
            return premultiply(mat, st)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
end
