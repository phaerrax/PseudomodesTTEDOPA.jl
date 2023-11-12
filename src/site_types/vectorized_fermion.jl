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

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vFermion")
    v = ITensors.state(sn, SiteType("Fermion"))
    return PseudomodesTTEDOPA.vec(kron(v, v'), gellmannbasis(2))
end
function vop(sn::StateName, ::SiteType"vFermion")
    sn = statenamestring(sn)
    on = sn[1] == 'v' ? sn[2:end] : sn
    return PseudomodesTTEDOPA.vec(try_op(OpName(on), SiteType("Fermion")), gellmannbasis(2))
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"Emp", st::SiteType"vFermion") = vstate(sn, st)
ITensors.state(sn::StateName"Occ", st::SiteType"vFermion") = vstate(sn, st)

function ITensors.state(::StateName"Up", st::SiteType"vFermion")
    return ITensors.state(StateName("Occ"), st)
end
function ITensors.state(::StateName"Dn", st::SiteType"vFermion")
    return ITensors.state(StateName("Emp"), st)
end

# States representing vectorised operators
# ----------------------------------------
ITensors.state(sn::StateName"vId", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"vN", st::SiteType"vFermion") = vop(sn, st)

ITensors.state(sn::StateName"vA", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"va", st::SiteType"vFermion") = vop(sn, st)

ITensors.state(sn::StateName"vAdag", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"vadag", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"vA†", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"va†", st::SiteType"vFermion") = vop(sn, st)

function ITensors.state(::StateName"vecId", ::SiteType"vFermion")
    return ITensors.state(StateName("vId"), SiteType("vFermion"))
end

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vFermion")
    return PseudomodesTTEDOPA.vec(x -> mat * x, gellmannbasis(2))
end
function postmultiply(mat, ::SiteType"vFermion")
    return PseudomodesTTEDOPA.vec(x -> x * mat, gellmannbasis(2))
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
