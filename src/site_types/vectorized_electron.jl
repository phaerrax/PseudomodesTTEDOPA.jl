# Space of electrons (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vElectron")

Create the Hilbert space for a site of type "vElectron", i.e. a mixed state describing a
site with a 1/2-spin degree of freedom.
The density matrix is represented in the generalised Gell-Mann basis of `Mat(ℂ⁴)`, composed
of Hermitian traceless matrices together with the identity matrix.
"""
ITensors.space(::SiteType"vElectron") = 16

# Shorthand notation:
elop(on::AbstractString) = ITensors.op(OpName(on), SiteType("Electron"))
function vstate(sn::AbstractString, ::SiteType"vElectron")
    v = ITensors.state(StateName(sn), SiteType("Electron"))
    return vec(kron(v, v'), gellmannbasis(4))
end
function vop(on::AbstractString, ::SiteType"vElectron")
    return vec(ITensors.op(OpName(on), SiteType("Electron")), gellmannbasis(4))
end

# States (actual ones)
# --------------------
ITensors.state(::StateName"Emp", st::SiteType"vElectron") = vstate("Emp", st)
ITensors.state(::StateName"Up", st::SiteType"vElectron") = vstate("Up", st)
ITensors.state(::StateName"Dn", st::SiteType"vElectron") = vstate("Dn", st)
ITensors.state(::StateName"UpDn", st::SiteType"vElectron") = vstate("UpDn", st)

# States representing vectorised operators
# ----------------------------------------
function ITensors.op(::OpName"Id", ::SiteType"Electron") # ITensors doesn't define this
    return Matrix(1.0I, 4, 4)
end

ITensors.state(::StateName"vId", st::SiteType"vElectron") = vop("Id", st)
ITensors.state(::StateName"vecId", st::SiteType"vElectron") = vop("Id", st)
ITensors.state(::StateName"vNup", st::SiteType"vElectron") = vop("Nup", st)
ITensors.state(::StateName"vNdn", st::SiteType"vElectron") = vop("Ndn", st)
ITensors.state(::StateName"vNtot", st::SiteType"vElectron") = vop("Ntot", st)
ITensors.state(::StateName"vNupNdn", st::SiteType"vElectron") = vop("NupNdn", st)

# Operators acting on vectorised spins
# ------------------------------------
function premultiply(mat, ::SiteType"vElectron")
    return PseudomodesTTEDOPA.vec(x -> mat * x, gellmannbasis(4))
end
function postmultiply(mat, ::SiteType"vElectron")
    return PseudomodesTTEDOPA.vec(x -> x * mat, gellmannbasis(4))
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the Electron site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vElectron"; kwargs...)
    name = strip(String(ITensors.name(on))) # Remove extra whitespace
    if name == "Id"
        return Matrix(1.0I, 16, 16)
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
            mat = try_op(OpName(on2), SiteType("Electron"); kwargs...)
            return postmultiply(mat, st)
        elseif on2 == ""
            mat = try_op(OpName(on1), SiteType("Electron"); kwargs...)
            return premultiply(mat, st)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
end
