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

# An element A ∈ Mat(ℂ⁴) is representeb by the a vector v such that
#     vᵢ = tr(Λᵢ A),
# while a linear map L : Mat(ℂ⁴) → Mat(ℂ⁴) by the matrix ℓ such that
#     ℓᵢⱼ = tr(Λᵢ L(Λⱼ)).

# Shorthand notation:
elop(on::AbstractString) = ITensors.op(OpName(on), SiteType("Electron"))
function vstate(sn::AbstractString, ::SiteType"vElectron")
    v = ITensors.state(StateName(sn), SiteType("Electron"))
    return vec(kron(v, v'), gellmannbasis(4))
end
function vop(on::AbstractString, ::SiteType"vElectron")
    return vec(ITensors.op(OpName(on), SiteType("Electron")), gellmannbasis(4))
end
premul(mat, ::SiteType"vElectron") = vec(x -> mat * x, gellmannbasis(4))
postmul(mat, ::SiteType"vElectron") = vec(x -> x * mat, gellmannbasis(4))

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
# If
#   L : Mat(ℂ⁴) ⊗ Mat(ℂ⁴) → Mat(ℂ⁴) ⊗ Mat(ℂ⁴)
# can be written as
#   L₁ ⊗ L₂ for Lᵢ : Mat(ℂ⁴) → Mat(ℂ⁴)
# then
#   ⟨êᵢ₁ ⊗ êᵢ₂, L(êⱼ₁ ⊗ êⱼ₂)⟩ = ⟨êᵢ₁, L₁(êⱼ₁)⟩ ⟨êᵢ₂, L₂(êⱼ₂)⟩.

function ITensors.op(::OpName"Id", ::SiteType"vElectron")
    return Matrix(1.0I, 16, 16)
end
ITensors.op(::OpName"⋅Id", st::SiteType"vElectron") = ITensors.op(OpName("Id"), st)
ITensors.op(::OpName"Id⋅", st::SiteType"vElectron") = ITensors.op(OpName("Id"), st)

# Number operators
ITensors.op(::OpName"Nup⋅", st::SiteType"vElectron") = premul(elop("Nup"), st)
ITensors.op(::OpName"⋅Nup", st::SiteType"vElectron") = postmul(elop("Nup"), st)

ITensors.op(::OpName"Ndn⋅", st::SiteType"vElectron") = premul(elop("Ndn"), st)
ITensors.op(::OpName"⋅Ndn", st::SiteType"vElectron") = postmul(elop("Ndn"), st)

ITensors.op(::OpName"Ntot⋅", st::SiteType"vElectron") = premul(elop("Ntot"), st)
ITensors.op(::OpName"⋅Ntot", st::SiteType"vElectron") = postmul(elop("Ntot"), st)

ITensors.op(::OpName"NupNdn⋅", st::SiteType"vElectron") = premul(elop("NupNdn"), st)
ITensors.op(::OpName"⋅NupNdn", st::SiteType"vElectron") = postmul(elop("NupNdn"), st)

# Jordan-Wigner string operator
ITensors.op(::OpName"F⋅", st::SiteType"vElectron") = premul(elop("F"), st)
ITensors.op(::OpName"⋅F", st::SiteType"vElectron") = postmul(elop("F"), st)

# Creation and annihilation operators (with and without strings, as needed)
ITensors.op(::OpName"Aup⋅", st::SiteType"vElectron") = premul(elop("Aup"), st)
ITensors.op(::OpName"⋅Aup", st::SiteType"vElectron") = postmul(elop("Aup"), st)

ITensors.op(::OpName"Aup†⋅", st::SiteType"vElectron") = premul(elop("Adagup"), st)
ITensors.op(::OpName"⋅Aup†", st::SiteType"vElectron") = postmul(elop("Adagup"), st)

function ITensors.op(::OpName"Aup†F⋅", st::SiteType"vElectron")
    return premul(elop("Adagup") * elop("F"), st)
end
function ITensors.op(::OpName"⋅Aup†F", st::SiteType"vElectron")
    return postmul(elop("Adagup") * elop("F"), st)
end

ITensors.op(::OpName"AupF⋅", st::SiteType"vElectron") = premul(elop("Aup") * elop("F"), st)
ITensors.op(::OpName"⋅AupF", st::SiteType"vElectron") = postmul(elop("Aup") * elop("F"), st)

ITensors.op(::OpName"Adn⋅", st::SiteType"vElectron") = premul(elop("Adn"), st)
ITensors.op(::OpName"⋅Adn", st::SiteType"vElectron") = postmul(elop("Adn"), st)

ITensors.op(::OpName"Adn†⋅", st::SiteType"vElectron") = premul(elop("Adagdn"), st)
ITensors.op(::OpName"⋅Adn†", st::SiteType"vElectron") = postmul(elop("Adagdn"), st)

ITensors.op(::OpName"FAdn⋅", st::SiteType"vElectron") = premul(elop("F") * elop("Adn"), st)
ITensors.op(::OpName"⋅FAdn", st::SiteType"vElectron") = postmul(elop("F") * elop("Adn"), st)

function ITensors.op(::OpName"FAdn†⋅", st::SiteType"vElectron")
    return premul(elop("F") * elop("Adagdn"), st)
end
function ITensors.op(::OpName"⋅FAdn†", st::SiteType"vElectron")
    return postmul(elop("F") * elop("Adagdn"), st)
end

function ITensors.op(::OpName"Adn†F⋅", st::SiteType"vElectron")
    return premul(elop("Adagdn") * elop("F"), st)
end
function ITensors.op(::OpName"⋅Adn†F", st::SiteType"vElectron")
    return postmul(elop("Adagdn") * elop("F"), st)
end
