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
function vstate(sn::AbstractString)
    v = ITensors.state(StateName(sn), SiteType("Electron"))
    return vec(kron(v,v'), gellmannbasis(4))
end
vop(on::AbstractString) = vec(ITensors.op(OpName(on), SiteType("Electron")), gellmannbasis(4))
premul(y) = vec(x -> y * x, gellmannbasis(4))
postmul(y) = vec(x -> x * y, gellmannbasis(4))

# States (actual ones)
# --------------------
ITensors.state(::StateName"Emp", ::SiteType"vElectron") = vstate("Emp")
ITensors.state(::StateName"Up", ::SiteType"vElectron") = vstate("Up")
ITensors.state(::StateName"Dn", ::SiteType"vElectron") = vstate("Dn")
ITensors.state(::StateName"UpDn", ::SiteType"vElectron") = vstate("UpDn")

# States representing vectorised operators
# ----------------------------------------
function ITensors.op(::OpName"Id", ::SiteType"Electron") # ITensors doesn't define this
    return Matrix(1.0I, 4, 4)
end

ITensors.state(::StateName"vId", ::SiteType"vElectron") = vop("Id")
ITensors.state(::StateName"vecId", ::SiteType"vElectron") = vop("Id")
ITensors.state(::StateName"vNup", ::SiteType"vElectron") = vop("Nup")
ITensors.state(::StateName"vNdn", ::SiteType"vElectron") = vop("Ndn")
ITensors.state(::StateName"vNtot", ::SiteType"vElectron") = vop("Ntot")

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
ITensors.op(::OpName"Nup⋅", ::SiteType"vElectron") = premul(elop("Nup"))
ITensors.op(::OpName"⋅Nup", ::SiteType"vElectron") = postmul(elop("Nup"))

ITensors.op(::OpName"Ndn⋅", ::SiteType"vElectron") = premul(elop("Ndn"))
ITensors.op(::OpName"⋅Ndn", ::SiteType"vElectron") = postmul(elop("Ndn"))

ITensors.op(::OpName"Ntot⋅", ::SiteType"vElectron") = premul(elop("Ntot"))
ITensors.op(::OpName"⋅Ntot", ::SiteType"vElectron") = postmul(elop("Ntot"))

# Jordan-Wigner string operator
ITensors.op(::OpName"F⋅", ::SiteType"vElectron") = premul(elop("F"))
ITensors.op(::OpName"⋅F", ::SiteType"vElectron") = postmul(elop("F"))

# Creation and annihilation operators (with and without strings, as needed)
ITensors.op(::OpName"Aup⋅", ::SiteType"vElectron") = premul(elop("Aup"))
ITensors.op(::OpName"⋅Aup", ::SiteType"vElectron") = postmul(elop("Aup"))

ITensors.op(::OpName"Aup†⋅", ::SiteType"vElectron") = premul(elop("Adagup"))
ITensors.op(::OpName"⋅Aup†", ::SiteType"vElectron") = postmul(elop("Adagup"))

ITensors.op(::OpName"Aup†F⋅", ::SiteType"vElectron") =  premul(elop("Adagup") * elop("F"))
ITensors.op(::OpName"⋅Aup†F", ::SiteType"vElectron") = postmul(elop("Adagup") * elop("F"))

ITensors.op(::OpName"AupF⋅", ::SiteType"vElectron") = premul(elop("Aup") * elop("F"))
ITensors.op(::OpName"⋅AupF", ::SiteType"vElectron") = postmul(elop("Aup") * elop("F"))

ITensors.op(::OpName"Adn⋅", ::SiteType"vElectron") = premul(elop("Adn"))
ITensors.op(::OpName"⋅Adn", ::SiteType"vElectron") = postmul(elop("Adn"))

ITensors.op(::OpName"Adn†⋅", ::SiteType"vElectron") = premul(elop("Adagdn"))
ITensors.op(::OpName"⋅Adn†", ::SiteType"vElectron") = postmul(elop("Adagdn"))

ITensors.op(::OpName"FAdn⋅", ::SiteType"vElectron") =  premul(elop("F") * elop("Adn"))
ITensors.op(::OpName"⋅FAdn", ::SiteType"vElectron") = postmul(elop("F") * elop("Adn"))

ITensors.op(::OpName"FAdn†⋅", ::SiteType"vElectron") =  premul(elop("F") * elop("Adagdn"))
ITensors.op(::OpName"⋅FAdn†", ::SiteType"vElectron") = postmul(elop("F") * elop("Adagdn"))

ITensors.op(::OpName"Adn†F⋅", ::SiteType"vElectron") =  premul(elop("Adagdn") * elop("F"))
ITensors.op(::OpName"⋅Adn†F", ::SiteType"vElectron") = postmul(elop("Adagdn") * elop("F"))
