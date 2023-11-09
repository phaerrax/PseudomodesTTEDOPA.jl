# Additional ITensor states and operators for the electron SiteType.
function ITensors.op(::OpName"Id", ::SiteType"Electron") # ITensors doesn't define this
    return Matrix(1.0I, 4, 4)
end

ITensors.op(::OpName"Ntot^2", st::SiteType"Electron") = ITensors.op(OpName("Ntot"), st)^2
ITensors.op(::OpName"ntot^2", st::SiteType"Electron") = ITensors.op(OpName("Ntot^2"), st)

# These aliases aren't defined in ITensors, but we use them a lot in our scripts,
# so we keep them for compatibility reasons.
ITensors.op(::OpName"NupNdn", st::SiteType"Electron") = ITensors.op(OpName("Nupdn"), st)
ITensors.op(::OpName"Aup†", st::SiteType"Electron") = ITensors.op(OpName("Adagup"), st)
ITensors.op(::OpName"Adn†", st::SiteType"Electron") = ITensors.op(OpName("Adagdn"), st)

# Here are some compositions of operators that are commonly used in defining exchange
# interactions with Electron site types: we need them, actually, for `vElectron` types,
# where multiplication between operators is trickier to use.
function ITensors.op(::OpName"Aup†F", st::SiteType"Electron")
    return ITensors.op(OpName("Adagup"), st) * ITensors.op(OpName("F"), st)
end
ITensors.op(::OpName"a†↑F", st::SiteType"Electron") = ITensors.op(OpName("Aup†F"), st)

function ITensors.op(::OpName"AupF", st::SiteType"Electron")
    return ITensors.op(OpName("Aup"), st) * ITensors.op(OpName("F"), st)
end
ITensors.op(::OpName"a↑F", st::SiteType"Electron") = ITensors.op(OpName("AupF"), st)

function ITensors.op(::OpName"FAdn", st::SiteType"Electron")
    return ITensors.op(OpName("F"), st) * ITensors.op(OpName("Adn"), st)
end
ITensors.op(::OpName"Fa↓", st::SiteType"Electron") = ITensors.op(OpName("FAdn"), st)

function ITensors.op(::OpName"FAdn†", st::SiteType"Electron")
    return ITensors.op(OpName("F"), st) * ITensors.op(OpName("Adagdn"), st)
end
ITensors.op(::OpName"Fa†↓", st::SiteType"Electron") = ITensors.op(OpName("FAdn†"), st)

function ITensors.op(::OpName"Adn†F", st::SiteType"Electron")
    return ITensors.op(OpName("Adagdn"), st) * ITensors.op(OpName("F"), st)
end
ITensors.op(::OpName"a†↓F", st::SiteType"Electron") = ITensors.op(OpName("Adn†F"), st)
