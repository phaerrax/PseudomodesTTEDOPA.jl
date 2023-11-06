# Additional ITensor states and operators for the electron SiteType.

function ITensors.op(::OpName"NupNdn", st::SiteType"Electron")
    return ITensors.op(OpName("Nup"), st) * ITensors.op(OpName("Ndn"), st)
end

ITensors.op(::OpName"Ntot^2", st::SiteType"Electron") = ITensors.op(OpName("Ntot"), st)^2
ITensors.op(::OpName"ntot^2", st::SiteType"Electron") = ITensors.op(OpName("Ntot^2"), st)
