# Additional ITensor states and operators for the electron SiteType.

function ITensors.op(::OpName"NupNdn", st::SiteType"Electron")
    return ITensors.op(OpName("Nup"), st) * ITensors.op(OpName("Ndn"), st)
end
