# Additional ITensor states and operators for the electron SiteType.

ITensors.op(::OpName"NupNdn", st::SiteType"Electron") = ITensors.op(OpName("Nup"), st) * ITensors.op(OpName("Ndn"), st)
