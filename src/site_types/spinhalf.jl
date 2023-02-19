# Additional ITensor states and operators for the spin-half SiteType.
#
# Note that while the following definitions use "S=1/2" to denote the site type, they hold
# for all other aliases of the spin-half type, such as "S=½", automatically.

# Null vector
ITensors.state(::StateName"0", ::SiteType"S=1/2") = [0; 0]
# Number operator, aka |↑⟩⟨↑|
ITensors.op(::OpName"N", st::SiteType"S=1/2") = ITensors.op(OpName("ProjUp"), st)
# Identity matrices
ITensors.op(::OpName"Id", ::SiteType"S=1/2") = Matrix(1.0I, 2, 2)
# Null matrix
ITensors.op(::OpName"0", ::SiteType"S=1/2") = zeros(2, 2)

# Ladder operators (aliases)
ITensors.op(::OpName"σ+", st::SiteType"S=1/2") = op(OpName("S+"), st)
ITensors.op(::OpName"σ-", st::SiteType"S=1/2") = op(OpName("S-"), st)
ITensors.op(::OpName"plus", st::SiteType"S=1/2") = op(OpName("S+"), st)
ITensors.op(::OpName"minus", st::SiteType"S=1/2") = op(OpName("S-"), st)
