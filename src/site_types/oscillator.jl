# Space of harmonic oscillators
# =============================

ITensors.alias(::SiteType"Osc") = SiteType"Qudit"()

"""
    ITensors.space(st::SiteType"Osc";
                   dim = 2,
                   conserve_qns = false,
                   conserve_number = false,
                   qnname_number = "Number")

Create the Hilbert space for a site of type "Osc".

Optionally specify the conserved symmetries and their quantum number labels.
"""
ITensors.space(st::SiteType"Osc"; kwargs...) = space(ITensors.alias(st); kwargs...)

# Forward "Osc" definitions to "Qudit" states.
# This way each state, operator and so on defined for the "Qudit" SiteType is already
# available for the "Osc" type.
ITensors.val(vn::ValName; st::SiteType"Osc") = val(vn, ITensors.alias(st))

function ITensors.state(sn::StateName, st::SiteType"Osc", s::Index; kwargs...)
    return state(sn, ITensors.alias(st), s; kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"Osc", dims::Int...; kwargs...)
    return op(on, ITensors.alias(st), dims...; kwargs...)
end

# ITensors.op functions all require an additional parameter, the dimension: for this reason,
# we need to compute ITensors.dim.(rs), i.e. the dimensions of the Indices of the operator,
# and append them to the function arguments.
function ITensors.op(on::OpName, st::SiteType"Osc", s1::Index, s_tail::Index...; kwargs...)
    rs = reverse((s1, s_tail...))
    dims = ITensors.dim.(rs)
    opmat = ITensors.op(on, st, dims...; kwargs...)
    return ITensors.itensor(opmat, prime.(rs)..., dag.(rs)...)
end

# Operators
# ---------
function ITensors.op(::OpName"Asum", st::SiteType"Osc", d::Int)
    return ITensors.op(OpName("A"), st, d) + ITensors.op(OpName("Adag"), st, d)
end
function ITensors.op(on::OpName"X", st::SiteType"Osc", d::Int)
    return ITensors.op(OpName("Asum"), st, d) / sqrt(2)
end
function ITensors.op(on::OpName"Y", st::SiteType"Osc", d::Int)
    return im / sqrt(2) *
           (ITensors.op(OpName("Adag"), st, d) - ITensors.op(OpName("A"), st, d))
end

# Aliases (for backwards compatibility)
ITensors.alias(::OpName"a-") = OpName"A"()
ITensors.alias(::OpName"a+") = OpName"Adag"()
ITensors.alias(::OpName"asum") = OpName"Asum"()

function ITensors.op(on::OpName"a+", st::SiteType"Osc", d::Int)
    return ITensors.op(ITensors.alias(on), st, d)
end
function ITensors.op(on::OpName"a-", st::SiteType"Osc", d::Int)
    return ITensors.op(ITensors.alias(on), st, d)
end
function ITensors.op(on::OpName"asum", st::SiteType"Osc", d::Int)
    return ITensors.op(ITensors.alias(on), st, d)
end
