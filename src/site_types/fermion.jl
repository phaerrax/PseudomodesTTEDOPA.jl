function ITensors.op(::OpName"Id", ::SiteType"Fermion")
    return Matrix(1.0I, 2, 2)
end

function ITensors.op!(Op::ITensor, ::OpName"A", ::SiteType"Fermion", s::Index)
  return Op[s' => 1, s => 2] = 1.0
end
function ITensors.op!(Op::ITensor, on::OpName"a", st::SiteType"Fermion", s::Index)
    return ITensors.op!(Op, OpName("A"), st, s)
end

function ITensors.op!(Op::ITensor, ::OpName"Adag", ::SiteType"Fermion", s::Index)
  return Op[s' => 2, s => 1] = 1.0
end
function ITensors.op!(Op::ITensor, on::OpName"adag", st::SiteType"Fermion", s::Index)
    return ITensors.op!(Op, OpName("Adag"), st, s)
end
function ITensors.op!(Op::ITensor, on::OpName"a†", st::SiteType"Fermion", s::Index)
    return ITensors.op!(Op, OpName("Adag"), st, s)
end
