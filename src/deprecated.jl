function ITensors.op(::OpName"Damping", ::SiteType"vOsc", d::Int; kwargs...)
    Base.depwarn(
        "`ITensors.op(::OpName\"Damping\", ::SiteType\"vOsc\", d::Int; kwargs...)` is
        deprecated; use `dissipator(n::Int, frequency::Real, temperature::Real)` instead,
        where `n` is the site where the dissipator should be applied. The newer function
        will return an OpSum object."
    )
    A = ITensors.op(OpName("A"), SiteType("Osc"), d)
    Adag = ITensors.op(OpName("Adag"), SiteType("Osc"), d)
    if temperature == 0
        fn = x -> A * x * Adag - 0.5 * (Adag * A * x + x * Adag * A)
    else
        avgn = 1 / expm(frequency / temperature)
        fn =
            x ->
                (avgn + 1) * (A * x * Adag - 0.5 * (Adag * A * x + x * Adag * A)) +
                avgn * (Adag * x * A - 0.5 * (A * Adag * x + x * A * Adag))
    end
    return vec(fn, gellmannbasis(d))
end

function ITensors.op(::OpName"Lindb+", ::SiteType"vOsc", d::Int)
    Base.depwarn(
        "`ITensors.op(::OpName\"Lindb+\", st::SiteType\"vOsc\", d::Int)` is deprecated;
        use `dissipator_loss(n::Int)` instead, where `n` is the site where the dissipator
        should be applied. The newer function will return an OpSum object."
    )
    A = ITensors.op(OpName("A"), SiteType("Osc"), d)
    Adag = ITensors.op(OpName("Adag"), SiteType("Osc"), d)
    fn = x -> A * x * Adag - 0.5 * (Adag * A * x + x * Adag * A)
    return vec(fn, gellmannbasis(d))
end
function ITensors.op(::OpName"Lindb-", ::SiteType"vOsc", d::Int)
    Base.depwarn(
        "`ITensors.op(::OpName\"Lindb-\", st::SiteType\"vOsc\", d::Int)` is deprecated;
        use `dissipator_gain(n::Int)` instead, where `n` is the site where the dissipator
        should be applied. The newer function will return an OpSum object."
    )
    A = ITensors.op(OpName("A"), SiteType("Osc"), d)
    Adag = ITensors.op(OpName("Adag"), SiteType("Osc"), d)
    fn = x -> Adag * x * A - 0.5 * (A * Adag * x + x * A * Adag)
    return vec(fn, gellmannbasis(d))
end

function mixedlindbladplus(s1::Index{Int64}, s2::Index{Int64})
    Base.depwarn("`mixedlindbladplus(s1::Index{Int64}, s2::Index{Int64}`) is deprecated;
                 use `mixedlindbladplus(n1::Int, n2::Int)` instead, where `n1` and `n2`
                 are the positions of `s1` and `s2` in the system. The newer function
                 will return an OpSum object.")
    return (
        op("a-⋅", s1) * op("⋅a+", s2) + op("a-⋅", s2) * op("⋅a+", s1) -
        0.5 * (op("a+⋅", s1) * op("a-⋅", s2) + op("a+⋅", s2) * op("a-⋅", s1)) -
        0.5 * (op("⋅a+", s1) * op("⋅a-", s2) + op("⋅a+", s2) * op("⋅a-", s1))
    )
end
function mixedlindbladminus(s1::Index{Int64}, s2::Index{Int64})
    Base.depwarn("`mixedlindbladminus(s1::Index{Int64}, s2::Index{Int64}`) is deprecated;
                 use `mixedlindbladminus(n1::Int, n2::Int)` instead, where `n1` and `n2`
                 are the positions of `s1` and `s2` in the system. The newer function
                 will return an OpSum object.")
    return (
        op("a+⋅", s1) * op("⋅a-", s2) + op("a+⋅", s2) * op("⋅a-", s1) -
        0.5 * (op("a-⋅", s1) * op("a+⋅", s2) + op("a-⋅", s2) * op("a+⋅", s1)) -
        0.5 * (op("⋅a-", s1) * op("⋅a+", s2) + op("⋅a-", s2) * op("⋅a+", s1))
    )
end
