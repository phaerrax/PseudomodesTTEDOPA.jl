# Application: GKSL equation

With the states and operators defined by this package, we can define a Lindbladian
operator acting on density matrices in a GKSL equation.
Let's take, as an example, a system of ``N=3`` harmonic oscillators described by the
Hamiltonian

```math
H=
\sum_{j=1}^{N}\omega_j a^\dagger_j a_j + 
\sum_{j=1}^{N-1}\lambda_j (a^\dagger_j a_{j+1} + a^\dagger_{j+1} a_j)
```

and where each oscillator is damped according to the dissipator operator

```math
\mathcal{D}(\rho)=
\sum_{j=1}^{N} \gamma_j \bigl( a_j \rho a^\dagger_j - \tfrac12 \rho a^\dagger_j a_j -
\tfrac12 a^\dagger_j a_j \rho \bigr).
```

First, create the site indices. In this example, we will use an arbitrary value for the
number of levels of the oscillators:

```julia
N = 3
s = siteinds("vOsc", 3; dim=5)
```

Create the OpSum object representing the unitary part of the evolution, ``-i[H,{-}]``:

```julia
ℓ = OpSum()
for j in 1:N
    ℓ += ω[j] * gkslcommutator("N", j)
end
for j in 1:(N - 1)
    ℓ +=
        λ[j] * (
            gkslcommutator("Adag", site1, "A", site2) +
            gkslcommutator("A", site1, "Adag", site2)
        )
end
```

Finally, add the dissipator:

```julia
for j in 1:N
    ℓ += γ, "A⋅ * ⋅Adag", n
    ℓ += -γ / 2, "N⋅", n
    ℓ += -γ / 2, "⋅N", n
end
```

This OpSum object can now, for example, be cast into an MPO and used to evolve a state with
the TDVP algorithm. Here is how ``L`` can be applied to a state ``ρ``:

```julia
ρ = MPS([state("ThermEq", s[j]; temperature=1.0, frequency=5.0) for j in 1:N])
L = MPO(ℓ, s)
apply(L, ρ)
```

Computing the expectation value of an operator ``A`` on a state ``ρ`` means
computing ``\operatorname{tr}(Aρ)``: with this package, this is done by contracting the
ITensor state representing the adjoint of the vectorized operator with the ITensor state
representing ``ρ``.
For example, if ``A`` is the number operator on site 3, i.e ``I\otimes I\otimes N``,
then we need the following code.

```julia
vec_n = MPS(s, [i == 3 ? "vN" : "vId" for i in 1:N])
result = dot(vec_n, ρ)
```

As a particular case, the trace of the state can be obtained by "measuring" the identity
operator:

```julia
vec_id = MPS(s, "vId")
trace = dot(vec_id, ρ)
```
