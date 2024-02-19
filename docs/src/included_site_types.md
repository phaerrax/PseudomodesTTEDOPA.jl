# SiteTypes included with PseudomodesTTEDOPA

This package extends the SiteType collection provided by ITensor, by defining types for
(vectorized) density matrices associated to built-in and operators acting on them.

For each of these site types, each operator from the "original", non-vectorized ITensor site
type is promoted as a pre- or post-multiplication operator.
In other words, if there exists an operator "A" for the "Fermion" type, then you will
find operators "A⋅" and "⋅A" (the dot here is a `\cdot`) available for the "vFermion"
site type.

Moreover, the `gkslcommutator` function allows you to easily define the unitary term (the
commutator) in the GKSL equation: given a list of operator names and site indices,
it returns an OpSum object representing the ``x\mapsto -i[A,x]`` map, for example
```julia-repl
julia> gkslcommutator("A", 1)
sum(
  0.0 - 1.0im A⋅(1,)
  0.0 + 1.0im ⋅A(1,)
)

julia> gkslcommutator("A", 1, "B", 3)
sum(
  0.0 - 1.0im A⋅(1,) B⋅(3,)
  0.0 + 1.0im ⋅A(1,) ⋅B(3,)
)
```

Operators are also vectorized, with the aim of computing expectation values: if
``\{a_i\}_{i=1}^{d}`` and ``\{r_i\}_{i=1}^{d}`` are, respectively, the coordinates of
the operator ``A`` and the state ``\rho`` with respect to the Gell-Mann basis, then the
expectation value ``\operatorname{tr}(\rho A)`` is given by ``\sum_{i=1}^{d}\overline{a_i}r_i``.
Below you can find some operators transformed this way: they are actually "states" in the
ITensor formalism, even if they represent operators.
Their operator names are usually obtained by prefixing a `v` to the original ITensor name.

## "vS=1/2" SiteType

Site indices with the "vS=1/2" site type represent the density matrix of a ``S=1/2`` spin,
namely its coefficients with respect to the basis of 2x2 Gell-Mann matrices.

Making a single "vS=1/2" site or collection of N "vS=1/2" sites
```
s = siteind("vS=1/2")
sites = siteinds("vS=1/2", N)
```


#### "vS=1/2" states

The available state names for "vS=1/2" sites are:
- `"Up"` spin in the ``|{\uparrow}\rangle\langle{\uparrow}|`` state
- `"Dn"` spin in the ``|{\downarrow}\rangle\langle{\downarrow}|`` state

#### "vS=1/2" vectorized operators

Vectorized operators or associated with "vS=1/2" sites can be made using the `state`
function, for example
```
Sx4 = state("vSx", s, 4)
N3 = state("vN", sites[3])
```

Spin operators:
- `"vId"` Identity operator ``I_2``
- `"vσx"` Pauli x matrix ``\sigma^x``
- `"vσy"` Pauli y matrix ``\sigma^y``
- `"vσz"` Pauli z matrix ``\sigma^z``
- `"vSx"` Spin x operator ``S^x = \frac{1}{2} \sigma_x``
- `"vSy"` Spin y operator ``S^y = \frac{1}{2} \sigma_y``
- `"vSz"` Spin z operator ``S^z = \frac{1}{2} \sigma_z``
- `"vN"` Number operator ``N = \frac{1}{2} (I_2+\sigma_z)``


## "Osc" SiteType

A SiteType representing a harmonic oscillator. It inherits definitions from the "Qudit"
type.

Available keyword arguments for customization:
- `dim` (default: 2): dimension of the index (number of oscillator levels)
For example:
```
sites = siteinds("Osc", N; dm=3)
```

#### "vOsc" Operators

Operators associated with "Osc" sites:

Single-oscillator operators:
- `"A"` (aliases: `"a"`) annihilation operator
- `"Adag"` (aliases: `"adag"`, `"a†"`) creation operator
- `"Asum"` equal to ``a+a^\dagger``
- `"X"` equal to ``\frac{1}{\sqrt{2}}(a^\dagger+a)``
- `"Y"` equal to ``\frac{i}{\sqrt{2}}(a^\dagger-a)``
- `"N"` (aliases: `"n"`) number operator


## "vOsc" SiteType

Available keyword arguments for customization:
- `dim` (default: 2): dimension of the index (number of oscillator levels)
For example:
```
sites = siteinds("vOsc", N; dim=4)
```

#### "vOsc" states

States associated with "vOsc" sites.
- "`n`", where `n` is an integer between `0` and `dim-1`, gives the Fock state
  ``|n\rangle\langle n|``
- `"ThermEq"`, with additional parameters `temperature` and `frequency`, gives the thermal
  equilibrium (Gibbs) state ``\operatorname{tr}(\exp(-\frac{\omega}{T}
N))\exp(-\frac{\omega}{T} N)``
- `"X⋅ThermEq"`, with additional parameters `temperature` and `frequency`, gives the thermal
  equilibrium state (as in the previous point) multiplied by ``X=\frac{1}{\sqrt{2}}(a^\dagger+a)`` on the left (this is not actually a state, but it is useful when computing the correlation function of the ``X`` operator)

Example:
```julia
s = siteind("vOsc", N; dim=4)
rho_eq = state("ThermEq", s; temperature=1.0, frequency=5.0)
fock_st = state("3", s)
```


#### "vOsc" Operators

Vectorized operators associated with "vOsc" sites.

Single-oscillator operators (see the operator for "Osc" sites for their meaning):
- `"vA"`
- `"vAdag"`
- `"vN"`
- `"vX"`
- `"vY"`
- `"vId"`


## "vFermion" SiteType

Site indices with the "vFermion" SiteType represent spinless fermion sites with the states
``|0\rangle``, ``|1\rangle``, corresponding to zero fermions or one fermion.

#### "vFermion" states

The available state names for "vFermion" sites are:
- `"Emp"` unoccupied fermion site
- `"Occ"` occupied fermion site

#### "vFermion" operators

Vectorized operators associated with "vFermion" sites:

- `"vId"` Identity operator
- `"vN"` Density operator
- `"vA"` (aliases: `"va"`) Fermion annihilation operator
- `"vAdag"` (aliases: `"vadag"`, `"vA†"`, `"va†"`) Fermion creation operator

Note that these creation and annihilation operators do not include Jordan-Wigner strings.


## "vElectron" SiteType

The states of site indices with the "vElectron" SiteType correspond to
``|0\rangle``, ``|{\uparrow}\rangle``, ``|{\downarrow}\rangle``, ``|{\uparrow\downarrow}\rangle``.

#### "vElectron" states

The available state names for "vElectron" sites are:
- `"Emp"` unoccupied electron site
- `"Up"` electron site occupied with one up electron
- `"Dn"` electron site occupied with one down electron
- `"UpDn"` electron site occupied with two electrons (one up, one down)

#### "vElectron" operators

Vectorized operators associated with "vElectron" sites:

- `"vId"` Identity operator
- `"vNtot"` Total density operator
- `"vNup"` Up density operator
- `"vNdn"` Down density operator
- `"vNupNdn"` Product of `n↑` and `n↓`

