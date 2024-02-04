# PseudomodesTTEDOPA.jl

Documentation for PseudomodesTTEDOPA.jl

### Spin chain operators
```@docs
exchange_interaction
exchange_interaction_adjoint
exchange_interaction′
spin_chain
spin_chain_adjoint
spin_chain′
```

### Pseudomode operators
```@docs
dissipator_loss
dissipator_gain
dissipator
mixedlindbladplus
mixedlindbladminus
```

### Chain-mapping
```@docs
getchaincoefficients
chainmapcoefficients
```

### Markovian closure operators
```@docs
Closure
```
```@docs
closure
```
```@docs
length(mc::Closure)
freqs(mc::Closure)
innercoups(mc::Closure)
outercoups(mc::Closure)
damps(mc::Closure)
freq(mc::Closure, j::Int)
innercoup(mc::Closure, j::Int)
outercoup(mc::Closure, j::Int)
damp(mc::Closure, j::Int)
```

```@docs
closure_op(mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int)
closure_op_adjoint(
    mc::Closure, sites::Vector{<:Index}, chain_edge_site::Int, gradefactor::Int
)
closure_op′
filled_closure_op
filled_closure_op_adjoint
filled_closure_op′
```
