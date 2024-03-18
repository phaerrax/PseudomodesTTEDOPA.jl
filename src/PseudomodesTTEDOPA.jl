module PseudomodesTTEDOPA

using Base.Filesystem
using CSV
using Combinatorics
using DataFrames
using ITensors
using IterTools
using JSON
using LinearAlgebra
using Measures
using Observers
using PolyChaos
using ProgressMeter
using QuadGK

import ITensors.state, ITensors.op, ITensors.space

export allequal,
    #canonicalbasis,
    #canonicalmatrix,
    chain,
    chop,
    consecutivepairs,
    construct_step_list,
    disablegrifqtech,
    embed_slice,
    filenamett,
    gellmannbasis,
    gellmannmatrix,
    groupresults,
    #isjson,
    levels,
    #linkdims,
    load_parameters,
    oscdimensions,
    parse_init_state_osc,
    partialtrace,
    #sitetypes,
    vec,
    vonneumannentropy,
    chain_L1_state,
    chain_basis_states,
    level_subspace_proj,
    parse_init_state,
    parse_spin_state,
    single_ex_state,
    sitenumber,
    jwstring

include("utils.jl")

export spin_current_op_list

include("current_operators.jl")

include("site_types/spinhalf.jl")
include("site_types/vectorized_spinhalf.jl")
include("site_types/fermion.jl")
include("site_types/vectorized_fermion.jl")
include("site_types/oscillator.jl")
include("site_types/electron.jl")
include("site_types/vectorized_electron.jl")
include("site_types/fermionic_dot3.jl")
include("site_types/vectorized_fermionic_dot3.jl")

export dissipator_loss, dissipator_gain, dissipator, mixedlindbladplus, mixedlindbladminus

include("site_types/vectorized_oscillator.jl")

export twositeoperators,
    localop,
    interactionop,
    evolve,
    fermioncurrent,
    forwardflux,
    backwardflux,
    dissipator_symmetric,
    dissipator_asymmetric,
    lindbladian_xy,
    hamiltonian_xy,
    gkslcommutator,
    gkslcommutator_itensor

include("operators.jl")

include("spin_chain.jl")
include("exchange_interaction.jl")
include("closure.jl")

export getchaincoefficients, thermalisedJ, chainmapcoefficients

include("tedopa.jl")

export mixedlindbladminus, mixedlindbladplus

include("deprecated.jl")

@warn """
This package was deprecated on 18 March 2024.
Its features have been split in three different packages:

* TEDOPA (https://github.com/phaerrax/TEDOPA.jl) for the chain-mapping
  functionality;
* LindbladVectorizedTensors (https://github.com/phaerrax/LindbladVectorizedTensors.jl)
  for the ITensors site type definitions and functions;
* MarkovianClosure (https://github.com/phaerrax/MarkovianClosure.jl) for the
  Markovian closure constructions.

This package will remain available, but it won't receive updates anymore.
"""

end
