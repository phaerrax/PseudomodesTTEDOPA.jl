module PseudomodesTTEDOPA

using Base.Filesystem
using CSV
using Combinatorics
using DataFrames
using ITensors
using JSON
using LinearAlgebra
using Measures
using Observers
using Plots
using PolyChaos
using ProgressMeter
using QuadGK

import ITensors.state,
       ITensors.op,
       ITensors.space

export allequal,
       #canonicalbasis,
       #canonicalmatrix,
       chain,
       chop,
       readablecolours,
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
       partialtrace,
       #sitetypes,
       vec,
       vonneumannentropy

include("utils.jl")

export groupplot,
       unifiedlogplot,
       unifiedplot

include("deprecated.jl")

export chain_L1_state,
       chain_basis_states,
       level_subspace_proj,
       parse_init_state,
       parse_spin_state,
       single_ex_state,
       spin_current_op_list

include("spin_chain_space.jl")

export mixedlindbladminus,
       mixedlindbladplus,
       osc_levels_proj,
       oscdimensions,
       parse_init_state_osc

include("harmonic_oscillator_space.jl")

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
       hamiltonian_xy

include("operators.jl")

export getchaincoefficients,
       thermalisedJ,
       chainmapcoefficients

include("tedopa.jl")

end
