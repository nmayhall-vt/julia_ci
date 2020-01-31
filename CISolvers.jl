module CISolvers
using Printf
using Parameters

import ConfigStrings
import Tools 

using IterativeSolvers, LinearMaps


@with_kw mutable struct State 
    #=
    Type to organize all the information related to optimizing 
    a dense vector over a space of determinants (ConfigString types)
    =#
# {{{
	no::Int = 0
	na::Int = 0  # number of alpha
	nb::Int = 0  # number of beta
        dima::Int = Tools.calc_nchk(no,na) 
        dimb::Int = Tools.calc_nchk(no,nb)
        dim::Int = dima*dimb
        converged::Bool = false
        restarted::Bool = false
        iteration::Int = 0
        algorithm::String = "direct"    #  options: direct/davidson
        n_roots::Int = 1
        vec_curr::Array{Number,1} = []
end
## }}}

end
