push!(LOAD_PATH, "./")
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


@with_kw struct Hamiltonian 
    #=
    Type to organize all the CI data 

        h0  is constant energy shift
        h1  is one body operator
        h2  is two body operator
    =#
# {{{
	h0::Real 
        h1::Array{Real,2}  
        h2::Array{Real,4}  
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

function compute_ab_terms!(v, sig, H)
    @assert(ndims(v) == ndims(sig))
    n_roots = size(v,2) 
    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=H.no, ne=H.na)
    ket_b = ConfigStrings.ConfigString(no=H.no, ne=H.nb)
    bra_a = ConfigStrings.ConfigString(no=H.no, ne=H.na)
    bra_b = ConfigStrings.ConfigString(no=H.no, ne=H.nb)

    ket_a_ca_lookup = ConfigStrings.fill_ca_lookup(ket_a)
    ket_b_ca_lookup = ConfigStrings.fill_ca_lookup(ket_b)

    ConfigStrings.reset!(ket_b)

    for Kb in 1:ket_b.max

        ConfigStrings.reset!(ket_a)
        for Ka in 1:ket_a.max
            K = Ka + (Kb-1) * ket_a.max

            #  <pq|rs> p'q'sr  --> (pr|qs) (a,b)
            for r in 1:ket_a.no
                for p in 1:ket_a.no
                    sign_a, La = ket_a_ca_lookup[Ka][p+(r-1)*ket_a.no]
                    if La == 0
                        continue
                    end
                    
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_ca_lookup[Kb][q+(s-1)*ket_b.no]
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * bra_a.max


                            Iprqs = H.h2[p,r,q,s]

                            for si in 1:n_roots
                                sig[K,si] += Iprqs * sign_a * sign_b * v[L,si]
                            end
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end

    return sig 
end

function matvec(v, H=H)
    #=  
    Function to compute the action of the Hamiltonian onto a vector
    =#
    sigma = compute_ab_terms(v,H) 
end


end
