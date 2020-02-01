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
        h1::Array{Float64,2}  
        h2::Array{Float64,4}  
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

function compute_ab_direct(H)
    #={{{=#

    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=H.no, ne=H.na)
    ket_b = ConfigStrings.ConfigString(no=H.no, ne=H.nb)
    bra_a = ConfigStrings.ConfigString(no=H.no, ne=H.na)
    bra_b = ConfigStrings.ConfigString(no=H.no, ne=H.nb)

    ket_a_ca_lookup = ConfigStrings.fill_ca_lookup(ket_a)
    ket_b_ca_lookup = ConfigStrings.fill_ca_lookup(ket_b)

    ConfigStrings.reset!(ket_b)
                           
    Hout = zeros(ket_a.max*ket_b.max,ket_a.max*ket_b.max)
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
                  
                    Lb = 1
                    sign_b = 1
                    L = 1 
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_ca_lookup[Kb][q+(s-1)*ket_b.no]
                            
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * bra_a.max

                            #sig[K,:] =  H.h2[p,r,q,s]
                            Hout[K,L] += H.h2[p,r,q,s] * sign_a * sign_b
                            #scr = sig[:,K] + H.h2[p,r,q,s] * sign_a * sign_b * v[:,L]
                            #sig[:,K] = scr 

                            #for si in 1:n_roots
                            #    sig[K,si] += Iprqs * sign_a * sign_b * v[L,si]
                            #end
                            continue
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end
    return Hout 
end
#=}}}=#

function compute_ab_terms(vin::Array{Float64,2}, H)
    #={{{=#

    v = transpose(vin)
    sig = 0*v
    if ndims(v)>1
        n_roots = size(v,2)
    else
        n_roots = 1
    end
    
    @assert(ndims(v) == ndims(sig))
    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=H.no, ne=H.na)
    ket_b = ConfigStrings.ConfigString(no=H.no, ne=H.nb)
    bra_a = ConfigStrings.ConfigString(no=H.no, ne=H.na)
    bra_b = ConfigStrings.ConfigString(no=H.no, ne=H.nb)

    ket_a_ca_lookup = ConfigStrings.fill_ca_lookup(ket_a)
    ket_b_ca_lookup = ConfigStrings.fill_ca_lookup(ket_b)

    ConfigStrings.reset!(ket_b)
                           
    scr = zeros(1,ket_a.max*ket_b.max)
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
                  
                    Lb = 1
                    sign_b = 1
                    L = 1 
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_ca_lookup[Kb][q+(s-1)*ket_b.no]
                            
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * bra_a.max

                            #sig[K,:] =  H.h2[p,r,q,s]
                            sig[:,K] += H.h2[p,r,q,s] * sign_a * sign_b * v[:,L]
                            #scr = sig[:,K] + H.h2[p,r,q,s] * sign_a * sign_b * v[:,L]
                            #sig[:,K] = scr 

                            #for si in 1:n_roots
                            #    sig[K,si] += Iprqs * sign_a * sign_b * v[L,si]
                            #end
                            continue
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end
    sig = transpose(sig)
    return sig 
end
#=}}}=#

function precompute_spin_diag_terms(H, e)
    #={{{=#

    #   Create local references to ci_strings
    ket = ConfigStrings.ConfigString(no=H.no, ne=e)
    bra = ConfigStrings.ConfigString(no=H.no, ne=e)

    ket_ca_lookup = ConfigStrings.fill_ca_lookup(ket)
    
    Hout = zeros(ket.max,ket.max)

    ConfigStrings.reset!(ket)

    for K in 1:ket.max

        #  hpq p'q 
        for p in 1:ket.no
            for q in 1:ket.no
                bra = deepcopy(ket)
                ConfigStrings.apply_annihilation!(bra,q)
                if bra.sign == 0
                    continue
                end
                ConfigStrings.apply_creation!(bra,p)
                if bra.sign == 0
                    continue
                end

                L = ConfigStrings.calc_linear_index(bra)

                term = H.h1[q,p]
                Hout[K,L] += term * bra.sign
            end
        end


        #  <pq|rs> p'q'sr -> (pr|qs) 
        for r in 1:ket.no
            for s in r+1:ket.no
                for p in 1:ket.no
                    for q in p+1:ket.no

                        bra = deepcopy(ket)

                        ConfigStrings.apply_annihilation!(bra,r) 
                        if bra.sign == 0
                            continue
                        end
                        ConfigStrings.apply_annihilation!(bra,s) 
                        if bra.sign == 0
                            continue
                        end
                        ConfigStrings.apply_creation!(bra,q) 
                        if bra.sign == 0
                            continue
                        end
                        ConfigStrings.apply_creation!(bra,p) 
                        if bra.sign == 0
                            continue
                        end
                        L = ConfigStrings.calc_linear_index(bra)
                        Ipqrs = H.h2[p,r,q,s]-H.h2[p,s,q,r]
                        Hout[K,L] += Ipqrs*bra.sign
                    end
                end
            end
        end
        ConfigStrings.incr!(ket)
    end
    return Hout
end
#=}}}=#

function matvec(v,H)
    #=  
    Function to compute the action of the Hamiltonian onto a vector
    =#
    sigma = compute_ab_terms(v,H) 
end


end
