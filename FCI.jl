push!(LOAD_PATH, "./")

module FCI
    

using LinearAlgebra 
using Printf
using Parameters

import ConfigStrings
import Tools 

using LinearMaps

struct Hamiltonian 
    #=
    Type to hold a second quantized Hamiltonian coefficients 

        h0  is constant energy shift
        h1  is one body operator
        h2  is two body operator
    =#
# {{{
	h0::Real 
        h1::Array{Float64,2}  
        h2::Array{Float64,4}  
end
## }}}

@with_kw struct Problem
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
end


function compute_spin_diag_terms_full!(H::Hamiltonian, P::Problem, Hmat)
    #={{{=#

    print(" Compute same spin terms.\n")
    @assert(size(Hmat,1) == P.dim)
    Hdiag_a = FCI.precompute_spin_diag_terms(H,P,P.na)
    Hdiag_b = FCI.precompute_spin_diag_terms(H,P,P.nb)
    Hmat .+= kron(Matrix(1.0I, P.dimb, P.dimb), Hdiag_a)
    Hmat .+= kron(Hdiag_b, Matrix(1.0I, P.dima, P.dima))
    
end
#=}}}=#


function compute_ab_terms_full!(H::Hamiltonian, P::Problem, Hmat)
    #={{{=#

    print(" Compute opposite spin terms.\n")
    @assert(size(Hmat,1) == P.dim)
    
    #v = transpose(vin)
    
    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    ket_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)
    bra_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    bra_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)

    
    ket_a_lookup = ConfigStrings.fill_ca_lookup(ket_a)
    ket_b_lookup = ConfigStrings.fill_ca_lookup(ket_b)
    
    ConfigStrings.reset!(ket_b)
                          
    for Kb in 1:ket_b.max

        ConfigStrings.reset!(ket_a)
        for Ka in 1:ket_a.max
            K = Ka + (Kb-1) * ket_a.max

            #  <pq|rs> p'q'sr  --> (pr|qs) (a,b)
            for r in 1:ket_a.no
                for p in 1:ket_a.no
                    sign_a, La = ket_a_lookup[Ka][p+(r-1)*ket_a.no]
                    if La == 0
                        continue
                    end
                  
                    Lb = 1
                    sign_b = 1
                    L = 1 
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_lookup[Kb][q+(s-1)*ket_b.no]
                            
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * bra_a.max
                            
                            Hmat[K,L] += H.h2[p,r,q,s] * sign_a * sign_b
                           
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end
    return  
end
#=}}}=#


function compute_ab_terms_full(H::Hamiltonian, P::Problem)
    #={{{=#

    #print(" Compute opposite spin terms. Shape of v: ", size(v), "\n")
    
    #v = transpose(vin)
    
    Hmat = zeros(Float64, P.dim, P.dim)
 
    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    ket_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)
    bra_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    bra_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)

    ket_a_lookup = ConfigStrings.fill_ca_lookup(ket_a)
    ket_b_lookup = ConfigStrings.fill_ca_lookup(ket_b)
    
    a_max = bra_a.max
    ConfigStrings.reset!(ket_b)
                          
    for Kb in 1:ket_b.max

        ConfigStrings.reset!(ket_a)
        for Ka in 1:ket_a.max
            K = Ka + (Kb-1) * ket_a.max

            #  <pq|rs> p'q'sr  --> (pr|qs) (a,b)
            for r in 1:ket_a.no
                for p in 1:ket_a.no
                    sign_a, La = ket_a_lookup[Ka][p+(r-1)*ket_a.no]
                    if La == 0
                        continue
                    end
                  
                    Lb = 1
                    sign_b = 1
                    L = 1 
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_lookup[Kb][q+(s-1)*ket_b.no]
                            
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * a_max
                           
                            Hmat[K,L] += H.h2[p,r,q,s] * sign_a * sign_b
                            continue
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end
    #sig = transpose(sig)
    return Hmat 
end
#=}}}=#


function compute_ab_terms(v, H::Hamiltonian, P::Problem)
    #={{{=#

    #print(" Compute opposite spin terms. Shape of v: ", size(v), "\n")
    @assert(size(v,1) == P.dim)
    
    #v = transpose(vin)
    
    sig = 0*v
  

    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    ket_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)
    bra_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    bra_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)

    ket_a_lookup = ConfigStrings.fill_ca_lookup(ket_a)
    ket_b_lookup = ConfigStrings.fill_ca_lookup(ket_b)
    
    a_max = bra_a.max
    ConfigStrings.reset!(ket_b)
                          
    n_roots = size(sig,2)
    scr = zeros(1,ket_a.max*ket_b.max)
    for Kb in 1:ket_b.max

        ConfigStrings.reset!(ket_a)
        for Ka in 1:ket_a.max
            K = Ka + (Kb-1) * ket_a.max

            #  <pq|rs> p'q'sr  --> (pr|qs) (a,b)
            for r in 1:ket_a.no
                for p in 1:ket_a.no
                    sign_a, La = ket_a_lookup[Ka][p+(r-1)*ket_a.no]
                    if La == 0
                        continue
                    end
                  
                    Lb = 1
                    sign_b = 1
                    L = 1 
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_lookup[Kb][q+(s-1)*ket_b.no]
                            
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * a_max
                           
                            #sig[K,:] += H.h2[p,r,q,s] * v[L,:]
                            #sig[K,:] += H.h2[p,r,q,s] * sign_a * sign_b * v[L,:]
                            for si in 1:n_roots
                                sig[K,si] += H.h2[p,r,q,s] * sign_a * sign_b * v[L,si]
                                #@views sig[K,si] .+= H.h2[p,r,q,s] * sign_a * sign_b * v[L,si]
                            end
                            continue
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end
    #sig = transpose(sig)
    return sig 
end
#=}}}=#


function compute_ab_terms(v, H::Hamiltonian, P::Problem, ket_a_lookup, ket_b_lookup)
    #={{{=#

    #print(" Compute opposite spin terms. Shape of v: ", size(v), "\n")
    @assert(size(v,1) == P.dim)
    
    #v = transpose(vin)
    
    sig = 0*v
  

    #   Create local references to ci_strings
    ket_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    ket_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)
    bra_a = ConfigStrings.ConfigString(no=P.no, ne=P.na)
    bra_b = ConfigStrings.ConfigString(no=P.no, ne=P.nb)

    a_max = bra_a.max
    ConfigStrings.reset!(ket_b)
                          
    n_roots = size(sig,2)
    scr = zeros(1,ket_a.max*ket_b.max)
    for Kb in 1:ket_b.max

        ConfigStrings.reset!(ket_a)
        for Ka in 1:ket_a.max
            K = Ka + (Kb-1) * ket_a.max

            #  <pq|rs> p'q'sr  --> (pr|qs) (a,b)
            for r in 1:ket_a.no
                for p in 1:ket_a.no
                    sign_a, La = ket_a_lookup[Ka][p+(r-1)*ket_a.no]
                    if La == 0
                        continue
                    end
                  
                    Lb = 1
                    sign_b = 1
                    L = 1 
                    for s in 1:ket_b.no
                        for q in 1:ket_b.no
                            sign_b, Lb = ket_b_lookup[Kb][q+(s-1)*ket_b.no]
                            
                            if Lb == 0
                                continue
                            end

                            L = La + (Lb-1) * a_max
                           
                            #sig[K,:] += H.h2[p,r,q,s] * v[L,:]
                            #sig[K,:] += H.h2[p,r,q,s] * sign_a * sign_b * v[L,:]
                            for si in 1:n_roots
                                sig[K,si] += H.h2[p,r,q,s] * sign_a * sign_b * v[L,si]
                                #@views sig[K,si] .+= H.h2[p,r,q,s] * sign_a * sign_b * v[L,si]
                            end
                            continue
                        end
                    end
                end
            end
            ConfigStrings.incr!(ket_a)

        end
        ConfigStrings.incr!(ket_b)
    end
    #sig = transpose(sig)
    return sig 
end
#=}}}=#


function precompute_spin_diag_terms(H::Hamiltonian, P::Problem, e)
    #={{{=#

    #   Create local references to ci_strings
    ket = ConfigStrings.ConfigString(no=P.no, ne=e)
    bra = ConfigStrings.ConfigString(no=P.no, ne=e)

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


function get_matvec_fn(ham::Hamiltonian, prb::Problem, HdiagA, HdiagB)
    #=
    Get function with takes a vector and returns action of H on that vector
    =#
#={{{=#
    ket_a = ConfigStrings.ConfigString(no=prb.no, ne=prb.na)
    ket_b = ConfigStrings.ConfigString(no=prb.no, ne=prb.nb)
    
    lookup_a = ConfigStrings.fill_ca_lookup(ket_a)
    lookup_b = ConfigStrings.fill_ca_lookup(ket_b)

    function (v)
        @time sig = compute_ab_terms(v, ham, prb, lookup_a, lookup_b)
        sig += HdiagA * v
        sig += HdiagB * v
        return sig 
    end
end
#=}}}=#


function get_map(ham::Hamiltonian, prb::Problem, HdiagA, HdiagB)
    #=
    Get LinearMap with takes a vector and returns action of H on that vector
    =#
    #={{{=#
    ket_a = ConfigStrings.ConfigString(no=prb.no, ne=prb.na)
    ket_b = ConfigStrings.ConfigString(no=prb.no, ne=prb.nb)
    
    lookup_a = ConfigStrings.fill_ca_lookup(ket_a)
    lookup_b = ConfigStrings.fill_ca_lookup(ket_b)

    function mymatvec(v)
        sig = compute_ab_terms(v, ham, prb)
        #sig = compute_ab_terms(v, ham, prb, lookup_a, lookup_b)
        sig += HdiagA * v
        sig += HdiagB * v
        return sig 
    end
    return LinearMap{Float64}(mymatvec, prb.dim; issymmetric=true, ismutating=false, ishermitian=true)
end
#=}}}=#

end
