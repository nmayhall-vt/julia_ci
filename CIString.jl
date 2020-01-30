module CIString
using Printf
using Parameters

function slater_det_energy(h0,h1,h2,na::Int,nb::Int)
	#Compute the energy of a Slater Det with na (nb) alpha (beta) electrons
#={{{=#
	E0 = h0 
	E1 = 0
	E2 = 0
	for i=1:na
		E1 += h1[i,i]
		for j=i:na
			E2 += h2[i,i,j,j]
			E2 -= h2[i,j,j,i]
		end
		for j=1:nb
			E2 += h2[i,i,j,j]
		end
	end
	for i=1:nb
		E1 += h1[i,i]
		for j=i:nb
			E2 += h2[i,i,j,j]
			E2 -= h2[i,j,j,i]
		end
	end
	return E0+E1+E2
end
#=}}}=#

function slater_det_energy(h0,h1,h2,na::Array{Int,1},nb::Array{Int,1})
	#Compute the energy of a Slater Det specified by 
	# alpha (beta) string, na (nb)
	#={{{=#
	E0 = h0 
	E1 = 0
	E2 = 0
	
	@assert(length(Set(na))==length(na))
	@assert(length(Set(nb))==length(nb))
	
	for ii=1:length(na)
		i = na[ii]
		E1 += h1[i,i]
		for jj=ii:length(na)
			j = na[jj]
			E2 += h2[i,i,j,j]
			E2 -= h2[i,j,j,i]
		end
		for jj=1:length(nb)
			j = nb[jj]
			E2 += h2[i,i,j,j]
		end
	end
	for ii=1:length(nb)
		i = nb[ii]
		E1 += h1[i,i]
		for jj=i:length(nb)
			j = nb[jj]
			E2 += h2[i,i,j,j]
			E2 -= h2[i,j,j,i]
		end
	end
	return E0+E1+E2
end
#=}}}=#

@with_kw mutable struct ConfigString
    #=
    Type to organize all the configuration string 
    =#
# {{{
	no::Int = 0
	ne::Int = 0
	sign::Int = 1
	lin_index::Int = 1
	config::Array{Int,1} = Vector(1:ne)
	#ca_lookup::Array{Array{Int,1},1}
        max::Int = calc_nchk(no,ne) 
end
## }}}

import Base: length 
function length(c::ConfigString)
        #=
	return number of strings
	=#
	return c.max
end


import Base: print 
function print(c::ConfigString)
	#= 
	Pretty print of an determinant string
	=#
	#={{{=#
        @printf("Index: %-10i NOrb: %-4i Dim: %-10i Sign: %2i ",c.lin_index, c.no, calc_nchk(c.no,c.ne), c.sign)
	print(c.config)
	print('\n')
end
#=}}}=#

function incr!(c::ConfigString)
	#= 
	Increment determinant string
	=#
#={{{=#
    if c.max == nothing
        calc_max!(c)
    end
    if c.lin_index == c.max
        return
    end
    c.lin_index += 1
    incr_comb!(c.config, c.no)
end
#=}}}=#

function calc_max!(c::ConfigString)
    #=
    Calculate dimension of space accessible to a ConfigString
    =#
    c.max = calc_nchk(c.no,c.ne)
end


function calc_nchk(n::Int,k::Int)
    #= 
    Calculate n choose k
    =#
    #={{{=#
    @assert(n>k)
    accum::BigInt = 1
    for i in 1:k
        accum = accum * (n-k+i) ÷ i
    end
    return accum
end
#=}}}=#

function incr_comb!(comb::Array{Int,1}, Mend::Int)
    #=
    For a given combination, form the next combination
    =#
#={{{=#
    N = length(comb)
    for i in N:-1:1
        if comb[i] < Mend - N + i 
            comb[i] += 1
            for j in i+1:N
                comb[j]=comb[j-1]+1
            end
            return
        end
    end
    return
end
#=}}}=#


function calc_linear_index!(c::ConfigString)
    #=
    Return linear index for lexically ordered __config string
    =#
#={{{=#
    c.lin_index = 1
    v_prev::Int = 0 

    for i in 1:c.ne 
        v = c.config[i]
        #todo: change mchn from function call to data lookup!
        for j in v_prev+1:v-1
            c.lin_index += calc_nchk(c.no-j,c.ne-i)
        end
        v_prev = v
    end
    return
end
#=}}}=#


function fill_ca_lookup!(c::ConfigString)
    #=
    Create an index table relating each string with all ia substitutions
        i.e., ca_lookup[Ka][c(p) + a(p)*n_p] = La
    =#
#={{{=#

    ket = CIString.ConfigString(no=c.no, ne=c.ne)
    bra = CIString.ConfigString(no=c.no, ne=c.ne)

    calc_max!(ket)
    
    for K in 1:ket.max
        Kv = Array{Int,1}
        for p in 1:ket.no
            for q in 1:ket.no
                bra = deepcopy(ket)

            end
        end
        incr!(ket)
    end
end
#=}}}=#


function reset!(c::ConfigString)
    c.config = Vector(1:ne)
    c.sign = 1
    c.lin_index = 1
end
function destroy_config!(c::ConfigString)
    c.config = Missing 
    c.sign = 0
    c.lin_index = 0
    c.max = 0
    c.ne = 0
    c.no = 0
end
function destroy_config!(c::ConfigString)
    c.config = Missing 
    c.sign = 0
    c.lin_index = 0
    c.max = 0
    c.ne = 0
    c.no = 0
end


function a!(c::ConfigString, orb_index::Int)
    #=
    apply annihilation operator a_i to current string
    where orb_index is i
    =#
    #={{{=#
    @assert(orb_index < c.no)
    if c.sign == 0
        return
    end
    found = -1
    for i in 1:c.ne
        if c.config[i] == orb_index
            found = i
            break
        end
    end
   
    
    if found == -1
        destroy_config!(c)
        return
    end

    if found % 2 != 0
        c.sign *= -1
    end

    deleteat!(c.config,found)
    
    c.ne -= 1
    
    #unset data that need to be recomputed
    c.max = 0
    c.lin_index = 0
end
#=}}}=#


function c!(c::ConfigString, orb_index::Int)
    #=
    apply creation operator a_i to current string
    where orb_index is i
    =#
    #={{{=#
    @assert(orb_index < c.no)
    if c.sign == 0
        return
    end
    
    insert_here = -1
    for i in 1:c.ne
        if c.config[i] > orb_index
            insert_here = i
            break
        elseif c.config[i] == orb_index
            destroy_conf!(c)
            return
        else
            insert_here += 1
        end
    end
   
    if insert_here % 2 != 0
        c.sign *= -1
    end

    insert!(c.config,insert_here,orb_index)
    
    c.ne += 1
    #unset data that need to be recomputed
    c.max = 0
    c.lin_index = 0
end
#=}}}=#


end

