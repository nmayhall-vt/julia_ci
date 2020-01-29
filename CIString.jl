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
	max::Int = 1
	lin_index::Int = 1
	config::Array{Int,1} = Vector(1:ne)
	#ca_lookup::Array{Array{Int,1},1}
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
	@printf("Index:%-12i NO:%4i Sign:%2i ",c.lin_index, c.no, c.sign)
	print(c.config)
	print('\n')
end
#=}}}=#

function incr!(c::ConfigString)
	#= 
	Increment determinant string
	=#
#={{{=#
	if c.lin_index == c.max-1
		return
        end
        c.lin_index += 1
        incr_comb!(c.config, c.no)
end
#=}}}=#


function incr_comb!(comb::Array{Int,1}, Mend::Int)
    #=
    For a given combination, form the next combination
    =#
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

end

