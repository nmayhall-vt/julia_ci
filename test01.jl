include("ConfigStrings.jl")
import .ConfigStrings


using TensorOperations
using NPZ
using Printf

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

ecore = npzread("ints_0b.npy")
ints_1b = npzread("ints_1b.npy")
ints_2b = npzread("ints_2b.npy")

n_orbs = size(ints_1b,1)
n_elec_a = 7
n_elec_b = 7 
print(" Number of Orbitals = ", size(ints_1b))
print(size(ints_1b[1:n_elec_a,1:n_elec_a]))
@printf("\n")

#@tensor begin
#	E1 = 2*ints_1b[1:n_elec_a,1:n_elec_a][a,a]
#	E2 = 2*ints_2b[1:n_elec_a,1:n_elec_a,1:n_elec_a,1:n_elec_a][a,a,b,b]
#	E2 -= 1*ints_2b[1:n_elec_a,1:n_elec_a,1:n_elec_a,1:n_elec_a][a,b,b,a]
#end
#@printf(" HF Energy = %12.8f\n", E1+E2+ecore)

ehf = slater_det_energy(ecore, ints_1b, ints_2b, n_elec_a, n_elec_b)
@printf(" HF Energy = %12.8f\n", ehf)

conf = Vector(1:n_elec_a)
ehf = slater_det_energy(ecore, ints_1b, ints_2b, conf, conf)
@printf(" HF Energy = %12.8f\n", ehf)

ket_a = ConfigStrings.ConfigString(no=n_orbs, ne=n_elec_a)
ehf = ConfigStrings.slater_det_energy(ecore, ints_1b, ints_2b, ket_a, ket_a)
@printf(" HF Energy = %12.8f\n", ehf)


ConfigStrings.print(ket_a)
ConfigStrings.calc_max!(ket_a)
for i in 1:min(ket_a.max,10)
    ConfigStrings.incr!(ket_a)
    ket_b = deepcopy(ket_a)
    ConfigStrings.calc_linear_index!(ket_b)
end

ConfigStrings.apply_annihilation!(ket_a,6)
ConfigStrings.calc_max!(ket_a)
ConfigStrings.calc_linear_index!(ket_a)
ConfigStrings.print(ket_a)
@assert(ket_a.config == [1, 2, 3, 4, 7, 8])
@assert(ket_a.sign == -1)
@assert(ket_a.max == 210)
@assert(ket_a.lin_index == 10)

ConfigStrings.apply_creation!(ket_a,5)
ConfigStrings.calc_linear_index!(ket_a)
ConfigStrings.calc_max!(ket_a)
ConfigStrings.print(ket_a)
@assert(ket_a.config == [1, 2, 3, 4, 5, 7, 8])
@assert(ket_a.sign == -1)
@assert(ket_a.max == 120)
@assert(ket_a.lin_index == 5)

ConfigStrings.apply_creation!(ket_a,6)
ConfigStrings.calc_linear_index!(ket_a)
ConfigStrings.calc_max!(ket_a)
ConfigStrings.print(ket_a)
@assert(ket_a.config == [1, 2, 3, 4, 5, 6, 7, 8])
@assert(ket_a.sign == 1)
@assert(ket_a.max == 45)
@assert(ket_a.lin_index == 1)

lookup = ConfigStrings.fill_ca_lookup(ket_a)
#for i in lookup
#    print(i)
#    print("\n")
#end
