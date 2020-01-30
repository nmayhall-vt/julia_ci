include("CIString.jl")
import .CIString


using TensorOperations
using NPZ
using Printf

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

ehf = CIString.slater_det_energy(ecore, ints_1b, ints_2b, n_elec_a, n_elec_b)
@printf(" HF Energy = %12.8f\n", ehf)

conf = Vector(1:n_elec_a)
ehf = CIString.slater_det_energy(ecore, ints_1b, ints_2b, conf, conf)
@printf(" HF Energy = %12.8f\n", ehf)

ket_a = CIString.ConfigString(no=n_orbs, ne=n_elec_a)
CIString.print(ket_a)
CIString.calc_max!(ket_a)
for i in 1:min(ket_a.max,10)
    CIString.incr!(ket_a)
    ket_b = deepcopy(ket_a)
    CIString.calc_linear_index!(ket_b)
end

CIString.apply_annihilation!(ket_a,6)
CIString.calc_max!(ket_a)
CIString.calc_linear_index!(ket_a)
CIString.print(ket_a)
@assert(ket_a.config == [1, 2, 3, 4, 7, 8])
@assert(ket_a.sign == -1)
@assert(ket_a.max == 210)
@assert(ket_a.lin_index == 10)

CIString.apply_creation!(ket_a,5)
CIString.calc_linear_index!(ket_a)
CIString.calc_max!(ket_a)
CIString.print(ket_a)
@assert(ket_a.config == [1, 2, 3, 4, 5, 7, 8])
@assert(ket_a.sign == -1)
@assert(ket_a.max == 120)
@assert(ket_a.lin_index == 5)

CIString.apply_creation!(ket_a,6)
CIString.calc_linear_index!(ket_a)
CIString.calc_max!(ket_a)
CIString.print(ket_a)
@assert(ket_a.config == [1, 2, 3, 4, 5, 6, 7, 8])
@assert(ket_a.sign == 1)
@assert(ket_a.max == 45)
@assert(ket_a.lin_index == 1)

lookup = CIString.fill_ca_lookup(ket_a)
#for i in lookup
#    print(i)
#    print("\n")
#end
