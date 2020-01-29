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
for i in 1:10
    CIString.incr!(ket_a)
    CIString.print(ket_a)
end
