push!(LOAD_PATH, "./")
import ConfigStrings
import CISolvers
import Tools

using FCI 

using LinearAlgebra 
using Random 
using Printf
using NPZ

using TensorOperations
using Arpack
using LinearMaps 
#using KrylovKit 
#using IterativeSolvers, LinearMaps, Preconditioners
using IterativeSolvers

ecore = npzread("ints_0b.npy")
ints_1b = npzread("ints_1b.npy")
ints_2b = npzread("ints_2b.npy")

n_elec_a = 1
n_elec_b = 1

n_orbs = size(ints_1b,1)
print(" Number of Orbitals = ", size(ints_1b))
print(size(ints_1b[1:n_elec_a,1:n_elec_a]))
@printf("\n")


Random.seed!(0); randstring()
function get_rand_mat(N)
    O = rand(N,N)-ones(N,N)*.5 - Diagonal(rand(N))*100
    O = O + transpose(O)
    return O
end

ham = FCI.Hamiltonian(ecore, ints_1b, ints_2b)
problem = FCI.Problem(no=size(ints_1b,1), na=n_elec_a, nb=n_elec_b)

A = get_rand_mat(problem.dim)

#H = FCIMap(A)

v = Matrix(1.0I, problem.dim, problem.dim)[:,1]


#mymatvec = FCI.matvec(A,v)
#
#A = LinearMap{Float64}(mymatvec, size(A,1); ismutating=false, ishermitian=true)
#
#@time e,v = eigs(A, nev = 1, which=:SR)
#println(e)

print(" Compute spin_diagonal terms\n")
@time Hdiag_a = FCI.precompute_spin_diag_terms(ham,problem,n_elec_a)
@time Hdiag_b = FCI.precompute_spin_diag_terms(ham,problem,n_elec_b)
print(" Kron them")
@time Hdiag_a = kron(Matrix(1.0I, problem.dimb, problem.dimb), Hdiag_a)
@time Hdiag_b = kron(Hdiag_b, Matrix(1.0I, problem.dima, problem.dima))
print(" done\n")
print(" done\n")

H = FCI.get_Map(ham, problem, Hdiag_a, Hdiag_b)

#@time e,v = powm(H; tol=1e-6, log=true, verbose=true) 
@time e,v = eigs(H, v0=v, nev = 1, which=:SR)

for ei in e
    @printf(" %12.8f\n",ei+ham.h0)
end
