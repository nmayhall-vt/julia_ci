push!(LOAD_PATH, "./")
import ConfigStrings
import CISolvers
import Tools

using FCI 

using LinearAlgebra 
using Random 
using Printf
using NPZ

using Profile 

using TensorOperations
using Arpack
using LinearMaps 
#using KrylovKit 
#using IterativeSolvers, LinearMaps, Preconditioners
#using IterativeSolvers

ecore = npzread("ints_0b.npy")
ints_1b = npzread("ints_1b.npy")
ints_2b = npzread("ints_2b.npy")

n_elec_a = 4
n_elec_b = 4

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

#A = get_rand_mat(problem.dim)

#H = FCIMap(A)

v = Matrix(1.0I, problem.dim, problem.dim)


#mymatvec = FCI.matvec(A,v)
#
#A = LinearMap{Float64}(mymatvec, size(A,1); ismutating=false, ishermitian=true)
#
#@time e,v = eigs(A, nev = 1, which=:SR)
#println(e)

Hmat = zeros(Float64, problem.dim, problem.dim)

print(" Compute spin_diagonal terms\n")
Hdiag_a = FCI.precompute_spin_diag_terms(ham,problem,problem.na)
Hdiag_b = FCI.precompute_spin_diag_terms(ham,problem,problem.nb)
print(" done\n")
print(" Kron them")
Hdiag_a = kron(Matrix(1.0I, problem.dimb, problem.dimb), Hdiag_a)
Hdiag_b = kron(Hdiag_b, Matrix(1.0I, problem.dima, problem.dima))
print(" done\n")

Hmap = FCI.get_map(ham, problem, Hdiag_a, Hdiag_b)


@time FCI.compute_spin_diag_terms_full!(ham, problem, Hmat)
#println(Hmat)
@time FCI.compute_ab_terms_full!(ham, problem, Hmat)
#@time Hmat += FCI.compute_ab_terms_full(ham, problem)
#Hmat += Hmap*v

#Hmat = .5*(Hmat + transpose(Hmat))
@time e,v = eigs(Hmat, v0=v[:,1], nev = 1, which=:SR)
e = real(e)
for ei in e
    @printf(" %12.8f  %12.8f\n",ei,ei+ham.h0)
end

exit()





#profile
function test(H,v)
    for i in 1:2
        a = H*v
    end
end
#@time test(H,v)
@time test(H,v)
@time test(H,v)
@time test(H,v)
#@profile test(H,v)
#@profview test(H,v)
#Profile.print(sortedby=:count)
#Profile.print()                      # Prints results from Profile.print()
#Profile.print()                      # Prints results from Profile.print()
exit()


#@time e,v = powm(H; tol=1e-6, log=true, verbose=true) 
@time e,v = eigs(H, v0=v, nev = 1, which=:SR)

for ei in e
    @printf(" %12.8f\n",ei+ham.h0)
end
