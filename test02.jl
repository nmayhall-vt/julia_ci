push!(LOAD_PATH, "./")
import ConfigStrings
import CISolvers
import Tools

using LinearAlgebra 
using Random 
using Printf
using NPZ

using Arpack
using IterativeSolvers, LinearMaps, Preconditioners

ecore = npzread("ints_0b.npy")
ints_1b = npzread("ints_1b.npy")
ints_2b = npzread("ints_2b.npy")


n_orbs = size(ints_1b,1)
n_elec_a = 2
n_elec_b = 2 
print(" Number of Orbitals = ", size(ints_1b))
print(size(ints_1b[1:n_elec_a,1:n_elec_a]))
@printf("\n")

H = CISolvers.Hamiltonian(h0=ecore, h1=ints_1b, h2=ints_2b,
                              no=n_orbs, na=n_elec_a, nb=n_elec_b)
#v = zeros(Real,10)
#v[1] = 1
#print(v)
n_roots = 1
v = zeros(H.dim,n_roots)
for i in 1:n_roots
    v[i,i] = 1
end
#normalize!(v)

sig = zeros(H.dim, n_roots)
v = CISolvers.compute_ab_terms!(v,sig,H)
exit()
ket_a = ConfigStrings.ConfigString(no=n_orbs, ne=n_elec_a)
ConfigStrings.print(ket_a)

state = CISolvers.State(no=n_orbs, na=n_elec_a, nb=n_elec_b;)
print(state)

#display(Tools.binom_coeff)

Random.seed!(0); randstring()

N = 8000
O = rand(N,N)-ones(N,N)*.5 - Diagonal(rand(N))*100
O = O + transpose(O)
function form_sigma(x; mat=O)
    return mat*x
end

A = LinearMap{Float64}(form_sigma, N; ismutating=false, ishermitian=true)


@time e,v = eigs(O, nev = 1, which=:SR)
@time e,v = eigs(O, nev = 1, which=:SR)
println(e)
@time e,v = eigs(A, nev = 1, which=:SR)
@time e,v = eigs(A, nev = 1, which=:SR)
println(e)

exit()

tol = 1e-5
@time e,v = powm(A; tol=tol, log=true, verbose=true) 
println(e)
@time x = lobpcg(A, false, 2; tol=tol, log=true) 
@time x = lobpcg(A, false, 2; tol=tol, log=true) 
println(x.trace)
e = x.λ
println(e)

@time ev = eigvals(O)[1:2]
println(ev)
