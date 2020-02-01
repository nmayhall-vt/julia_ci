push!(LOAD_PATH, "./")
import ConfigStrings
import CISolvers
import Tools

using LinearAlgebra 
using Random 
using Printf
using NPZ

using TensorOperations
using Arpack
using IterativeSolvers, LinearMaps, Preconditioners

ecore = npzread("ints_0b.npy")
ints_1b = npzread("ints_1b.npy")
ints_2b = npzread("ints_2b.npy")


n_orbs = size(ints_1b,1)
n_elec_a = 1
n_elec_b = 1 
print(" Number of Orbitals = ", size(ints_1b))
print(size(ints_1b[1:n_elec_a,1:n_elec_a]))
@printf("\n")

H = CISolvers.Hamiltonian(h0=ecore, h1=ints_1b, h2=ints_2b,
                              no=n_orbs, na=n_elec_a, nb=n_elec_b)
#v = zeros(Real,10)
#v[1] = 1
#print(v)
n_roots = 2
v = zeros(H.dim,n_roots)
for i in 1:n_roots
    v[i,i] = 1
end
#normalize!(v)

state = CISolvers.State(no=n_orbs, na=n_elec_a, nb=n_elec_b;)
sig = zeros(H.dim, n_roots)
sig = CISolvers.compute_ab_terms(v,H)

@printf(" ecore: %12.8f\n",ecore)

Hdiag_a = CISolvers.precompute_spin_diag_terms(H,state.na)
Hdiag_b = CISolvers.precompute_spin_diag_terms(H,state.nb)

# v(I*J,s) -> v(I,J,s)
vin = Matrix(1.0I, state.dim, state.dim)

#vin = reshape(vin, state.dima, state.dimb, size(vin,2))
#println(size(vin))

#Hmat = zeros(state.dima, state.dimb, state.dim)
#@tensor begin
#    Hmat[a,c,d]  = Hdiag_a[a,b] * vin[b,c,d] 
#    Hmat[a,c,d] += Hdiag_b[c,a] * vin[a,a,d] 
#end
#vin  = reshape(vin, state.dima*state.dimb, size(vin,2))
#Hmat = reshape(Hmat, state.dima*state.dimb, size(vin,2))


Hmat  = kron(Hdiag_b, Matrix(1.0I, state.dima, state.dima))
Hmat += kron(Matrix(1.0I, state.dimb, state.dimb), Hdiag_a)

Hmat += CISolvers.matvec(vin , H)
@time e,v = eigs(Hmat, nev = 10, which=:SR)
e = e .+ ecore
for ei in e
    @printf(" %12.8f  %12.8fi\n", real(ei), imag(ei))
end


exit()

function matvec_wrapper(H)
    function (v)
        return CISolvers.matvec(v,H)
    end
end
wrapped = matvec_wrapper(H)
A = LinearMap{Float64}(wrapped, state.dim; ismutating=false, ishermitian=true)

@time e,v = eigs(A, nev = 10, which=:SR)
for ei in e
    @printf(" %12.8f  %12.8fi\n", real(ei), imag(ei))
end

Hmat = CISolvers.matvec(Matrix(1.0I, state.dim, state.dim) , H)
@time e,v = eigs(Hmat, nev = 10, which=:SR)
for ei in e
    @printf(" %12.8f  %12.8fi\n", real(ei), imag(ei))
end


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
