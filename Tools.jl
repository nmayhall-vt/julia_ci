module Tools

function calc_nchk(n::Int,k::Int)
    #= 
    Calculate n choose k
    =#
    #={{{=#
    @assert(n>=k)
    accum::BigInt = 1
    for i in 1:k
        accum = accum * (n-k+i) ÷ i
    end
    return accum
end
#=}}}=#


N = 20
binom_coeff = Array{BigInt,2}(undef,N,N)
for i in 1:N
    for j in i:N
        binom_coeff[j,i] = calc_nchk(j,i)
    end
end


end
