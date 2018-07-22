
#include("mem_variable.jl")
#include("mem_ftn_def.jl")




function p2ABC(p::Array{Float64},dim::Int64)

    A=zeros(Complex128,dim,dim)
    B=zeros(Complex128,dim,dim)
    C=zeros(Complex128,dim,dim)
    for i=1:dim
        A[i,i]=p[i]
        B[i,i]=p[dim+i]
    end
    for i=1:dim-1
        C[i,i]=p[2*dim+i]
    end
    C-=trace(C)/dim
    
     s=div(dim*(dim-1),2)
     offd=3*dim-1
     A += vec2utri( p[offd+1:(offd)+s],   p[offd+s+1:offd+2s]  ,dim)
     B += vec2utri( p[offd+2s+1:offd+3s], p[offd+3s+1:offd+4s] ,dim)
     C += vec2utri( p[offd+4s+1:offd+5s], p[offd+5s+1:offd+6s] ,dim)

    A= Hermitian(A)
    B= Hermitian(B)
    C= Hermitian(C)
    return (A,B,C)
end


function vec2utri(u::Array{Float64},v::Array{Float64},dim::Int64)
    A=zeros(Complex128,dim,dim)
    t=1
    for i=1:dim
        for j=i+1:dim
            A[i,j]=u[t]+im*v[t]
            t+=1
        end
    end
    return A
end


function utr2vec(A::Hermitian{Complex{Float64},Array{Complex{Float64},2}},dim::Int64)
    s=div((dim*(dim-1)),2)
    u=zeros(Float64,2s)
    t=1
    Ar= real(A)
    Ai= imag(A)
    for i=1:dim
        for j=i+1:dim
            u[t]=Ar[i,j]
            u[s+t]=Ai[i,j]
            t+=1
        end
    end
    return u
end



function gaussianPotential!(imagFreqFtn::strImagFreqFtn, realFreqFtn::strRealFreqFtn, kernel::strKernel, numeric::strNumeric, xinit, Aw::Array{Array{Complex128,2}})
Egrid = numeric.Egrid
dim=size(imagFreqFtn.moments1)[1]
##=
function getABC2(F::Array{Float64}, p::Array{Float64})
    (A,B,C)=p2ABC(p,dim)

    E0=Array{Float64}(numeric.Egrid)

    local_trace_value =  Array{Float64}(Egrid)
    for w=1:Egrid
#        realFreqFtn.H_extern[w] = Hermitian( A * kernel.moment[3,w] + B* kernel.moment[2,w] + C * kernel.moment[1,w]  ) 
        realFreqFtn.H_extern[w] =  - Hermitian( A * numeric.ERealAxis[w]^2 + B* numeric.ERealAxis[w]+ C  ) 
        E0[w]=eigmin(Hermitian(realFreqFtn.H_extern[w]))
    end
    chempot = minimum(E0)
    for w=1:Egrid
        Aw[w], local_trace_value[w]       = density_matrix_constructor( realFreqFtn.H_extern[w] -chempot *eye(realFreqFtn.H_extern[w]))
#        Aw[w] *= numeric.dERealAxis[w]
#        local_trace_value[w] *= numeric.dERealAxis[w]
        Aw[w] = Hermitian(Aw[w])
    end
    #sum_rule
    TotalPartitionFtn =  kernel.moment[1,:]'local_trace_value
    for w=1:numeric.Egrid
      Aw[w] *= 1./TotalPartitionFtn
    end

    mw = kernel.moment * Aw
    x=Hermitian(mw[1])
    y=Hermitian(mw[2])
    z=Hermitian(mw[3])

    x=Hermitian(imagFreqFtn.moments1-  mw[1])
    y=Hermitian(imagFreqFtn.moments2-  mw[2])
    z=Hermitian(imagFreqFtn.moments3-  mw[3])
    x[dim,dim]=0

    for i=1:dim
      F[i]     = y[i,i]
      F[dim+i] = z[i,i]
    end
    for i=1:dim-1
      F[2dim+i] = x[i,i]
    end

    s=div(dim*(dim-1),2)
    offd = 3dim-1
    F[offd+1:offd+2s] = utr2vec(y,dim)
    F[offd+2s+1:offd+4s] = utr2vec(z,dim)
    F[offd+4s+1:offd+6s] = utr2vec(x,dim)

end
# =#



result=nlsolve(getABC2, xinit, xtol = 1e-5, ftol=1e-10 )
xfinal = result.zero


(A, B, C) = p2ABC(xfinal,dim)

println( result  )
println("Model gaussian A: ", A)
println("Model gaussian B: ", B)
println("Model gaussian C: ", C)

return xfinal

end
