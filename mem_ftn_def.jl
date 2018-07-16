###########################################
# Jae-Hoon Sim, KAIST 2018.01.
###########################################
#############################################################################################
# (start) define Ftns
#############################################################################################

include("mem_variable.jl")
include("mem_gaussianPot.jl")







global cubic_coeff_transf 
function get_total_energy(imagFreqFtn ,Aw::Array{Array{Complex128,2}}, kernel, numeric)
        NumSubOrbit        = size(Aw[1])[1]
        Energy=0
        EnergyMatrix=zeros(Complex128, NumSubOrbit, NumSubOrbit)
        for n=1:numeric.N_Matsubara
           GreenReproduced=zeros(Complex128, NumSubOrbit, NumSubOrbit)
           for w=1:numeric.Egrid
              GreenReproduced += (kernel.Kernel[n,w]*Aw[w])
           end
           EnergyMatrix = (imagFreqFtn.GreenFtn[n] - GreenReproduced)
           Energy    +=    (vecnorm(EnergyMatrix)^2)
        end
        moments =  kernel.moment * Aw
        temp = [ imagFreqFtn.moments1, imagFreqFtn.moments2, imagFreqFtn.moments3, imagFreqFtn.moments4, imagFreqFtn.moments5 ] - moments
        Energy_tail =  (real(trace(copy(temp') * kernel.gamma * temp)))
        Energy  += Energy_tail

        return Energy
end





function matrixFtn_norm( Aw::Array{Array{Complex128,2}} , weight::Array{Float64})
            # 1-norm, \int dw abs(A) = sum(abs(A))
            lengthOfVector = length(Aw)
            normVector_temp = Array{Float64}(lengthOfVector)
            for w=1:lengthOfVector
                  normVector_temp[w] =  (norm(Aw[w])) *weight[w]
            end
            return ( sum(normVector_temp))
end




function construct_model_spectrum(Green_Inf::Hermitian{Complex{Float64},Array{Complex{Float64},2}}, numeric, NumSubOrbit::Int64,  center::Float64,  eta::Float64 , shape::String)
    Spectral_default_model =  Array{Array{Complex128,2}}(numeric.Egrid)
   for w = 1:numeric.Egrid
       E = numeric.ERealAxis[w]
       if(shape == "G")
         Spectral_default_model[w] = Hermitian( (Green_Inf *( exp(-0.5 *((E-center)/eta)^2) + (1e-30)    ) ) )
       elseif(shape == "L")
         Spectral_default_model[w] = (Green_Inf * 1/( (E-center)^2 + eta^2) +1e-30*(eye(Green_Inf)) )
       end
   end
   model_norm = trace(sum(Spectral_default_model))

#no
#   Spectral_default_model  /= model_norm
#yes
  Spectral_default_model  ./= model_norm
  return Spectral_default_model
end







function Hamiltonian_zeroSet( imagFreqFtn, kernel,
                              Egrid::Int64,Spectral_default_model::Array{Array{Complex128,2}}, auxiliary_inverse_temp::Float64)
                                   
    Hamiltonian_out = -auxiliary_inverse_temp * ( (kernel.Kernel_dagger)*(imagFreqFtn.GreenFtn) ) ;    
    # d(G-KA)^2 =  -K'(G - KA) , second term,i.e. K'K,  = interaction term, 

#    temp = kernel.gamma * [ imagFreqFtn.moments1, imagFreqFtn.moments2, imagFreqFtn.moments3, imagFreqFtn.moments4, imagFreqFtn.moments5 ]
    for w=1:Egrid
#        Hamiltonian_out[w] -=  auxiliary_inverse_temp *   (kernel.moment[:,w])'  * temp
        Hamiltonian_out[w] -=  logm( (Spectral_default_model[w]) + 1e-15*eye(Spectral_default_model[w]) )
        Hamiltonian_out[w] +=  Hamiltonian_out[w]'    ;
    end   
    Hamiltonian_out *=      (0.5);
    for w=1:Egrid
      Hamiltonian_out[w] = Hermitian(Hamiltonian_out[w])
    end

    return Hamiltonian_out
end

function Hamiltonian_update(Aw::Array{Array{Complex128,2}},  Hamiltonian_zero::Array{Array{Complex128,2}}, numeric::strNumeric, auxiliary_inverse_temp::Float64, kernel::strKernel)
    Hamiltonian_out = deepcopy(Hamiltonian_zero)
    Egrid = numeric.Egrid
    NumSubOrbit = size(Aw[1])[1]

    ATran = zeros(Complex128, Egrid);
    for i=1:NumSubOrbit*NumSubOrbit
        for w=1:Egrid
            ATran[w] = Aw[w][i]
        end
        h_updateTran  = kernel.Kernel*ATran
        h_updateTran  = auxiliary_inverse_temp * (kernel.Kernel_dagger*(h_updateTran))::Array{Complex128}
        for w=1:Egrid
            Hamiltonian_out[w][i] += h_updateTran[w] 
        end
    end

    for w=1:Egrid
        Hamiltonian_out[w] +=  Hamiltonian_out[w]'
        Hamiltonian_out[w] = Hermitian(Hamiltonian_out[w]/2.0)
    end

    return Hamiltonian_out #rw
end

@everywhere function density_matrix_constructor(Hamiltonian_w::Array{Complex128,2} )
    Aw_pmap = Hermitian(expm(Hermitian( -(Hamiltonian_w)))) 
    tracevalue =   real(trace(Aw_pmap))
    return  Aw_pmap, tracevalue 
end

function spectralFunction_constructor(Hamiltonian::Array{Array{Complex128,2}}, NumSubOrbit::Int64, numeric::strNumeric, kernel::strKernel)

    local_trace_value = Array{Float64}(numeric.Egrid)
    Aw_out =           Array{Array{Complex128,2}}(numeric.Egrid)
        E0=Array{Float64}(numeric.Egrid)
        if(NumSubOrbit ==1)
          for w=1:numeric.Egrid
              E0[w] = real(Hermitian(Hamiltonian[w])[1])
          end
        else 
          for w=1:numeric.Egrid
              E0[w]=eigmin(Hermitian(Hamiltonian[w]))
          end
        end
        chempot = minimum(E0)
        for w=1:numeric.Egrid
            Hamiltonian[w] -= chempot*eye(Hamiltonian[w])
        end

        #Construct A(w) from input H
        for  w=1:numeric.Egrid
            Aw_out[w], local_trace_value[w] = density_matrix_constructor(Hamiltonian[w])
            Aw_out[w] = Hermitian(Aw_out[w])
        end
        #sum_rule
        TotalPartitionFtn =  kernel.moment[1,:]'local_trace_value
        for w=1:numeric.Egrid
          Aw_out[w] *= 1./TotalPartitionFtn
        end
        return Aw_out
end




function find_optimal_mixing_weight!(mixing_local::Float64,    criteria_Aw::Float64, criteria_Aw_prev100::Float64)
                                                                    
   if(criteria_Aw < criteria_Aw_prev100)
          mixing_local  *= 1.02
   elseif(criteria_Aw > criteria_Aw_prev100 ) 
          mixing_local  /= 1.05
   end
   mixing_local = min(mixing_local, mixing_max)
   mixing_local = max(mixing_local, mixing_min)
   return mixing_local 
end



function Aw_Iteration(realFreqFtn::strRealFreqFtn, imagFreqFtn::strImagFreqFtn, kernel::strKernel, auxiliary_inverse_temp::Float64, NumSubOrbit::Int64, numeric::strNumeric, mixing::Float64)
    #Start iteration loop
    MixingInformation = pulayInfo(numeric.pulay_mixing_history ,numeric.pulay_mixing_start,  numeric.pulay_mixing_step,
                        Array{Array{Array{Complex128,2}}}(0), Array{Array{Array{Complex128,2}}}(0), Array{Array{Array{Complex128,2}}}(0),Array{Float64}(0), "simple" )
    pulayInitializer!(MixingInformation)
    converg = false; 
    iter=0;


    criteria_Aw = 0;
    totEin = 0;

    criteria_Aw_Initial =0.0
    criteria_Aw_prev100 = 0.0
    totE_prev100 = 0.0

    mixing_local = mixing
    Aw_out =           Array{Array{Complex128,2}}(numeric.Egrid)
    Aw_in =            Array{Array{Complex128,2}}(numeric.Egrid)
    Hamiltonian_in =            Array{Array{Complex128,2}}(numeric.Egrid)

    #Initial set
    Aw_in = deepcopy(realFreqFtn.Aw)
    Hamiltonian_zero =  Hamiltonian_zeroSet(imagFreqFtn, kernel, numeric.Egrid   ,realFreqFtn.Spectral_default_model, auxiliary_inverse_temp)
    #Start iteration
    while true
        #new Hamiltonian, H[A(w)]
        Hamiltonian_in =  Hamiltonian_update(   Aw_in, Hamiltonian_zero, numeric, auxiliary_inverse_temp, kernel) 
        Aw_out = spectralFunction_constructor(Hamiltonian_in,  NumSubOrbit, numeric,kernel)

        if iter==0
            #estimate criteria
            resdAw = Aw_in - Aw_out
            criteria_Aw_prev100 =  matrixFtn_norm(resdAw, numeric.dERealAxis) ;
        end

        if iter%30 == 0  || iter == 1
            #estimate criteria
            resdAw = Aw_in - Aw_out
            criteria_Aw =  matrixFtn_norm(resdAw, numeric.dERealAxis) ;

              if(  criteria_Aw < 1e-10*NumSubOrbit     )
                converg = true
                realFreqFtn.Hw[:]  = Hamiltonian_in
                break
              elseif(   ( (iter > NumIter &&  criteria_Aw > criteria_Aw_prev100) || criteria_Aw > criteria_Aw_prev100 *10)    )
                converg = false
                break
              end
            mixing_local = find_optimal_mixing_weight!(mixing_local, criteria_Aw, criteria_Aw_prev100)
            criteria_Aw_prev100 = criteria_Aw 
        end 
        #Prepare next iteration
        Aw_in = pulayMixing!(iter, mixing_local, Aw_in,   Aw_out-Aw_in, MixingInformation)
        iter += 1
    end #end of iteration
    return ( Aw_in  ,   converg,   criteria_Aw, iter)
end









function  mem_annealing( kernel,  realFreqFtn, imagFreqFtn,   auxTempRenomalFactorInitial, write_information, fname_out)

    start_temperature = kernel.auxiliary_inverse_temp_range[1]
    end_temperature =   kernel.auxiliary_inverse_temp_range[2]
    logEnergy =  Array{Float64}(0)
    logAlpha  =  Array{Float64}(0)
    alpha_list = Array{Float64}(0)
    mixing =  0.0::Float64
    Aw = deepcopy(realFreqFtn.Aw)
    auxTempRenomalFactor = deepcopy(auxTempRenomalFactorInitial)
    NumSubOrbit = size(imagFreqFtn.GreenFtn[1])[1]
    
    
    
    
    
    
    mixing = deepcopy(numeric.mixingInitial)
    auxiliary_inverse_temp = deepcopy(start_temperature);
    
    f=open("$(workDirect)/information.out","a")
    write(f,"\n" )
    close(f)   
    
    auxiliary_inverse_temp_prev = deepcopy(auxiliary_inverse_temp)
    
    
    
    
    mixingPrev = deepcopy(mixing)
    mixingNext =  deepcopy(mixing)
    iter=0;
    while  auxiliary_inverse_temp  <=  end_temperature
        mixing = mixingNext
        auxTempRenomalFactor = min(1+ ((auxTempRenomalFactor-1)*1.01)   , auxTempRenomalFactorInitial)
        auxiliary_inverse_temp *= auxTempRenomalFactor;  
    
    
    
        (Aw,  converg,   normAwRD, iter) = Aw_Iteration(realFreqFtn, imagFreqFtn, kernel, auxiliary_inverse_temp, NumSubOrbit, numeric, mixing)
        mixingNext =  min( mixing*1.001, numeric.mixing_max)
        if !(converg)
          mixingNext =   max( mixing / 1.1, numeric.mixing_min)
          auxTempRenomalFactor = max(1+ ((auxTempRenomalFactor-1)/1.1), 1.001)
        end
        trial = 0
        while !(converg)
          trial+=1
          println("not converged... $(auxiliary_inverse_temp_prev)\t$(auxiliary_inverse_temp)\tmixing:$(mixing); Iter:$(iter)  normAw:$(normAwRD) trial:$(trial)")
          auxiliary_inverse_temp = 0.05 * (auxiliary_inverse_temp) + 0.95* (auxiliary_inverse_temp_prev);
        (Aw, converg,   normAwRD,iter) = Aw_Iteration(realFreqFtn, imagFreqFtn,kernel, auxiliary_inverse_temp, NumSubOrbit, numeric,mixing)
          if trial==10 break end
        end
        if trial==10 break end
    
       realFreqFtn.Aw[:] =  Aw
    
    
        auxiliary_inverse_temp_prev = deepcopy(auxiliary_inverse_temp)
    
        Energy = get_total_energy( imagFreqFtn, realFreqFtn.Aw, kernel,numeric)
        ##=  comment:write
        if write_information
            push!(logEnergy, log(Energy))
            push!(logAlpha, log(auxiliary_inverse_temp))
            push!(alpha_list, auxiliary_inverse_temp)
        
            f=open("$(workDirect)/information.out","a")
            write(f,"$auxiliary_inverse_temp\t\t$(Energy) ; $(iter)\t$(mixing)\t$(auxTempRenomalFactor)   \n" )
            close(f)   
#           for w=1:numeric.Egrid
#                     Aw[w] =  copy(bestUnitary) * realFreqFtn.Aw[w]  * (bestUnitary')
#           end
           write_spectral_ftn(NumSubOrbit, imagFreqFtn.Normalization, numeric, Aw, fname_out,  "")
           # comment:Reproduce
            GreenFtn_rep = kernel.Kernel *  Aw ;
            moments_rep  = kernel.moment *  Aw ;
#            println(" $(-real(trace(moments_rep[1])))/x +$(real(trace(moments_rep[3])))/x**3"  )
    
            
            if inversion_method
              for iw=1:numeric.N_Matsubara
                 z = ((2.*iw-1)*pi/inverse_temp)im;
                 GreenFtn_rep[iw] *= trace(eye(GreenFtn_rep[iw]))
                 GreenFtn_rep[iw] = (    z*eye(GreenFtn_rep[iw])-inv(GreenFtn_rep[iw])    )
              end
            end
    
            GreenFtn_rep *= imagFreqFtn.Normalization ;
            for i = 1:NumSubOrbit
                for j=1:NumSubOrbit
                f=open( "$(workDirect)/reproduce_$(i+startOrbit)_$(j+startOrbit).out","w")
                for iw=1:numeric.N_Matsubara
                    z = ((2.*iw-1)*pi/inverse_temp)im;
                    ftn = GreenFtn_rep[iw][i,j]/sqrt(Energy_weight[iw]) + imagFreqFtn.GreenConst[i,j] 
                    write(f,"$(imag(z))  $(real(ftn))    $(imag(ftn))  \n" )
                end
                end
                close(f)
            end
    
    
            f=open("$(workDirect)/cubic_coeff.out","w")
             temp = kernel.cubic_coeff_transf * Aw
             for w=1:numeric.Egrid-1
               for subinex = 0:19
                E = (20-subinex)/20 * numeric.ERealAxis[w] + subinex/20 * numeric.ERealAxis[w+1]
                interpol_val =  temp[4*(w-1)+1][1,1]*(E-numeric.ERealAxis[w])^3 + 
                                temp[4*(w-1)+2][1,1]*(E-numeric.ERealAxis[w])^2 + 
                                temp[4*(w-1)+3][1,1]*(E-numeric.ERealAxis[w])^1 + 
                                temp[4*(w-1)+4][1,1]*(E-numeric.ERealAxis[w])^0
                
                interpol_val = real(interpol_val)*imagFreqFtn.Normalization
                write(f,"$(E) $(interpol_val)\n")  
                end
             end
             close(f)
        end #write_information
        # =#
    end #while_ (auxiliary_inverse_temp  <=  end_temperature)

#    for w=1:numeric.Egrid
#              realFreqFtn.Aw[w] =  bestUnitary * realFreqFtn.Aw[w]  * bestUnitary'
#    end
    return  logEnergy,  logAlpha
end













function construct_EnergyGrid( EwinOuterLeft::Float64, EwinOuterRight::Float64, EwinInnerLeft::Float64, EwinInnerRight::Float64, EgridInner::Int64)
#Real freq. grid [PRE 94, 023303,  wol =0]
ERealAxis =  Array{Float64}(0);
dERealAxis =  Array{Float64}(0);

DeltaW    = (EwinInnerRight - EwinInnerLeft)/EgridInner #>0
coreDelta = DeltaW  / (coreDense)

##left

EL = EwinInnerLeft
EL1 = EwinOuterLeft
ELN = EL - DeltaW

bL = EL + sqrt(DeltaW * (EL-EL1))
Nfrac = (bL - EL)/DeltaW
NL = ceil(Nfrac)

EL1 = EL - (DeltaW*NL)^2 / DeltaW
bL  = EL + sqrt(DeltaW * (EL-EL1))

duL = - DeltaW / ((EL-bL)*(EL-DeltaW-bL))  #<0

##right
ER = EwinInnerRight
ER1 = EwinOuterRight
ERN = ER + DeltaW

bR = ER - sqrt(DeltaW*(ER1-ER))
Nfrac  = (ER-bR)/DeltaW
NR = ceil(Nfrac)

ER1 = ER + (DeltaW*NR)^2 / DeltaW 
bR = ER - sqrt(DeltaW*(ER1-ER)) # = ER-DeltaW*NR

duR = DeltaW/((ER-bR)*(ER+DeltaW-bR))  #>0


#write ERealAxis
for i=1:NL
  push!(ERealAxis,  1/(i*duL) +bL  )  
end

for i=1:EgridInner+1
  push!(ERealAxis,  EwinInnerLeft + (i-1) * DeltaW)   #(i=1)=> EwinInnerLeft, ...,  i=EgridInner+1=> EwinInnerRight
  if(abs(ERealAxis[end]) <= coreWindow  &&  abs(ERealAxis[end]+DeltaW) <= coreWindow )   
    for j=1:coreDense-1
      push!(ERealAxis,  ERealAxis[end]+coreDelta)
    end
  end
end

for i=NR:-1:1
  push!(ERealAxis,  1/(i*duR) + bR)   
end

Egrid = length(ERealAxis)
println("Egrid: ",Egrid)   #Egrid = num of energy point including endpoint , not space number

push!(dERealAxis, (ERealAxis[2] - ERealAxis[1])/2 )
for i=2:Egrid-1;
   push!(dERealAxis, (ERealAxis[i+1] - ERealAxis[i-1])/2 )
end
push!(dERealAxis,  (ERealAxis[Egrid] - ERealAxis[Egrid-1])/2)


return ERealAxis, dERealAxis, Egrid;
end



function read_matsubara_GreenFtn!( inverse_temp::Float64, numeric, Energy_weight::Array{Float64}, startOrbit::Int64, NumSubOrbit::Int64, fname::String)
N_Matsubara = numeric.N_MatsubaraDATA

GreenFtn = Array{Array{Complex128,2}}(N_Matsubara)
GreenConstFull =  zeros(Complex128, NumFullOrbit,NumFullOrbit)
for jj = 1:size(GreenFtn)[1]
  GreenFtn[jj]= zeros(Complex128, NumSubOrbit,NumSubOrbit)
end

#=
#Read Matsubara Green's Ftn
=#
#A = readdlm(fname,' ');
#A = readdlm(fname,',');
A = readdlm(fname);

B=convert(Array{Int64,2},A[:,1:3]);
C = A[:,4:end];
A = [];
iw=0;
iwmax=0;


for jj = 1:size(B)[1]
    iw=B[jj,1];
     i=B[jj,2]+1;
     j=B[jj,3]+1;
    if startOrbit<i<=NumSubOrbit+startOrbit && startOrbit<j<=NumSubOrbit+startOrbit
        isub= i-startOrbit
        jsub= j-startOrbit
        if  0<=iw<N_Matsubara
            GreenFtn[iw+1][isub,jsub]=C[jj,1]+C[jj,2]im
       end
    end
    if  iw>=N_Matsubara-5 && i<=NumFullOrbit && j<=NumFullOrbit
         GreenConstFull[i,j] += C[jj,1]+C[jj,2]im
    end
    if iwmax < iw
            iwmax = iw
    end
end #file lines

GreenConstFull = GreenConstFull/5.0

if(iwmax+1!=N_Matsubara)
    println("Please check N_Matsubara, $(iwmax+1)")
end


#find moments

NumMom=4
momentstest1 = Array{Array{Complex128,2}}(0)
momentstest2 = Array{Array{Complex128,2}}(0)
momentstest3 = Array{Array{Complex128,2}}(0)
momentstest4 = Array{Array{Complex128,2}}(0)
for w_Asymto = 11:N_Matsubara-10
  ASlen = N_Matsubara-w_Asymto   # num of Mat. freq. in the asymtotic region
  KM =  zeros(Float64,ASlen*2,NumMom)
  Gasym =  Array{Array{Complex128,2}}(ASlen*2)
  for jj = 1:ASlen
    wn = N_Matsubara-ASlen+jj
    z = ((2.*wn-1)*pi/inverse_temp);
    KM[2*jj-1,1] =   1.0        #const
    KM[2*jj,  2] =  -1/z        #norm
    KM[2*jj-1,3] =  -1/z^2      #<w>
    KM[2*jj,  4] =   1/z^3      #<w^2>
  
    Gasym[2*jj-1] = Hermitian((GreenFtn[wn] + GreenFtn[wn]')/2);
    Gasym[2*jj  ] = Hermitian((GreenFtn[wn] - GreenFtn[wn]')/(2im));
  
  end
  KMsvdf = svdfact(KM, thin=false)
  moments = (  ((KMsvdf[:Vt]')[:,1:NumMom]  * diagm(1./ KMsvdf[:S]) * (KMsvdf[:U]')[1:NumMom,:])   * Gasym   )
  push!(momentstest1 , Hermitian(moments[3])^(1/1))
  push!(momentstest2 , Hermitian(moments[4])^(1/2))
end
x=Array{Float64}(0)
for w_Asymto = 12:N_Matsubara-11
push!(x,  + norm( momentstest2[w_Asymto-10+1] - momentstest2[w_Asymto-10] )^2 + norm( momentstest2[w_Asymto-10-1] - momentstest2[w_Asymto-10] )^2
     )
end
w_Asymto = indmin(x) + 11
println("Asymto iwn = $((2*w_Asymto-1)*pi/inverse_temp) at n=$(w_Asymto)")

ASlen = N_Matsubara-w_Asymto
KM =  zeros(Float64,ASlen*2,NumMom)
Gasym =  Array{Array{Complex128,2}}(ASlen*2)
for jj = 1:ASlen
  wn = N_Matsubara-ASlen+jj
  z = ((2.*wn-1)*pi/inverse_temp);
  KM[2*jj-1,1] =   1.0
  KM[2*jj,  2] =  -1/z
  KM[2*jj-1,3] =  -1/z^2
  KM[2*jj,  4] =   1/z^3

  Gasym[2*jj-1] = Hermitian((GreenFtn[wn] + GreenFtn[wn]')/2);
  Gasym[2*jj  ] = Hermitian((GreenFtn[wn] - GreenFtn[wn]')/(2im));

end
KMsvdf = svdfact(KM, thin=false)
moments = (  ((KMsvdf[:Vt]')[:,1:NumMom]  * diagm(1./ KMsvdf[:S]) * (KMsvdf[:U]')[1:NumMom,:])   * Gasym   )
moments[1] = Hermitian(moments[1])
moments[2] = Hermitian(moments[2])
moments[3] = Hermitian(moments[3])
moments[4] = Hermitian(moments[4])



println("const   :  ",      real(trace(moments[1])))
println("1/(iw)  :  ",      real(trace(moments[2])))
println("1/(iw)^2:  ",      real(trace(moments[3])))
println("1/(iw)^3:  ",     (real(trace(moments[4]))))

N = real(trace(moments[2]))
c = real(trace(moments[3]))/N
s = real(trace(moments[4]))/N

s= sqrt(s-c^2)
println("center:", c)
println("stddev:", s)



#set N_Matsubara for low freq. 
if N_Matsubara > 500
  temp=[]
  for w = 1:N_Matsubara
      z = ((2.*w-1)*pi/inverse_temp)im;
      push!(temp,
      (norm( GreenFtn[w]  - (moments[1] + moments[2]/z + moments[3]/(z^2) + moments[4]/(z^3) ))    )/(norm(GreenFtn[w]))
      )
  end
  numeric.N_Matsubara = indmin(temp)
  println("Asymto iwn = $((2*numeric.N_Matsubara-1)*pi/inverse_temp) at n=$(numeric.N_Matsubara)")
  GreenFtn = GreenFtn[1:numeric.N_Matsubara]
else
  numeric.N_Matsubara = N_Matsubara
end



#Best unitary transformation
global bestUnitary= eye(NumSubOrbit);
#if basisTransform
#  eigValue, eigVect = eig(Hermitian(moments[2]));
#  bestUnitary = deepcopy(eigVect)
#end

#for iw = 1:numeric.N_Matsubara
# GreenFtn[iw] = copy(bestUnitary') *  GreenFtn[iw] * bestUnitary
#end
#moments = [ copy(bestUnitary')* moments[j] *bestUnitary for j=1:length(moments)]


if !(offDiagonal)
     for iw = 1:numeric.N_Matsubara
     GreenFtn[iw] = diagm(diag(GreenFtn[iw]))
     end
     moments = [diagm(diag(moments[j])) for j=1:length(moments)]
end
moments = [ (moments[j] + moments[j]')/2 for j=1:length(moments)]



#Renomalize...
Normalization = real(trace(moments[2]))
println("Normalization:$Normalization")
println("read Giw[1]   : $(-real((moments[2][1,1])))/x +$(real((moments[4][1,1])))/x**3")
println("read Giw_trace: $(-real(trace(moments[2])))/x +$(real(trace(moments[4])))/x**3")

for iw = 1:numeric.N_Matsubara
        z = ((2.*iw-1)*pi/inverse_temp)im;
        GreenFtn[iw] -= moments[1]
        GreenFtn[iw] /= Normalization
        GreenFtn[iw] *=  sqrt(Energy_weight[iw])
        if inversion_method
             GreenFtn[iw] =  inv(z*eye(GreenFtn[iw]) - GreenFtn[iw]) 
             GreenFtn[iw] /= trace(eye(GreenFtn[iw]))
        end
end

#println("<1>  :",Hermitian(moments[2])/Normalization)
#println("<w>  :",Hermitian(moments[3])/Normalization)
#println("<w^2>:",Hermitian(moments[4])/Normalization)
#println("<stddev>:", eig(Hermitian(moments[4])/Normalization - (Hermitian(moments[3])/Normalization)^2)[1])

temp=strImagFreqFtn(
     GreenFtn
    ,GreenConstFull
    ,Hermitian(moments[1])
    ,Normalization
     ,bestUnitary
    ,Hermitian(moments[2])/Normalization
    ,Hermitian(moments[3])/Normalization
    ,Hermitian(moments[4])/Normalization
    ,zeros(moments[1])
    ,zeros(moments[1])
)

return temp 
end





function KK_relation( Aw::Array{Array{Complex128,2}},  numeric::strNumeric ) 
  NumSubOrbit = size(Aw[1])[1]
  delta = 1e-10
  Aw_RealPart = Array{Array{Complex128,2}}(numeric.Egrid)
  for jj = 1:size(Aw_RealPart)[1]
    Aw_RealPart[jj]= zeros(Complex128, NumSubOrbit,NumSubOrbit)
  end
  for w1 = 1:numeric.Egrid
      omega = numeric.ERealAxis[w1];
      for w2= 1:numeric.Egrid-1

           # Note : 1/(w-w'+im*eta) = P 1/(w-w') -i*pi*deltaFtn
           # ==> P/(w-w') = Real( 1/(w-w'+im*eta) ) 
            Ea = numeric.ERealAxis[w2]
            Eb = numeric.ERealAxis[w2+1]  
            f2a = -pi*Aw[w2]  
            f2b = -pi*Aw[w2+1]
            a = (f2b-f2a)/(Eb-Ea)
            b = (Eb*f2a - Ea*f2b)/(Eb-Ea)
            Aw_RealPart[w1] += a*(Eb + (omega-im*delta)*log(Eb-omega+im*delta)) + b*log(Eb-omega+im*delta)
            Aw_RealPart[w1] -= a*(Ea + (omega-im*delta)*log(Ea-omega+im*delta)) + b*log(Ea-omega+im*delta)

      end
      Aw_RealPart[w1] += Aw_RealPart[w1]'
      Aw_RealPart[w1] /= 2
  end
  Aw_RealPart/=pi
  return Aw_RealPart
end

#############################################################################################
# (end) define Ftns
#############################################################################################



function construct_Kernel_inCubicSpline( numeric::strNumeric)

   N_Matsubara = numeric.N_Matsubara
   
   kernel = strKernel(
    auxiliary_inverse_temp_range
   ,zeros(Complex128, N_Matsubara,   numeric.Egrid)
   ,zeros(Complex128, numeric.Egrid, N_Matsubara)
   ,zeros(Float64, 5, numeric.Egrid)
   ,zeros(Complex128, numeric.Egrid, numeric.Egrid)
   ,zeros(Complex128, 4*(numeric.Egrid-1),numeric.Egrid)
   ,zeros(Complex128, 5,5)
   )


   num_seg_coeff = size(kernel.cubic_coeff_transf)[1]
   dSegment=[] 
   for j=1:numeric.Egrid-1
     push!(dSegment, numeric.ERealAxis[j+1] - numeric.ERealAxis[j]) 
   end
   T = zeros(Float64, num_seg_coeff,  numeric.Egrid)
   B = zeros(Float64, num_seg_coeff, num_seg_coeff)
   Kernel_from_cubic= zeros(Complex128, N_Matsubara, num_seg_coeff)



   #right hand side of the Eq. A2~A5
   for j=1:numeric.Egrid-1
     T[4*(j-1)+1, j  ] =1.0 
     T[4*(j-1)+2, j+1] =1.0 

     T[4*(j-1)+3, j  ] =0.0
     T[4*(j-1)+4, j  ] =0.0
     T[4*(j-1)+3, j+1] =0.0
     T[4*(j-1)+4, j+1] =0.0
   end


   #left hand side of the Eq. A2~A5
   for j=1:numeric.Egrid-2
       s = 4*(j-1)
       B[s+1, s+4] = 1

       B[s+2, s+1] = (dSegment[j])^3
       B[s+2, s+2] = (dSegment[j])^2
       B[s+2, s+3] = (dSegment[j])^1
       B[s+2, s+4] = 1

       B[s+3, s+1] = 3*(dSegment[j])^2
       B[s+3, s+2] = 2*(dSegment[j])^1
       B[s+3, s+3] = 1
       B[s+3, s+7] = -1

       B[s+4, s+1] = 6*(dSegment[j])^1
       B[s+4, s+2] = 2
       B[s+4, s+6] = -2
   end
   j=numeric.Egrid-1
   s = 4*(j-1)
   B[s+1, s+4] = 1

   B[s+2, s+1] = (dSegment[j])^3
   B[s+2, s+2] = (dSegment[j])^2
   B[s+2, s+3] = (dSegment[j])^1
   B[s+2, s+4] = 1




#   #boundary conditions:
#   B[s+3, 3] = 1    #S'(w1) = 0
#   #w_N * S'(w_N) = 0
#   B[s+4, s+1] = ERealAxis[j+1] *3*(dSegment[j])^2
#   B[s+4, s+2] = ERealAxis[j+1] *2*(dSegment[j])^1
#   B[s+4, s+3] = ERealAxis[j+1] *1



   #boundary conditions:
   B[s+3, 2] = 2    #S''(w1) = 0
   #S''(w_N) = 0
   B[s+4, s+1] = 6*dSegment[j]
   B[s+4, s+2] = 2




    #integration of the cubic polynomials
   f1(iwn, wj, w ) =  -(iwn-wj)^2 * w - (iwn-wj)*w^2/2 -w^3/3 -(iwn-wj)^3 * log(-iwn + wj +w)
   f2(iwn, wj, w ) =  -(iwn-wj)*w - w^2 /2                    -(iwn-wj)^2 * log(-iwn + wj +w)
   f3(iwn, wj, w ) =  -w                                      -(iwn-wj)   * log(-iwn + wj +w)
   f4(iwn, wj, w ) =                                          -             log(-iwn + wj +w)

   for n=1:N_Matsubara
    iwn = ((2.*n-1)*pi/inverse_temp)im;
    for j=1:numeric.Egrid-1
       wj = numeric.ERealAxis[j] 
       Kernel_from_cubic[n,4*(j-1)+1] = f1(iwn, wj, dSegment[j]) - f1(iwn, wj, 0)
       Kernel_from_cubic[n,4*(j-1)+2] = f2(iwn, wj, dSegment[j]) - f2(iwn, wj, 0)
       Kernel_from_cubic[n,4*(j-1)+3] = f3(iwn, wj, dSegment[j]) - f3(iwn, wj, 0)
       Kernel_from_cubic[n,4*(j-1)+4] = f4(iwn, wj, dSegment[j]) - f4(iwn, wj, 0)
    end
   end

   kernel.cubic_coeff_transf[:,:] = (inv(B) * T)[:,:]
   kernel.Kernel[:,:] = (Kernel_from_cubic*  inv(B) * T)


   Gamma1 = 0.0
   Gamma2 = 0.0
   Gamma3 = 0.0
   Gamma4 = 0.0
   Gamma5 = 0.0
   
   for iw = N_Matsubara:1e4*N_Matsubara
    z = ((2.*iw-1)*pi/inverse_temp);
    Gamma1 +=      (1/z  )^2
    Gamma2 +=      (1/z^2)^2
    Gamma3 +=      (1/z^3)^2
   end
   
   
   println("Gamma1:",Gamma1)
   println("Gamma2:",Gamma2)
   println("Gamma3:",Gamma3)

   
   kernel.gamma[1,1] = Gamma1 
   kernel.gamma[2,2] = Gamma2 
   kernel.gamma[3,3] = Gamma3 

   kernel.gamma[4,4] = Gamma4
   kernel.gamma[5,5] = Gamma5 
   kernel.gamma[1,3] = Gamma2 
   kernel.gamma[3,1] = Gamma2 


    kernel.gamma = Hermitian(kernel.gamma)
    
   
   
   
#moment[j.:] Aw[:] = \int w^(j-1) A(w).,               
#A(w) \approx (Aw[l+1]-Aw[l])/dSegment[l]*(w-E[l])+Aw[l] 
#        =   Aw[l+1]*(w-E[l])/dSegment[l] - Aw[l]*[(w -E[l])/dSegment[l]-1];
#        =   Aw[l+1]*(w-E[l])/dSegment[l] - Aw[l]*[(w-{E[l]+dSegment[l]}]/dSegment[l];
#        =   Aw[l+1]*(w-E[l])/dSegment[l] - Aw[l]*[(w- E[l+1]]/dSegment[l];

   for l=1:numeric.Egrid-1
    for j=1:5
      b =  ( numeric.ERealAxis[l+1]^(j+1) - numeric.ERealAxis[l]^(j+1) )/(j+1) - numeric.ERealAxis[l]   *  ( numeric.ERealAxis[l+1]^(j) - numeric.ERealAxis[l]^(j) )/(j) 
      c = -( numeric.ERealAxis[l+1]^(j+1) - numeric.ERealAxis[l]^(j+1) )/(j+1) + numeric.ERealAxis[l+1] *  ( numeric.ERealAxis[l+1]^(j) - numeric.ERealAxis[l]^(j) )/(j) 
      kernel.moment[j,l+1] +=   (b/dSegment[l])
      kernel.moment[j,l  ] +=   (c/dSegment[l])

    end
   end
   
   
   
   for l=1:numeric.Egrid
     for n=1:N_Matsubara
       iwn = ((2.*n-1)*pi/inverse_temp)im;
       kernel.Kernel_dagger[l,n] = 1/(-iwn - numeric.ERealAxis[l] )
     end
   end
#   for n=1:N_Matsubara
#      kernel.Kernel[n,:] *=  sqrt(Energy_weight[n])
#   end


   return kernel
end







# #=
function construct_Kernel_linear!( numeric::strNumeric)

   N_Matsubara = numeric.N_Matsubara
   ERealAxis = numeric.ERealAxis
  dERealAxis = numeric.dERealAxis


   kernel = strKernel(
    auxiliary_inverse_temp_range
   ,zeros(Complex128, N_Matsubara, numeric.Egrid)
   ,zeros(Float64, 5, numeric.Egrid)
   ,zeros(Complex128, numeric.Egrid, numeric.Egrid)
   ,zeros(Complex128, 4*(numeric.Egrid-1),numeric.Egrid)
   ,zeros(Complex128, 5,5)
   )




# #=
   num_seg_coeff = 2*(numeric.Egrid-1)
   dSegment=[] 
   for j=1:numeric.Egrid-1
     push!(dSegment, ERealAxis[j+1] - ERealAxis[j]) 
   end
   T = zeros(Float64, num_seg_coeff,  numeric.Egrid)
   B = zeros(Float64, num_seg_coeff, num_seg_coeff)
   Kernel_from_cubic= zeros(Complex128, N_Matsubara, num_seg_coeff)


   for j=1:numeric.Egrid-1
     T[2*(j-1)+1, j]   =1.0
     T[2*(j-1)+2, j+1] =1.0
   end

   for j=1:numeric.Egrid-1
       s = 2*(j-1)
       B[s+1, s+2] = 1

       B[s+2, s+1] = (dSegment[j])^1
       B[s+2, s+2] = 1
   end


   f3(iwn, wj, w ) =  -w                                      -(iwn-wj)   * log(-iwn + wj +w)
   f4(iwn, wj, w ) =                                          -             log(-iwn + wj +w)

   for n=1:N_Matsubara
    iwn = ((2.*n-1)*pi/inverse_temp)im;
    for j=1:numeric.Egrid-1
       wj = ERealAxis[j] 
       Kernel_from_cubic[n,2*(j-1)+1] = f3(iwn, wj, dSegment[j]) - f3(iwn, wj, 0)
       Kernel_from_cubic[n,2*(j-1)+2] = f4(iwn, wj, dSegment[j]) - f4(iwn, wj, 0)
    end
   end
# =#

kernel.Kernel[:,:] = Kernel_from_cubic*  inv(B) * T


 #=
   for iw = 1:N_Matsubara
       z = ((2.*iw-1)*pi/inverse_temp)im;
       for w= 1:numeric.Egrid
           if w==1
            El   =  ERealAxis[w]
           else
            E_m1 =  ERealAxis[w-1]
            El   = (ERealAxis[w] + E_m1 )/2
           end
           E0   =  ERealAxis[w]
           if w==numeric.Egrid
             Er   =  ERealAxis[w]
           else
             E_p1 =  ERealAxis[w+1]
             Er   = (ERealAxis[w] + E_p1 )/2   #Here, Er-El = dERealAxis[i] = (ERealAxis[i+1] - ERealAxis[i-1])/2
           end
           kernel.Kernel[iw,w]    +=               log((El-z)/(Er-z))                     / dERealAxis[w]
           #left
           if(w==1)
#             Kernel[iw,w]  +=   ((E0^2/z^2 ) * log((E0)/((E0)-z)) - E0^2/z* (1/(E0))) /dERealAxis[w]
           else
             kernel.Kernel[iw,w-1] += ( (El-E0) + (z-E0)*log((El-z)/(E0-z)) )/ ( E_m1- E0) / dERealAxis[w-1]
             kernel.Kernel[iw,w]   -= ( (El-E0) + (z-E0)*log((El-z)/(E0-z)) )/ ( E_m1- E0) / dERealAxis[w]
           end
           #right
           if(w==numeric.Egrid)
#             Kernel[iw,w] +=   (-(E0^2/z^2 ) * log((E0)/(E0-z)) + E0^2/z* (1/E0))   /dERealAxis[w]
           else
             kernel.Kernel[iw,w+1] += ( (E0-Er) + (z-E0)*log((E0-z)/(Er-z)))/ ( E_p1- E0) / dERealAxis[w+1]
             kernel.Kernel[iw,w]   -= ( (E0-Er) + (z-E0)*log((E0-z)/(Er-z)))/ ( E_p1- E0) / dERealAxis[w]
           end
       end
      kernel.Kernel[iw,:] *=  sqrt(Energy_weight[iw])
   end
# =#
end
# =#
