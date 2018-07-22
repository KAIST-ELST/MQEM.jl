###########################################
# Jae-Hoon Sim, KAIST 2018.01.
###########################################
#############################################################################################
# (start) define Ftns
#############################################################################################

#include("mem_variable.jl")
#include("mem_gaussianPot.jl")

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
                    normVector_temp[w] =  (norm(Aw[w]))^2 *weight[w]
              end
              return sqrt( sum(normVector_temp))
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
                                     
	  Gmoment = [ imagFreqFtn.moments1, imagFreqFtn.moments2, imagFreqFtn.moments3, imagFreqFtn.moments4, imagFreqFtn.moments5 ]
	  Hamiltonian_out  = -auxiliary_inverse_temp * ( (kernel.Kernel_dagger)*(imagFreqFtn.GreenFtn) ) ;    
 	  Hamiltonian_out += -auxiliary_inverse_temp * ( (kernel.moment_dagger)*(kernel.gamma*Gmoment) ) ;    

      # d(G-KA)^2 =  -K'(G - KA) , second term,i.e. K'K,  = interaction term, 
  
      for w=1:Egrid
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
  
      A_transpose = zeros(Complex128, Egrid);
      for i=1:NumSubOrbit*NumSubOrbit
          for w=1:Egrid
              A_transpose[w] = Aw[w][i]
          end
          h_updateTran  = kernel.Kernel*A_transpose
          h_updateTran  = auxiliary_inverse_temp * (kernel.Kernel_dagger*(h_updateTran))::Array{Complex128}

	  h_update_moment = kernel.gamma * (kernel.moment * A_transpose)
	  h_update_moment = auxiliary_inverse_temp *  (kernel.moment_dagger * h_update_moment)
          for w=1:Egrid
		  Hamiltonian_out[w][i] += h_updateTran[w]  + h_update_moment[w]
          end
      end
  
      for w=1:Egrid
          Hamiltonian_out[w] +=  Hamiltonian_out[w]'
          Hamiltonian_out[w] = Hermitian(Hamiltonian_out[w]/2.0)
      end
  
      return Hamiltonian_out #rw
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
	      Aw_out[w] = Hermitian(expm(Hermitian( -Hamiltonian[w]) ))
	      local_trace_value[w] = trace(Aw_out[w])
          end
          #sum_rule
          TotalPartitionFtn =  kernel.moment[1,:]'local_trace_value
          for w=1:numeric.Egrid
            Aw_out[w] *= 1./TotalPartitionFtn
          end
          return Aw_out
  end
  
  
  
  
  function find_optimal_mixing_weight!(mixing_local::Float64,    criteria_Aw::Float64, criteria_Aw_prev100::Float64, mixing_parm::mixing_parm_)
                                                                      
     if(criteria_Aw < criteria_Aw_prev100)
            mixing_local  *= 1.02
     elseif(criteria_Aw > criteria_Aw_prev100 ) 
            mixing_local  /= 1.05
     end
     mixing_local = min(mixing_local, mixing_parm.mixing_max)
     mixing_local = max(mixing_local, mixing_parm.mixing_min)
     return mixing_local 
  end
  
  
  
  function Aw_Iteration(realFreqFtn::strRealFreqFtn, imagFreqFtn::strImagFreqFtn, kernel::strKernel, auxiliary_inverse_temp::Float64, NumSubOrbit::Int64, numeric::strNumeric, mixing::Float64, mixing_parm::mixing_parm_)
      #Start iteration loop
      MixingInformation = pulayInfo(mixing_parm.pulay_mixing_history ,mixing_parm.pulay_mixing_start,  mixing_parm.pulay_mixing_step,
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
  
          if iter%10 == 0  || iter == 1
              #estimate criteria
              resdAw = Aw_in - Aw_out
              criteria_Aw =  matrixFtn_norm(resdAw, numeric.dERealAxis) ;
  
                if(  criteria_Aw < 1e-11*NumSubOrbit     )
	          E_in = get_total_energy(imagFreqFtn, Aw_in, kernel, numeric)
	          E_out = get_total_energy(imagFreqFtn, Aw_out, kernel, numeric)
		  if( abs(E_in-E_out) < 1e-8 )
			  converg = true
			  realFreqFtn.Hw[:] = Hamiltonian_in
			  break
		  end
                elseif(   ( (iter > numeric.NumIter &&  criteria_Aw > criteria_Aw_prev100) || criteria_Aw > criteria_Aw_prev100 *10)    )
                  converg = false
                  break
                end
              mixing_local = find_optimal_mixing_weight!(mixing_local, criteria_Aw, criteria_Aw_prev100, mixing_parm)
              criteria_Aw_prev100 = criteria_Aw 
          end 
          #Prepare next iteration
          Aw_in = pulayMixing!(iter, mixing_local, Aw_in,   Aw_out-Aw_in, numeric.dERealAxis, MixingInformation)
          iter += 1
      end #end of iteration
      return ( Aw_in  ,   converg,   criteria_Aw, iter)
  end
  
  
  
  
  
  
  
  
  
  function  mem_annealing( kernel::strKernel,  realFreqFtn::strRealFreqFtn, imagFreqFtn::strImagFreqFtn, numeric::strNumeric,   mem_fit_parm::mem_fit_parm_, mixing_parm::mixing_parm_, write_information::Bool, fname_out, data_info::data_info_, startOrbit::Int64)
  
      start_temperature = mem_fit_parm.auxiliary_inverse_temp_range[1]
      end_temperature =   mem_fit_parm.auxiliary_inverse_temp_range[2]
      logEnergy =  Array{Float64}(0)
      logAlpha  =  Array{Float64}(0)
      alpha_list = Array{Float64}(0)
      mixing =  0.0::Float64
      Aw = deepcopy(realFreqFtn.Aw)
      auxTempRenomalFactor = deepcopy(mem_fit_parm.auxTempRenomalFactorInitial)
      NumSubOrbit = size(imagFreqFtn.GreenFtn[1])[1]
      workDirect =  data_info.workDirect
      
      mixing = deepcopy(numeric.mixingInitial)
      auxiliary_inverse_temp = deepcopy(start_temperature);
      
      information_file=open("$(workDirect)/information.out","a")
      write(information_file,"\n" )
      close(information_file)   
      
      auxiliary_inverse_temp_prev = deepcopy(auxiliary_inverse_temp)
      
      
      
      
      mixingPrev = deepcopy(mixing)
      mixingNext =  deepcopy(mixing)
      iter=0;
      while  auxiliary_inverse_temp  <=  end_temperature
          mixing = mixingNext
          auxTempRenomalFactor = min(1+ ((auxTempRenomalFactor-1)*1.01)   , mem_fit_parm.auxTempRenomalFactorInitial)
          auxiliary_inverse_temp *= auxTempRenomalFactor;  
      
      
      
          (Aw,  converg,   normAwRD, iter) = Aw_Iteration(realFreqFtn, imagFreqFtn, kernel, auxiliary_inverse_temp, NumSubOrbit, numeric, mixing, mixing_parm)
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
          (Aw, converg,   normAwRD,iter) = Aw_Iteration(realFreqFtn, imagFreqFtn,kernel, auxiliary_inverse_temp, NumSubOrbit, numeric,mixing, mixing_parm)
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
          
              information_file=open("$(workDirect)/information.out","a")
              write(information_file,"$auxiliary_inverse_temp\t\t$(Energy) ; $(iter)\t$(mixing)\t$(auxTempRenomalFactor)   \n" )
              close(information_file)   
             write_spectral_ftn(NumSubOrbit, imagFreqFtn.Normalization, numeric, Aw, fname_out,  "")
             # comment:Reproduce
              GreenFtn_rep = kernel.Kernel *  Aw ;
              moments_rep  = kernel.moment *  Aw ;
      
              GreenFtn_rep *= imagFreqFtn.Normalization ;
              for i = 1:NumSubOrbit
                  for j=1:NumSubOrbit
                  f=open( "$(workDirect)/reproduce_$(i+startOrbit)_$(j+startOrbit).out","w")
                  for iw=1:numeric.N_Matsubara
                      z = ((2.*iw-1)*pi/data_info.inverse_temp)im;
                      ftn = GreenFtn_rep[iw][i,j] + imagFreqFtn.GreenConst[i,j] 
                      write(f,"$(imag(z))  $(real(ftn))    $(imag(ftn))  \n" )
                  end
                  for iw=numeric.N_Matsubara+1:numeric.N_MatsubaraDATA
                      z = ((2.*iw-1)*pi/data_info.inverse_temp)im;
		      ftn =   moments_rep[1]/z +  moments_rep[2]/ z^2 + moments_rep[3] /z^3+ moments_rep[4] /z^4
		      ftn += imagFreqFtn.GreenConst
		      ftn = ftn[i,j]
                      write(f,"$(imag(z))  $(real(ftn))    $(imag(ftn))  \n" )
                  end
                  close(f)
                  end
              end
      
#      
#              f=open("$(workDirect)/cubic_coeff.out","w")
#               temp = kernel.cubic_coeff_transf * Aw
#               for w=1:numeric.Egrid-1
#                 for subinex = 0:19
#                  E = (20-subinex)/20 * numeric.ERealAxis[w] + subinex/20 * numeric.ERealAxis[w+1]
#                  interpol_val =  temp[4*(w-1)+1][1,1]*(E-numeric.ERealAxis[w])^3 + 
#                                  temp[4*(w-1)+2][1,1]*(E-numeric.ERealAxis[w])^2 + 
#                                  temp[4*(w-1)+3][1,1]*(E-numeric.ERealAxis[w])^1 + 
#                                  temp[4*(w-1)+4][1,1]*(E-numeric.ERealAxis[w])^0
#                  
#                  interpol_val = real(interpol_val)*imagFreqFtn.Normalization
#                  write(f,"$(E) $(interpol_val)\n")  
#                  end
#               end
#               close(f)
          end #write_information
          # =#

      end #while_ (auxiliary_inverse_temp  <=  end_temperature)
  
      return  logEnergy,  logAlpha
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  #function construct_EnergyGrid( EwinOuterLeft::Float64, EwinOuterRight::Float64, EwinInnerLeft::Float64, EwinInnerRight::Float64, EgridInner::Int64)
  function construct_EnergyGrid( real_freq_grid_info::real_freq_grid_info_)
  	EwinInnerRight = real_freq_grid_info.EwinInnerRight
  	EwinInnerLeft = real_freq_grid_info.EwinInnerLeft
  	EwinOuterRight = real_freq_grid_info.EwinOuterRight
  	EwinOuterLeft = real_freq_grid_info.EwinOuterLeft
  	coreDense     = real_freq_grid_info.coreDense
  	coreWindow    = real_freq_grid_info.coreWindow
  	EgridInner    = real_freq_grid_info.EgridInner
  
  
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
  
  
  
  function read_matsubara_GreenFtn!( data_info::data_info_, numeric::strNumeric, startOrbit::Int64, NumSubOrbit::Int64)
  N_Matsubara = numeric.N_MatsubaraDATA
  fname = "$(data_info.workDirect)/$(data_info.inputFile)"
  inverse_temp = data_info.inverse_temp
  NumFullOrbit = data_info.NumFullOrbit
  
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
	  println("Please check N_Matsubara, $(iwmax+1) != $(N_Matsubara)")
  end
  
  
  #find moments
  
  NumMom=4
  moments    = Array{Array{Complex128,2}}(4)
  momentstrace = Array{Array{Float64}}(4)
  momentstest  = Array{Array{Array{Complex128,2}}}(4)
  momentstrace[1] = []
  momentstrace[2] = []
  momentstrace[3] = []
  momentstrace[4] = []
  momentstest[1] = []
  momentstest[2] = []
  momentstest[3] = []
  momentstest[4] = []

  w_Asymto_start=1
  w_Asymto_end=N_Matsubara-10
  w_Asymto_range = w_Asymto_end - w_Asymto_start +1
  for w_Asymto = w_Asymto_start:w_Asymto_end
    ASlen = N_Matsubara-w_Asymto+1   # num of Mat. freq. in the asymtotic region 
    KM =  zeros(Float64,ASlen*2,NumMom)
    Gasym =  Array{Array{Complex128,2}}(ASlen*2)
    for jj = 1:ASlen
	wn = (w_Asymto-1)+jj
	z = ((2.*wn-1)*pi/inverse_temp);
	KM[2*jj-1,1] =   1.0        #const ,
	KM[2*jj,  2] =  -1/z        #norm  , >0
	KM[2*jj-1,3] =  -1/z^2      #<w>
	KM[2*jj,  4] =   1/z^3      #<w^2> , >0
    
	Gasym[2*jj-1] = Hermitian((GreenFtn[wn] + GreenFtn[wn]')/2);
	Gasym[2*jj  ] = Hermitian((GreenFtn[wn] - GreenFtn[wn]')/(2im));

    end
    KMsvdf = svdfact(KM, thin=false)
    moments = (  ((KMsvdf[:Vt]')[:,1:NumMom]  * diagm(1./ KMsvdf[:S]) * (KMsvdf[:U]')[1:NumMom,:])   * Gasym   )
    push!(momentstest[1], Hermitian(moments[1]))   #const
    push!(momentstest[2], Hermitian(moments[2]))   #norm
    push!(momentstest[3], Hermitian(moments[3]))   #<w>
    push!(momentstest[4], Hermitian(moments[4]))   #<w^2>

    push!(momentstrace[1] ,  trace(momentstest[1][end]))  #const
    push!(momentstrace[2] ,  trace(momentstest[2][end]))  #norm
    push!(momentstrace[3] ,  trace(momentstest[3][end]))  #<w>
    push!(momentstrace[4] ,  trace(momentstest[4][end]))  #<w^2>
  end

  Nv =  div(N_Matsubara ,20) #why 20?
  for mom=1:4
     var_j=[]
     for j = Nv+1:w_Asymto_range-Nv
	     push!(var_j , var(momentstrace[mom][j-Nv:j+Nv])  )
     end
     j0 = indmin(var_j) + Nv
     moments[mom] =  mean(momentstest[mom][j0-Nv:j0+Nv])
  end




  println(moments[1])
  println(moments[2])
  println(moments[3])
  println(moments[4])
  
  
  
  
  println("const   :  ",      real(trace(moments[1])))
  println("1/(iw)  :  ",      real(trace(moments[2])))
  println("1/(iw)^2:  ",      real(trace(moments[3])))
  println("1/(iw)^3:  ",     (real(trace(moments[4]))))
  
  N = real(trace(moments[2]))
  c = real(trace(moments[3]))/N
  s = real(trace(moments[4]))/N
  
  s= (s-c^2)
  println("center  :", c)
  println("variance:", s)
  
  
  
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
#  global bestUnitary= eye(NumSubOrbit);
  #if basisTransform
  #  eigValue, eigVect = eig(Hermitian(moments[2]));
  #  bestUnitary = deepcopy(eigVect)
  #end
  
  #for iw = 1:numeric.N_Matsubara
  # GreenFtn[iw] = copy(bestUnitary') *  GreenFtn[iw] * bestUnitary
  #end
  #moments = [ copy(bestUnitary')* moments[j] *bestUnitary for j=1:length(moments)]
  
  
  #if !(offDiagonal)
  #     for iw = 1:numeric.N_Matsubara
  #     GreenFtn[iw] = diagm(diag(GreenFtn[iw]))
  #     end
  #     moments = [diagm(diag(moments[j])) for j=1:length(moments)]
  #end
#  moments = [ (moments[j] + moments[j]')/2 for j=1:length(moments)]
  
  
  
  #Renomalize...
  Normalization = real(trace(moments[2]))
  println("Normalization:$Normalization")
  println("read Giw[1]   : $(-real((moments[2][1,1])))/x +$(real((moments[4][1,1])))/x**3")
  println("read Giw_trace: $(-real(trace(moments[2])))/x +$(real(trace(moments[4])))/x**3")
  
  for iw = 1:numeric.N_Matsubara
          z = ((2.*iw-1)*pi/inverse_temp)im;
          GreenFtn[iw] -= moments[1]
          GreenFtn[iw] /= Normalization
  #        if inversion_method
  #             GreenFtn[iw] =  inv(z*eye(GreenFtn[iw]) - GreenFtn[iw]) 
  #             GreenFtn[iw] /= trace(eye(GreenFtn[iw]))
  #        end
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
  
  
  
  function construct_Kernel_inCubicSpline( numeric::strNumeric, data_info::data_info_ )
  
     N_Matsubara = numeric.N_Matsubara
     
     kernel = strKernel(
      zeros(Complex128, N_Matsubara,   numeric.Egrid)
     ,zeros(Complex128, numeric.Egrid, N_Matsubara)
     ,zeros(Float64, 5, numeric.Egrid)
     ,zeros(Float64, numeric.Egrid, 5)
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
  
  
  
  
     #integration of the cubic polynomials  \int (1/(iwn - w) *  (w-wj)^3 +  (w-wj)^2 +  (w-wj) + 1 ) = f1+f2+f3+f4
     f1(iwn, wj, w ) =  -(iwn-wj)^2 * w - (iwn-wj)*w^2/2 -w^3/3 -(iwn-wj)^3 * log(-iwn + wj +w)
     f2(iwn, wj, w ) =  -(iwn-wj)*w - w^2 /2                    -(iwn-wj)^2 * log(-iwn + wj +w)
     f3(iwn, wj, w ) =  -w                                      -(iwn-wj)   * log(-iwn + wj +w)
     f4(iwn, wj, w ) =                                          -             log(-iwn + wj +w)
  
     for n=1:N_Matsubara
      iwn = ((2.*n-1)*pi/data_info.inverse_temp)im;
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
      z = ((2.*iw-1)*pi/data_info.inverse_temp);
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

#     kernel.gamma[4,4] = Gamma4
#     kernel.gamma[5,5] = Gamma5 
#     kernel.gamma[1,3] = Gamma2 
#     kernel.gamma[3,1] = Gamma2 
  
  
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
         iwn = ((2.*n-1)*pi/data_info.inverse_temp)im;
         kernel.Kernel_dagger[l,n] = 1/(-iwn - numeric.ERealAxis[l] )
       end
     end
     for l=1:numeric.Egrid
      for j=1:5
	      kernel.moment_dagger[l  ,j] =  numeric.ERealAxis[l]^(j-1)
      end
     end
  
  
     return kernel
  end
  
  
  
  
  
  
  
   #=
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
  end
  # =#


#  module default_model
         using NLsolve
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
#                 Aw[w], local_trace_value[w]       = density_matrix_constructor(  )
		 Aw[w] = Hermitian(expm(Hermitian( - (realFreqFtn.H_extern[w] -chempot *eye(realFreqFtn.H_extern[w]))  )))
		 local_trace_value[w] = trace(Aw[w])
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


#         export gaussianPotential!
#end

