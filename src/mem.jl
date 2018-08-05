# fit to find optimal alpha
# https://github.com/JuliaNLSolvers/LsqFit.jl
#Pkg.clone("https://github.com/wildart/TOML.jl.git")
###########################################
# Jae-Hoon Sim, KAIST 2018.01.
#2018.07.13
#2018.08.03
###########################################



using LsqFit
import TOML



using MQEM


#input  param
input_toml = TOML.parsefile(ARGS[1])
function input_ftn(var::String)
  if haskey(input_toml, var) 
     println("Input : $(var):",input_toml[var] )
    return input_toml[var]
  else
     println("Input error: ",var )
     exit()
  end
end
function input_ftn(var::String, val)
  if haskey(input_toml, var) 
     println("Input : $(var):",input_toml[var] )
    return input_toml[var]
  else
     println("Input default: $(var):",val )
    return val
  end
end
####################################################
workDirect      =  input_ftn("workDirect", pwd())
inputFile       =  input_ftn("inputFile")
NumFullOrbit        = input_ftn("NumOrbit")
num_of_subblock =  input_ftn("num_of_subblock", 1)
start_cluster    = input_ftn("start_cluster",0)
N_Matsubara2     = input_ftn("N_Matsubara")
inverse_temp    = input_ftn("inverse_temp")
Asymtotic_HighFreq  = input_ftn("Asymtotic_HighFreq", true)



data_info = data_info_(workDirect,
		       inputFile,
		       NumFullOrbit,
		       num_of_subblock,
		       start_cluster,
		       N_Matsubara2,
		       inverse_temp,
		       Asymtotic_HighFreq)

EwinOuterRight  = input_ftn("EwinOuterRight", 10.0)
EwinOuterLeft   = input_ftn("EwinOuterLeft" ,-10.0)
EwinInnerRight  = input_ftn("EwinInnerRight" , 3.0)
EwinInnerLeft   = input_ftn("EwinInnerLeft"  ,-3.0)
coreWindow      = input_ftn("coreWindow"  , 0.0)
coreDense       = input_ftn("coreDense"   , 2)
EgridInner           = input_ftn("Egrid", 100)

real_freq_grid_info= real_freq_grid_info_(EwinOuterRight,
					   EwinOuterLeft,
					   EwinInnerRight,
					   EwinInnerLeft,
					   coreWindow,
					   coreDense,
					   EgridInner
					   )

default_model  = input_ftn("default_model" ,"g_mat")
Model_range_right = input_ftn("Model_right",EwinOuterRight)
Model_range_left  = input_ftn("Model_left" ,EwinOuterLeft)
auxiliary_inverse_temp_range = input_ftn("auxiliary_inverse_temp_range")
auxTempRenomalFactorInitial  = input_ftn("auxTempRenomalFactorInitial", 0.01)

mem_fit_parm= mem_fit_parm_(default_model,
			    Model_range_right,
			    Model_range_left,
			     auxiliary_inverse_temp_range,
			     auxTempRenomalFactorInitial
			     )



NumIter              = input_ftn("NumIter",       10000)
mixingInitial        = input_ftn("mixingInitial", 0.1)
mixing_max           = input_ftn("mixing_max", 0.3)
mixing_min           = input_ftn("mixing_min", 1e-6)
pulay_mixing_start   = input_ftn("pulay_mixing_start", 100)
pulay_mixing_step    = input_ftn("pulay_mixing_step", 1)
pulay_mixing_history = input_ftn("pulay_mixing_history", 1)

mixing_parm= mixing_parm_(NumIter,
			   mixingInitial,
			   mixing_max,
			   mixing_min,
			   pulay_mixing_start,
			   pulay_mixing_step,
			   pulay_mixing_history
			   )

inputInfo = strInputInfo(data_info,
			 real_freq_grid_info,
			 mem_fit_parm,
			 mixing_parm
			 )

#############################################################################################
# (Start) Initialize real space variables 
#############################################################################################



#ERealAxis, dERealAxis, Egrid2 = construct_EnergyGrid(EwinOuterLeft, EwinOuterRight, EwinInnerLeft, EwinInnerRight, EgridInner)
ERealAxis, dERealAxis, Egrid2 = construct_EnergyGrid(inputInfo.real_freq_grid_info)
println(ERealAxis);

numeric = strNumeric(
  Egrid2
 ,ERealAxis
 ,dERealAxis
 ,N_Matsubara2
 ,N_Matsubara2
 ,NumIter
 ,mixing_min
 ,mixing_max
 ,mixingInitial
 ,pulay_mixing_history
 ,pulay_mixing_step
 ,pulay_mixing_start
)


phyParm = strPhyParm(
           inverse_temp
          ,NumFullOrbit
          ,div(NumFullOrbit,num_of_subblock)::Int64
          ,num_of_subblock
          )


realFreqFtn  = strRealFreqFtn(
  Array{Array{Complex128,2}}(numeric.Egrid)
 ,Array{Array{Complex128,2}}(numeric.Egrid)
 ,Array{Array{Complex128,2}}(numeric.Egrid)
 ,Array{Array{Complex128,2}}(numeric.Egrid)
 ,Array{Array{Complex128,2}}(numeric.Egrid)
)
for w=1:numeric.Egrid
   realFreqFtn.H_extern[w] = Hermitian(zeros(Complex128, phyParm.NumSubOrbit, phyParm.NumSubOrbit))
end



#############################################################################################
# (end) Initialize variables 
#############################################################################################

#etc...
#Energy_weight = zeros(N_Matsubara2)
#for jj = 1:N_Matsubara2
#  wval = ((2.*jj-1)*pi/inverse_temp)
#  Energy_weight[jj] =  1.0
#end

println("sub-space dim=$(phyParm.NumSubOrbit)")




#############################################################################################
# #=  mainloop
#############################################################################################





#Initialize matsubara information

####start main
f=open("$(workDirect)/information.out","w");write(f,"\n" );close(f)
for cluster=start_cluster:num_of_subblock-1
    global startOrbit = phyParm.NumSubOrbit*cluster
    println("Block = $(startOrbit)")
    
    # =#
    #file
#    global fname = "$(workDirect)/$(inputFile)"
    global fname_out = Array{String}(phyParm.NumSubOrbit,phyParm.NumSubOrbit)
    global fname_contniuedSpectrum = Array{String}(phyParm.NumSubOrbit,phyParm.NumSubOrbit)
    global fname_reproduce = Array{String}(phyParm.NumSubOrbit,phyParm.NumSubOrbit)
    for i = 0:phyParm.NumSubOrbit-1
        for j=0:phyParm.NumSubOrbit-1
            fname_out[i+1,j+1] = "$(workDirect)/Density_of_state_$(i+startOrbit)_$(j+startOrbit).dat"
      fname_reproduce[i+1,j+1] = "$(workDirect)/reproduce_$(i+startOrbit)_$(j+startOrbit).out"
      fname_contniuedSpectrum[i+1,j+1] = "$(workDirect)/realFreq_Sw.dat_$(i+startOrbit+1)_$(j+startOrbit+1)"
        end
    end
    
    
    GreenConstFull =  zeros(Complex128, NumFullOrbit,NumFullOrbit)
    (imagFreqFtn) = read_matsubara_GreenFtn!( data_info,  numeric, startOrbit, phyParm.NumSubOrbit)
    kernel = construct_Kernel_inCubicSpline(numeric, data_info)
    sigmaFlat = (EwinInnerRight - EwinInnerLeft)
    
    sigma = min(( trace(imagFreqFtn.moments3)-trace(imagFreqFtn.moments2)^2),  ((EwinOuterRight - EwinOuterLeft)/4)^2)
    if sigma>0
       sigma = sqrt(sigma)
    else 
      default_model = "f"
    end
    
    
    if(default_model=="f") #flat default_model
	    realFreqFtn.Spectral_default_model = 
	    construct_model_spectrum(imagFreqFtn.moments1, numeric,inputInfo.mem_fit_parm,  phyParm.NumSubOrbit, 
				     trace(imagFreqFtn.moments2), sigmaFlat,"F", kernel) 

    elseif(default_model =="g")
	    realFreqFtn.Spectral_default_model =
	    construct_model_spectrum(imagFreqFtn.moments1, numeric,inputInfo.mem_fit_parm,  phyParm.NumSubOrbit,
				     trace(imagFreqFtn.moments2), trace(imagFreqFtn.moments3)/phyParm.NumSubOrbit,"G", kernel)

    elseif(default_model=="g_mat")  #gaussian default_model
      Ainit= Hermitian(-1.0/(2*sigma^2) *eye(imagFreqFtn.moments2))
      Binit= Hermitian(-( trace(imagFreqFtn.moments2) / sigma^2)*eye(imagFreqFtn.moments2))
      s= 3*phyParm.NumSubOrbit +2*3*(div(phyParm.NumSubOrbit*(phyParm.NumSubOrbit-1),2)) -1     # -1 come from the fact that C=trace less
        # First term = diagonal for A,B,C
        # sencdterm = real, imag foa off-diagoal of the  A,B,C
    
      xinit=zeros(Float64,s)
    
      for i=1:phyParm.NumSubOrbit
           xinit[i]                     = Ainit[i,i]
           xinit[phyParm.NumSubOrbit+i] = Binit[i,i]
      end
     
      xinit = gaussianPotential!(imagFreqFtn, realFreqFtn, kernel, numeric, xinit, realFreqFtn.Spectral_default_model)
      moments_rep  = kernel.moment * realFreqFtn.Spectral_default_model
      println("Model: $(-real(trace(moments_rep[1])))/x +$(real(trace(moments_rep[3])))/x**3")
    end
    
    temp = deepcopy(realFreqFtn.Spectral_default_model)
    write_spectral_ftn(phyParm.NumSubOrbit, imagFreqFtn.Normalization, numeric, temp, fname_out, "_model")
    temp = []
    


 #Main part!########################################################################
 tic()
    f=open("$(workDirect)/information.out","a");write(f,"\nsub_block$(cluster)\n" );close(f)
 # #=
 #  #= get aux. temperature spectrum
   realFreqFtn.Aw = deepcopy(realFreqFtn.Spectral_default_model)
   (logEnergy, logAlpha) = mem_annealing( kernel,  realFreqFtn, imagFreqFtn,  numeric, mem_fit_parm, mixing_parm, true, fname_out, data_info, startOrbit)
 #  =#
 toc()
 




 tic()
 #  #= Find optimal alpha
   function hyp_tan_model(x,p)
    return  -p[1] * tanh.((x-p[2])/p[3]) + p[4]   #note : -tanh(x) = -2/(e^(2x) +1 ) -1,    p[3]= 2*temperature in FD distribution
   end
   p0=[ 1 , mean(logAlpha), (logAlpha[end] -logAlpha[1])/2 , mean(logEnergy)  ]
   fit = curve_fit( hyp_tan_model, logAlpha,  logEnergy,  p0  )
   p= fit.param

 #  alpha_optimal =  exp(p[3] *log((1+sqrt(5))/2) + p[2])   #second derivative of tanh(x) has maximum at log((1+sqrt(5))/2)
   alpha_optimal =  exp( 1.5 *p[3] + p[2])
   println("fitted logE(loga): -$(p[1])*tanh((x-($(p[2])))/$(p[3])) + $(p[4])")


 #  #= find optimal spectrum
   if alpha_optimal > mem_fit_parm.auxiliary_inverse_temp_range[2]   alpha_optimal=deepcopy(mem_fit_parm.auxiliary_inverse_temp_range[2]) 
   else
      mem_fit_parm.auxiliary_inverse_temp_range[2] = alpha_optimal
      realFreqFtn.Aw = deepcopy(realFreqFtn.Spectral_default_model)
      (logEnergy, logAlpha) = mem_annealing( kernel,  realFreqFtn, imagFreqFtn,numeric,   mem_fit_parm, mixing_parm, false, fname_out, data_info, startOrbit)
   end
   energy_optimal = get_total_energy( imagFreqFtn, realFreqFtn.Aw, kernel,numeric);
   energy_tail_optimal = get_total_energy_for_tail( imagFreqFtn, realFreqFtn.Aw, kernel,numeric);
   println("optimal alpha: $(alpha_optimal),  log: $(log(alpha_optimal))")
   println("Total_energy_at_optimal_alpha: $(energy_optimal) tail_contribution: $(energy_tail_optimal) with_grid $(numeric.Egrid) \n")
 #  =#
 toc()
    
    
    
    
    ##################################################################################
    #KK relation
    Aw_RealPart = KK_relation( realFreqFtn.Aw, numeric)
    write_results(phyParm.NumSubOrbit, fname_out, fname_contniuedSpectrum, realFreqFtn.Aw, numeric.ERealAxis, imagFreqFtn.Normalization, imagFreqFtn.GreenConst, Aw_RealPart, numeric.Egrid)
    
    for Exchangecluster=cluster+1:num_of_subblock-1
      for i = 1:phyParm.NumSubOrbit
          for j=1:phyParm.NumSubOrbit
              ifull=i+startOrbit;
              jfull= j+ (phyParm.NumSubOrbit*Exchangecluster)
              f=open("$(workDirect)/realFreq_Sw.dat_$(ifull)_$(jfull)","w")
              for w=1:numeric.Egrid
                  E = numeric.ERealAxis[w]
                  Ftnij =  imagFreqFtn.GreenConstFull[ifull,jfull]
                  write(f, "$E $(real(Ftnij)) $(imag(Ftnij))\n")
              end
              close(f)
              f=open("$(workDirect)/realFreq_Sw.dat_$(jfull)_$(ifull)","w")
              for w=1:numeric.Egrid
                  E = numeric.ERealAxis[w]
                  Ftnji = imagFreqFtn.GreenConstFull[jfull,ifull]
                  write(f, "$E $(real(Ftnji)) $(imag(Ftnji))\n")
              end
              close(f)
          end
       end
    end

##=  mainloop end 
end
#############################################################################################
