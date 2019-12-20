struct strPhyParm
 inverse_temp::Float64
 NumFullOrbit::Int64
 NumSubOrbit::Int64
 num_of_subblock::Int64
end


struct strImagFreqFtn
    GreenFtn
    GreenConstFull
    GreenConst::Array{ComplexF64,2}
    Normalization::Float64
    moments1
    moments2
    moments3
    moments4
    moments5
end


mutable struct strRealFreqFtn
 Spectral_default_model::Array{Array{ComplexF64,2}}
 logSpectral_default_model::Array{Array{ComplexF64,2}}
 Aw::Array{Array{ComplexF64,2}}
 Hw::Array{Array{ComplexF64,2}}
 H_extern::Array{Array{ComplexF64,2}}
end


mutable struct strKernel
 Kernel::Array{ComplexF64,2}
 Kernel_dagger::Array{ComplexF64,2}
 moment::Array{Float64,2}
 moment_dagger::Array{Float64,2}
 smooth::Array{ComplexF64,2}
 cubic_coeff_transf::Array{ComplexF64,2}
 gamma::Array{ComplexF64,2}
 partial::Array{ComplexF64,2}
end

mutable struct strNumeric
 Egrid::Int64
 ERealAxis::Array{Float64}
 dERealAxis::Array{Float64}
 N_MatsubaraDATA::Int64
 N_Matsubara::Int64
 NumIter::Int64
 mixing_min::Float64
 mixing_max::Float64
 mixingInitial::Float64
 pulay_mixing_history::Int64
 pulay_mixing_step::Int64
 pulay_mixing_start::Int64
end


 struct data_info_
   workDirect      
   inputFile       
   NumFullOrbit    
   num_of_subblock 
   start_cluster   
   N_Matsubara2    
   inverse_temp    
   Asymtotic_HighFreq
 end

 struct real_freq_grid_info_
  EwinOuterRight  
  EwinOuterLeft   
  EwinInnerRight  
  EwinInnerLeft   
  coreWindow      
  coreDense       
  EgridInner      
 end
 struct mem_fit_parm_
  default_model  
  Model_range_right
  Model_range_left
  auxiliary_inverse_temp_range 
  auxTempRenomalFactorInitial  
 end
 struct mixing_parm_
  NumIter              
  mixingInitial        
  mixing_max           
  mixing_min           
  pulay_mixing_start   
  pulay_mixing_step    
  pulay_mixing_history 
 end

struct strInputInfo
 data_info::data_info_
 real_freq_grid_info::real_freq_grid_info_
 mem_fit_parm::mem_fit_parm_
 mixing_parm::mixing_parm_
end


