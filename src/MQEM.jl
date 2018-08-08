module MQEM
 include("mem_variable.jl")
 include("MQEM_ftn.jl")
 include("mem_pulay.jl")
 include("mem_write_data.jl")


# using .default_model

 export construct_EnergyGrid
 export data_info_, real_freq_grid_info_, mem_fit_parm_, mixing_parm_, strInputInfo
 export strPhyParm, strRealFreqFtn, strNumeric, strKernel
 export read_matsubara_GreenFtn!, construct_Kernel_inCubicSpline, write_spectral_ftn, search_alpha, get_total_energy, get_total_energy_for_tail, KK_relation, write_results, find_default_parm
 export construct_model_spectrum, gaussianPotential!
end
