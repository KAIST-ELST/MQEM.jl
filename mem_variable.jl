struct strPhyParm
 inverse_temp::Float64
 NumFullOrbit::Int64
 NumSubOrbit::Int64
 num_of_subblock::Int64
end


struct strImagFreqFtn
    GreenFtn
    GreenConstFull
    GreenConst::Array{Complex128,2}
    Normalization::Float64
    bestUnitary
    moments1
    moments2
    moments3
    moments4
    moments5
end


mutable struct strRealFreqFtn
 Spectral_default_model::Array{Array{Complex128,2}}
 logSpectral_default_model::Array{Array{Complex128,2}}
 Aw::Array{Array{Complex128,2}}
 Hw::Array{Array{Complex128,2}}
 H_extern::Array{Array{Complex128,2}}
end


mutable struct strKernel
 auxiliary_inverse_temp_range::Array{Float64}
 Kernel::Array{Complex128,2}
 Kernel_dagger::Array{Complex128,2}
 moment::Array{Float64,2}
 smooth::Array{Complex128,2}
 cubic_coeff_transf::Array{Complex128,2}
 gamma::Array{Complex128,2}
# svdInfo
# svdmax::Int64
# svdKernel
# svdinteraction_V
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
