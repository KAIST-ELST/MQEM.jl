###########################################
# Jae-Hoon Sim, KAIST 2018.01.
###########################################
#############################################################################################
# (start) define Ftns
#############################################################################################
function write_results( NumSubOrbit::Int64, fname_out::Array{String}, fname_contniuedSpectrum::Array{String},Aw::Array{Array{Complex128,2}}, ERealAxis::Array{Float64}, Normalization::Float64, GreenConst::Array{Complex128}, Aw_RealPart::Array{Array{Complex128,2}}, Egrid::Int64)

   for i = 1:NumSubOrbit
       for j=1:NumSubOrbit
           f=open("$(fname_out[i,j])_optimal","w")     
           for w=1:Egrid
               E = ERealAxis[w]
               write(f,"$E  $(real(Aw[w][i,j])*Normalization )  $(imag(Aw[w][i,j]) *Normalization )\n" )
           end
           close(f)
       end
   end


  MatrixValuedOutPut= Array{Array{Complex128,2}}(Egrid)
  for w=1:Egrid
        MatrixValuedOutPut[w] = Aw_RealPart[w] -pi*(Aw[w])im 
        MatrixValuedOutPut[w] *=  Normalization;
        MatrixValuedOutPut[w] += GreenConst;
  end
  for i = 1:NumSubOrbit
      for j=1:NumSubOrbit
          f=open(fname_contniuedSpectrum[i,j],"w")
          for w=1:Egrid
              E = ERealAxis[w]
              write(f, "$E $(real(MatrixValuedOutPut[w][i,j])) $(imag(MatrixValuedOutPut[w][i,j]))\n")
          end
          close(f)
      end
  end
end





function write_spectral_ftn(NumSubOrbit::Int64,Normalization::Float64,  numeric, SpecFtn::Array{Array{Complex128,2}}, kernel::strKernel,fname_out::Array{String} ,Deco::String  )

######
#####   num_seg_coeff = 2*(numeric.Egrid-1)
#####
#####   T = zeros(Float64, num_seg_coeff,  numeric.Egrid)
#####   B = zeros(Float64, num_seg_coeff, num_seg_coeff)
#####   Kernel_from_cubic= zeros(Complex128, numeric.Egrid, num_seg_coeff)
#####
#####
#####   for j=1:numeric.Egrid-1
#####     s = 2*(j-1)
#####     T[s+1, j]   = 1.0
#####     T[s+2, j+1] = 1.0
#####   end
#####
#####   for j=1:numeric.Egrid-1
#####       s = 2*(j-1)
#####       B[s+1, s+2] = 1
#####
#####       B[s+2, s+1] = (numeric.ERealAxis[j+1]-numeric.ERealAxis[j])^1
#####       B[s+2, s+2] = 1
#####   end
#####
######   f1(wj, w) = 1.0/2.0 *(w-wj)^2
#####   f1(wj, w) = 0.0
#####   f2(wj, w) = w
######   f3(iwn, wj, w ) =  -w   -(iwn-wj)   * log(-iwn + wj +w)
######   f4(iwn, wj, w ) =       -             log(-iwn + wj +w)
#####   for n=1:numeric.Egrid
#####    for j=1:numeric.Egrid-1
#####       s = 2*(j-1)
######       wj = numeric.ERealAxis[j]  - numeric.ERealAxis[n]
######       Kernel_from_cubic[n,2*(j-1)+1] = (f3(-igamma, wj, dSegment[j]) - f3(-igamma, wj, 0))  -  (f3(igamma, wj, dSegment[j]) - f3(igamma, wj, 0))
######       Kernel_from_cubic[n,2*(j-1)+2] = (f4(-igamma, wj, dSegment[j]) - f4(-igamma, wj, 0))  -  (f4(igamma, wj, dSegment[j]) - f4(igamma, wj, 0))
#####
#####       En=numeric.ERealAxis[n]
#####       Ej1=numeric.ERealAxis[j+1]
#####       Ej=numeric.ERealAxis[j]
#####
#####       upper_cut = minimum([En+blur_width, Ej1])
#####       lower_cut = maximum([En-blur_width, Ej])
#####       if upper_cut > lower_cut
#####         Kernel_from_cubic[n,s+1] =  f1(Ej, upper_cut) - f1(Ej,lower_cut)
#####         Kernel_from_cubic[n,s+2] =  f2(Ej, upper_cut) - f2(Ej,lower_cut)
#####       end
#####    end
#####   end
#####
#####   println( inv(B)*T*SpecFtn )
#####   println 
#####   println( inv(B)*T*SpecFtn )
#####   ###################################################

 FtnWrite = kernel.smooth * SpecFtn *Normalization
# FtnWrite =  SpecFtn *Normalization
 for i = 1:NumSubOrbit
     for j=1:NumSubOrbit
         f=open("$(fname_out[i,j])$(Deco)","w")     
         for w=1:numeric.Egrid
             E = numeric.ERealAxis[w]
             write(f,"$E  $(real(FtnWrite[w][i,j]) )  $(imag(FtnWrite[w][i,j]) )\n" )
         end
         close(f)
     end
 end

end
