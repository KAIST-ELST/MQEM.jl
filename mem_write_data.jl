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
#               write(f,"$E  $(real(Aw[w][i,j])/dERealAxis[w] *Normalization )  $(imag(Aw[w][i,j])/dERealAxis[w] *Normalization )\n" )
               write(f,"$E  $(real(Aw[w][i,j])*Normalization )  $(imag(Aw[w][i,j]) *Normalization )\n" )
           end
           close(f)
       end
   end


  MatrixValuedOutPut= Array{Array{Complex128,2}}(Egrid)
  for w=1:Egrid
        MatrixValuedOutPut[w] = Aw_RealPart[w] -pi*(Aw[w])im 
        if inversion_method
            MatrixValuedOutPut[w] *= trace(eye(MatrixValuedOutPut[w]))
            MatrixValuedOutPut[w] =  (    z*eye(MatrixValuedOutPut[w])-inv(MatrixValuedOutPut[w])    )
        end
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





function write_spectral_ftn(NumSubOrbit::Int64,Normalization::Float64,  numeric, SpecFtn::Array{Array{Complex128,2}},fname_out::Array{String} ,Deco::String  )
# FtnWrite = smooth * SpecFtn *Normalization
 FtnWrite =  SpecFtn *Normalization
 for i = 1:NumSubOrbit
     for j=1:NumSubOrbit
         f=open("$(fname_out[i,j])$(Deco)","w")     
         for w=1:numeric.Egrid
             E = numeric.ERealAxis[w]
#             write(f,"$E  $(real(FtnWrite[w][i,j])/numeric.dERealAxis[w]  )  $(imag(FtnWrite[w][i,j])/numeric.dERealAxis[w] )\n" )
             write(f,"$E  $(real(FtnWrite[w][i,j]) )  $(imag(FtnWrite[w][i,j]) )\n" )
         end
         close(f)
     end
 end

end
