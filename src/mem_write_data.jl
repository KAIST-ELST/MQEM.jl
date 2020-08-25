###########################################
# Jae-Hoon Sim, KAIST 2018.01.
###########################################
#############################################################################################
# (start) define Ftns
#############################################################################################
function write_results( NumSubOrbit::Int64, fname_out::Array{String}, fname_contniuedSpectrum::Array{String}, fname_reproduce::Array{String}, kernel,
                         Aw::Array{Array{ComplexF64,2}}, ERealAxis::Array{Float64}, Normalization::Float64, GreenConst::Array{ComplexF64}, Gw_RealPart::Array{Array{ComplexF64,2}},
                        Egrid::Int64, inverse_temp, startOrbit::Int64, cluster::Int64)
    
     #Aw
     for i = 1:NumSubOrbit
         for j=1:NumSubOrbit
             f=open("$(fname_out[i,j])","w")     
             for w=1:Egrid
                 E = ERealAxis[w]
                 write(f,"$E  $(real(Aw[w][i,j])*Normalization )  $(imag(Aw[w][i,j]) *Normalization )\n" )
             end
             close(f)
         end
     end
    
    
    #G_retarded(w)
    MatrixValuedOutPut = Array{Array{ComplexF64,2}}(undef,Egrid)
    for w=1:Egrid
          MatrixValuedOutPut[w] = Gw_RealPart[w] -pi*(Aw[w])im 
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

    #reporduced G(iwn)
    GreenFtn_rep = kernel.Kernel *  Aw ;
    GreenFtn_rep *= Normalization ;
    for i = 1:NumSubOrbit
        for j=1:NumSubOrbit
        f=open( fname_reproduce[i,j],"w")
        for iw=1:length(GreenFtn_rep)
            z = ((2.0*iw-1)*pi/inverse_temp)im;
            ftn = GreenFtn_rep[iw][i,j] + GreenConst[i,j] 
            write(f,"$(imag(z))  $(real(ftn))    $(imag(ftn))  \n" )
        end
        end
        close(f)
    end
    
    f=open( "Sw_SOLVER.full_fromRetardedSw.dat_$(cluster)_$(cluster)" ,"w")
    for iw=1:length(GreenFtn_rep)
        for i = 1:NumSubOrbit
            for j=1:NumSubOrbit
                ftn = GreenFtn_rep[iw][i,j] + GreenConst[i,j] 
                if abs(ftn)>1e-5
                    write(f,"$(iw-1) $( startOrbit+i-1)  $(startOrbit+j-1) $(real(ftn))    $(imag(ftn))  \n" )
                end
            end
        end
    end
    close(f)
end





function write_spectral_ftn(NumSubOrbit::Int64,Normalization::Float64,  numeric, SpecFtn::Array{Array{ComplexF64,2}}, kernel::strKernel,fname_out::Array{String} ,Deco::String  )


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

