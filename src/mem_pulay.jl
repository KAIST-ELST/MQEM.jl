####################################################
#(start) def:Pulay Mixing
####################################################


function matrixFtn_innerProduct( Aw::Array{Array{Complex128,2}} , Bw::Array{Array{Complex128,2}}, weight::Array{Float64})
            lengthOfVector = length(Aw)
            normVector_temp = Array{Complex128}(lengthOfVector)
            for w=1:lengthOfVector
                  normVector_temp[w] = trace( ( Aw[w]') * Bw[w] )  * weight[w]
            end
            return sum(normVector_temp)
end

#Data store
type pulayInfo
    mixingMemory::Int64
    mixingStart::Int64
    mixingStep::Int64
    inputHistory::Array{Array{Array{Complex128,2}}}
    residHistory::Array{Array{Array{Complex128,2}}}
    additionHistory::Array{Array{Array{Complex128,2}}}
    PulayMatrix::Array{Float64}
    previousMixing::String
end

#Initializer
function pulayInitializer!(MixInfo::pulayInfo)
    MixInfo.inputHistory = Array{Array{Array{Complex128,2}}}(MixInfo.mixingMemory);
    MixInfo.residHistory = Array{Array{Array{Complex128,2}}}(MixInfo.mixingMemory);
#    MixInfo.additionHistory = Array{Array{Array{Complex128,2}}}(MixInfo.mixingMemory);
    MixInfo.PulayMatrix =  zeros(Float64, MixInfo.mixingMemory+1, MixInfo.mixingMemory+1);
    MixInfo.PulayMatrix[1, 2:end]=-1;
    MixInfo.PulayMatrix[2:end, 1]=-1;
end



#Pulay Mixing
#module pulay_DIIS_mixing
function pulayMixing!(iter::Int64, mixing::Float64,
         FtnIn::Array{Array{Complex128,2}},  FtnResid::Array{Array{Complex128,2}}, weight,
         MixInfo::pulayInfo)
    
    if iter<MixInfo.mixingStart #Simple Mixing
        FtnMixed = FtnIn      + mixing*FtnResid
    elseif  iter>=MixInfo.mixingStart #Pulay matrix update
        removeComp = MixInfo.mixingMemory  
        maxValue = 0
           
        #construct Aw, shift elements
        MixInfo.inputHistory[2:removeComp] = MixInfo.inputHistory[1:removeComp-1]
        MixInfo.inputHistory[1] = FtnIn
 
        MixInfo.residHistory[2:removeComp] = MixInfo.residHistory[1:removeComp-1]
        MixInfo.residHistory[1] = FtnResid

 
        MixInfo.PulayMatrix[:, 3:removeComp+1] = MixInfo.PulayMatrix[:, 2:removeComp-1+1]
        MixInfo.PulayMatrix[3:removeComp+1, :] = MixInfo.PulayMatrix[2:removeComp-1+1, :]
        for i=1:min(iter-MixInfo.mixingStart+1, MixInfo.mixingMemory)
            MixInfo.PulayMatrix[2,i+1] = real( matrixFtn_innerProduct( MixInfo.residHistory[1] , MixInfo.residHistory[i],   weight ))
            MixInfo.PulayMatrix[i+1,2] = MixInfo.PulayMatrix[2,i+1]
        end


        dim = min(MixInfo.mixingMemory+1, iter-MixInfo.mixingStart+1);
        temp = MixInfo.PulayMatrix[1:dim, 1:dim];

        #To obtain best guessing  Aw
        if (iter%abs(MixInfo.mixingStep)==0) &&  MixInfo.mixingMemory>1  &&
           iter >= MixInfo.mixingStart + MixInfo.mixingMemory  &&
            det(Symmetric(temp)) != 0 

            #Pulay interpol
            #Linear Solver
#            rhs=zeros(MixInfo.mixingMemory+1)
            rhs=zeros(dim)
            rhs[1]=-1
            alpha = \(Symmetric(temp),rhs )
            alpha = alpha[2:end]
            FtnOpt = sum(alpha .* MixInfo.inputHistory[1:dim-1])
            FtnOptResd = sum(alpha.*MixInfo.residHistory[1:dim-1])
            if any(isnan, alpha) || any(isinf, alpha)
               FtnOpt = FtnIn
               FtnOptResd = FtnResid
            end
            if(MixInfo.mixingStep <0) 
                  FtnMixed = (1-mixing) * MixInfo.inputHistory[1] + mixing*(  FtnOpt + mixing * FtnOptResd  )
            else  FtnMixed =                                                  FtnOpt + mixing * FtnOptResd 
            end
            MixInfo.previousMixing = "pulay"
        else
            #Simple mixing
            FtnMixed = FtnIn + mixing*FtnResid
            MixInfo.previousMixing = "simple"
        end
    end
    return  FtnMixed
end
####################################################
# (end) def:Pulay Mixing
####################################################

