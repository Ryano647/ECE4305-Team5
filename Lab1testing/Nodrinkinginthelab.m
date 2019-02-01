tas = 1:(numSamples-1)/((frameSize*samplesPerSymbol)-1):numSamples;
jksd=[];
for ert= 1:numFrames
    ids = 1;
    for ids= 1:frameSize
        hjk = 1;
        if ids ==1
            jksd(((ert-1)*frameSize*samplesPerSymbol)+1)= real(allFilteredData(1,ert+1));
            for hjk =1:(samplesPerSymbol-1)
                jksd = [jksd,0];
            end
        else
            jksd((samplesPerSymbol*(ids-1)+((ert-1)*frameSize*samplesPerSymbol)+1))=real(allFilteredData(ids,ert+1));
            hds(ids)=((samplesPerSymbol*(ids-1))+((ert-1)*frameSize*samplesPerSymbol)+1);
            for hjk =1:(samplesPerSymbol-1)
                %
                jksd = [jksd,0];
            end
        end 
    end
end