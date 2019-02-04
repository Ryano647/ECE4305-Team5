zFittedTxData = zeros(1,frameSize*samplesPerSymbol);
zFittedRxData = zeros(1,frameSize*samplesPerSymbol);
for w= 1:frameSize
    zFittedTxData(samplesPerSymbol*(w-1)+1) = modulatedData(w);
    zFittedRxData(samplesPerSymbol*(w-1)+1) = allDownsampledRxData(w);
end

