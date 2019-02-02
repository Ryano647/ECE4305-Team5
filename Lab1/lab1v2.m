close all

%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 8;
frameSize = 8;
numFrames = 4;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterSymbolSpan = 4;
rollOff = 0.2; %default 0.2

%% Visuals
% cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
%     'Name','Baseband');
% cdPost = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
%     'SymbolsToDisplaySource','Property',...
%     'SymbolsToDisplay',frameSize/2,...
%     'Name','Baseband with Timing Offset');
% cdPre.Position(1) = 50;
% cdPost.Position(1) = cdPre.Position(1)+cdPre.Position(3)+10;% Place side by side

%% Impairments
snr = 200;
timingOffset = samplesPerSymbol*0.01; % Samples

%% Generate symbols
data = randi([0 modulationOrder-1], numSamples, 1);
%data = [1,0,0,0,0,0,0,0]';
% allData = cat(2,allData,data);
mod = comm.DBPSKModulator();
modulatedData = mod.step(data);

%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'RolloffFactor', rollOff);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', samplesPerSymbol,...
    'RolloffFactor', rollOff); % Set to filterUpsample/2 when introducing timing estimation
RxFltRef = clone(RxFlt);

%% Add noise source
chan = comm.AWGNChannel( ...
    'NoiseMethod',  'Signal to noise ratio (SNR)', ...
    'SNR',          snr, ...
    'SignalPower',  1, ...
    'RandomStream', 'mt19937ar with seed');

%% Add delay
varDelay = dsp.VariableFractionalDelay;

%% Setup visualization object(s)
sa = dsp.SpectrumAnalyzer('SampleRate',sampleRateHz,'ShowLegend',true);

%% Model of error
% Add timing offset to baseband signal
filteredData = [];
allFilteredDatawoff = [];
allFilteredData = [];
%timeIndex = 0;
for k=1:frameSize:(numSamples)
    

%     if k == (25*frameSize+1)
%         djssd=1;
%     elseif k == (50*frameSize+1)
%         djssd=2;
%     elseif k == (75*frameSize+1)
%         djssd=3;
%     end
    
    timeIndex = (k:k+frameSize-1).';
    
    % Filter signal
    filteredTXData = step(TxFlt, modulatedData(timeIndex));
    
    % Pass through channel
    noisyData = step(chan, filteredTXData);
    
    % Time delay signal
    offsetData = step(varDelay, noisyData, k/frameSize*timingOffset); % Variable delay
    
    % Filter signal
    filteredData = step(RxFlt, offsetData);

    filteredDataRef = step(RxFltRef, noisyData);
    if k == 1 
        allFilteredDatawoff= filteredData;
        allFilteredData = filteredDataRef;
        allFilteredTXData = filteredTXData;
        allnoisyData = noisyData;
    end
    if k >1
    allFilteredDatawoff = cat(2,allFilteredDatawoff,filteredData);
    allFilteredData = cat(2,allFilteredData,filteredDataRef);
    allFilteredTXData = cat(2,allFilteredTXData,filteredTXData);
    allnoisyData = cat(2,allnoisyData,noisyData);
    end
    % Visualize Error
%     step(cdPre,filteredDataRef);
%     step(cdPost,filteredData);pause(0.1); 
    
end

allFilteredTXData = real(allFilteredTXData);
allFilteredData = real(allFilteredData);
allnoisyData = real(allnoisyData);

%% Transmit and Receive Plots
%By R O'brian

timeIndex = timeIndex / (1e6);


%Code for creating the correct amountdata points to plot Tx and Rx plots
%Sean Brady

%Constants
jksd=[];
longTx = [];
longNoise = [];
%creating a time variable to match the size we need to plot Tx and Rx
%filtered data
tas = 1:(numSamples-1)/((frameSize*samplesPerSymbol*numFrames)-1):numSamples;

for ert= 1:numFrames
    ids = 1;
    %Getting all of filteredTXData and noisydata in 1 array
    for yui= 1:(samplesPerSymbol*frameSize)
        longTx=[longTx, allFilteredTXData(yui,ert)];
        longNoise = [longNoise, allnoisyData(yui,ert)];
    end 
    %creating a array with properly spaced filteredData points to plot with
    %filterTX data
    for ids= 1:frameSize
        hjk = 1;
        if ids ==1
            
            jksd(((ert-1)*frameSize*samplesPerSymbol)+1)...
                = real(allFilteredData(1,ert));
            
            for hjk =1:(samplesPerSymbol-1)
                jksd = [jksd,0];
            end
        else
            jksd((samplesPerSymbol*(ids-1)+((ert-1)*frameSize*samplesPerSymbol)+1))=real(allFilteredData(ids,ert));
            hds(ids)=((samplesPerSymbol*(ids-1))+((ert-1)*frameSize*samplesPerSymbol)+1);
            for hjk =1:(samplesPerSymbol-1)
                jksd = [jksd,0];
            end
        end 
    end
end



 stem(tas, jksd);
 hold on;
 plot(tas, longTx,'-o','MarkerIndices',1:length(filteredTXData));
 hold on;
 plot(tas, longNoise,'-o','MarkerIndices',1:length(noisyData));
 hold off;
 title('Transmit and Receive Plot')
 xlabel('Time (ms)')
 ylabel('Amplitude')




