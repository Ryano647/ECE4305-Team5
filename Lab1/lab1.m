
%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 8;
frameSize = 2^10;
numFrames = 200;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterSymbolSpan = 4;
rollOff = 0.2; %default 0.2

%% Visuals
cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');
cdPost = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset');
cdPre.Position(1) = 50;
cdPost.Position(1) = cdPre.Position(1)+cdPre.Position(3)+10;% Place side by side

%% Impairments
snr = 15;
timingOffset = samplesPerSymbol*0.01; % Samples

%% Generate symbols
data = randi([0 modulationOrder-1], numSamples*2, 1);
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
    'DecimationFactor', samplesPerSymbol); % Set to filterUpsample/2 when introducing timing estimation
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
allFilteredData = [];
timeIndex;
for k=1:frameSize:(numSamples - frameSize)
    
    timeIndex = (k:k+frameSize-1).';
    % Filter signal
    filteredTXData = step(TxFlt, modulatedData(timeIndex));
    
    % Pass through channel
    noisyData = step(chan, filteredTXData);
    
    % Time delay signal
    offsetData = step(varDelay, noisyData, k/frameSize*timingOffset); % Variable delay
    
    % Filter signal
    filteredData = step(RxFlt, offsetData);
    if k == 1 
        allFilteredData= filteredData;
    end
    filteredDataRef = step(RxFltRef, noisyData);
    allFilteredData = cat(2,allFilteredData,filteredData);
    
    % Visualize Error
    step(cdPre,filteredDataRef);
    step(cdPost,filteredData);pause(0.1); 
    
end

plot(timeIndex, real(allFilteredData(:,1)));



