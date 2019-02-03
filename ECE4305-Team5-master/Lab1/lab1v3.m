
%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 9;
frameSize = 8;
numFrames = 3;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterSymbolSpan = 4;

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
snr = 2000;
timingOffset = samplesPerSymbol*0.01; % Samples

%% Generate symbols
data = randi([0 modulationOrder-1], numSamples, 1);
mod = comm.BPSKModulator();
demod = comm.BPSKDemodulator();
modulatedData = mod.step(data);

%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', 1);% Set to filterUpsample/2 when introducing timing estimation
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
downsampledRxData = [];
allDownsampledRxData = [];
allDemodulatedData = [];
Phase_Shift = 0;
for k=1:frameSize:(numSamples)
    downsampledRxData = [];
    timeIndex = (k:k+frameSize-1).';
    
    % Filter signal Note* Size=frameSize*SamplesPerSymbol
    filteredTXData = step(TxFlt, modulatedData(timeIndex));
    
    % Pass through channel
    noisyData = step(chan, filteredTXData);
    
    % Time delay signal
    offsetData = step(varDelay, noisyData, k/frameSize*timingOffset); % Variable delay
    
    % Filter signal
    filteredData = step(RxFlt, offsetData);
    filteredDataRef = step(RxFltRef, noisyData);
    
    %downsample filtered data
    %by Sean Brady
    for i =1:frameSize
        downsampledRxData = [downsampledRxData,filteredDataRef((i-1)*samplesPerSymbol+1)];
        if i ==8
            allDownsampledRxData = [allDownsampledRxData,downsampledRxData];
        end
    end
    
    %PLL
    %by Sean Brady
    
    %asl = symbol index
    
    asl =5;
    error_sig = real(filteredData(1+(samplesPerSymbol*(asl-1))...
        +(samplesPerSymbol-1)/2)*1/sampleRateHz+Phase_Shift)*(sign(real(...
        filteredData(1+samplesPerSymbol*(asl'-1))*1/sampleRateHz+Phase_Shift))-sign(real(...
        filteredData(1+samplesPerSymbol*(asl))*1/sampleRateHz+Phase_Shift)));
    
    %Demod
    %by Sean Brady
    
    demodulatedData = demod.step(downsampledRxData');
    allDemodulatedData = [allDemodulatedData,demodulatedData];
    
    
    % Visualize Error
%     step(cdPre,filteredDataRef);
%     step(cdPost,filteredData);pause(0.1); 
    
end
timeIndex = timeIndex / (1e6);

%Producing Fig 3, Just lastest frame
tas = 1:(frameSize-1)/((frameSize*samplesPerSymbol)-1):frameSize;
zFittedTxData = zeros(1,frameSize*samplesPerSymbol);
zFittedRxData = zeros(1,frameSize*samplesPerSymbol);
for w= 1:frameSize
    zFittedTxData(samplesPerSymbol*(w-1)+1) = modulatedData(w+frameSize*(numFrames-1));
    zFittedRxData(samplesPerSymbol*(w-1)+1) = allDownsampledRxData(w+frameSize*(numFrames-1));
end

stem(tas,real(zFittedTxData));
hold on;
stem(tas,real(zFittedRxData)); %1:length(noisyData) 
hold on;
plot(tas,real(filteredDataRef)); %1:length(filteredTXData)
hold on;
plot(tas,real(filteredTXData));
hold off;
title('Transmit and Receive Plot')
xlabel('Time (ms)')
ylabel('Amplitude')

