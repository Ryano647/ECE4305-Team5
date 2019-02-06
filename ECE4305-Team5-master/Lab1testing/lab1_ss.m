
%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 8;
frameSize = 2^10;
numFrames = 200;
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
cdPLL = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset after PLL');
cdPre.Position(1) = 50;
cdPost.Position(1) = cdPre.Position(1)+cdPre.Position(3)+10;% Place side by side
cdPLL.Position(1) = cdPost.Position(1)+cdPost.Position(3)+10;
%% Impairments
snr = 15;
timingOffset = samplesPerSymbol*0.01; % Samples

%% Generate symbols
data = randi([0 modulationOrder-1], numSamples*2, 1);
mod = comm.DBPSKModulator();
modulatedData = mod.step(data);

%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', samplesPerSymbol);% Set to filterUpsample/2 when introducing timing estimation
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
pllData = [];

pll = comm.SymbolSynchronizer(...
    'SamplesPerSymbol', samplesPerSymbol, ...
    'NormalizedLoopBandwidth', 0.01, ...
    'DampingFactor', 1.0, ...
    'DetectorGain', 10, ...
    'TimingErrorDetector', 'Zero-Crossing (decision-directed)');


for k=1:frameSize:(numSamples)
    
    timeIndex = (k:k+frameSize-1).';
    
    % Filter signal
    filteredTXData = step(TxFlt, modulatedData(timeIndex));
    
    % Pass through channel
    noisyData = step(chan, filteredTXData);
    
    % Time delay signal
    offsetData = step(varDelay, noisyData, k/frameSize*timingOffset); % Variable delay
    
    release(pll);
    % Filter signal
    filteredData = step(RxFlt, offsetData);
    filteredDataRef = real(step(RxFltRef, noisyData));
    pllData = pll(filteredData);
    
    
    % Visualize Error
    step(cdPre,filteredDataRef);
    step(cdPost,filteredData);
    step(cdPLL,pllData);
    pause(0.1); 
end
