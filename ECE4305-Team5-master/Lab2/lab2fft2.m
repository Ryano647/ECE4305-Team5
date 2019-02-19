%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 1;
frameSize = 2^10;
numFrames = 100;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterUpsample = 4; %original value 4
filterSymbolSpan = 8;
fftOrder = 2^10;
fftIndex = 1;
NormalizedLoopBandwidth = 0.000001;
DampingFactor = [0.9];

%% Configure LF and PI
LoopFilter = dsp.IIRFilter('Structure', 'Direct form II transposed', ...
'Numerator', [1 0], 'Denominator', [1 -1]);
Integrator = dsp.IIRFilter('Structure', 'Direct form II transposed', ...
'Numerator', [0 1], 'Denominator', [1 -1]);

%% Calculate range estimates
    NormalizedPullInRange = min(1, 2*pi*sqrt(2)*DampingFactor*...
    NormalizedLoopBandwidth);
    MaxFrequencyLockDelay = (4*NormalizedPullInRange^2)/...
        (NormalizedLoopBandwidth)^3;
    MaxPhaseLockDelay = 1.3/(NormalizedLoopBandwidth);

%% Impairments
snr = 15;
frequencyOffsetHz = 1e5; % Offset in hertz
phaseOffset = 0; % Radians
frequencyOffsetHz = (frequencyOffsetHz) * 1;

%% Generate symbols
data = randi([0 samplesPerSymbol], numSamples, 1);
mod = comm.QPSKModulator();
modulatedData = mod.step(data);

%% Add TX Filter
TxFlt = comm.RaisedCosineTransmitFilter('OutputSamplesPerSymbol', filterUpsample, 'FilterSpanInSymbols', filterSymbolSpan);
filteredData = step(TxFlt, modulatedData);

%% Add noise
noisyData = awgn(filteredData,snr);%,'measured');

%% Setup visualization object(s)
sa = dsp.SpectrumAnalyzer('SampleRate',sampleRateHz,'ShowLegend',true);

%% Calculate coefficients for FFC
    PhaseRecoveryLoopBandwidth = NormalizedLoopBandwidth*samplesPerSymbol;
    PhaseRecoveryGain = samplesPerSymbol;
    PhaseErrorDetectorGain = log2(4); DigitalSynthesizerGain = -1;
    theta = PhaseRecoveryLoopBandwidth/...
    ((DampingFactor + 0.25/DampingFactor)*samplesPerSymbol);
    delta=1+ 2*DampingFactor*theta + theta*theta;
    % G1
    ProportionalGain = (4*DampingFactor*theta/delta)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);
    % G3
    IntegratorGain = (4/samplesPerSymbol*theta*theta/delta)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);

%% Model of error
% Add frequency offset to baseband signal

% Precalculate constant(s)
normalizedOffset = 1i.*2*pi*frequencyOffsetHz./sampleRateHz;
sigNoMod = [];
offsetEstimates = zeros(floor(length(noisyData)/fftOrder),1);
indexToHz = sampleRateHz/(modulationOrder*fftOrder);
corrected = [];

 
Phase= 0; previousSample = complex(0);
LoopFilter.release();Integrator.release();

offsetData = zeros(size(noisyData));
outputr = zeros(size(noisyData)); 

for k=1:frameSize:numSamples*filterUpsample
    
    LoopFilter.release();Integrator.release();
    % Create phase accurate vector
    timeIndex = (k:k+frameSize-1).';
    fftIndex = (k:k+fftOrder-1);
    freqShift = exp(normalizedOffset*timeIndex + phaseOffset);
%     freqShift = cos(normalizedOffset + phaseOffset);
%     freqShift = exp(normalizedOffset*timeIndex.'); % textbook example
    
    % Offset data and maintain phase between frames
    offsetData(timeIndex) = (noisyData(timeIndex).*freqShift);
%     offsetData = filteredData.*freqShift; % textbook example
     sigNoMod = offsetData(timeIndex).^modulationOrder;
    
     freqHist = abs(fft(sigNoMod));
    % Determine most likely offset
    [~,maxInd] = max(freqHist);
    offsetInd = maxInd - 1;
    if maxInd>=fftOrder/2 % Compensate for spectrum shift
    offsetInd = offsetInd - fftOrder;
    end
    % Convert to Hz from normalized frequency index
    offsetEstimates(k) = offsetInd * indexToHz;
    corrected = noisyData - offsetData;
    for t = 1:frameSize  
    output(timeIndex(t)) = corrected(timeIndex(t))*exp(1i*Phase);
    % PED
    phErr(timeIndex(t)) = sign(real(previousSample)).*imag(previousSample)...
        - sign(imag(previousSample)).*real(previousSample);
    % Loop Filter
    loopFiltOut = step(LoopFilter,phErr(timeIndex(t))*IntegratorGain);
    % Direct Digital Synthesizer
    DDSOut = step(Integrator,phErr(timeIndex(t))*ProportionalGain + loopFiltOut);
    Phase = DigitalSynthesizerGain * DDSOut;
    previousSample = output(timeIndex(t));
    end

    % Visualize Error
    step(sa,[noisyData(timeIndex),output(timeIndex)]);pause(0.01); %#ok<*UNRCH>
    
    %,offsetData(timeIndex),corrected(timeIndex), 
end
e= zeros(size(output));
esum = 0;
dsum = 0;
n=0;
for k= 1:length(output)
    e(k)= (real(output(k)) -real(noisyData(k)))^2 +(imag(output(k)) -imag(noisyData(k)))^2;
    d(k) = (real(output(k))^2 +imag(output(k))^2);
    esum= esum + e(k);
    dsum = dsum + d(k);
    n = n+1;
end
eVM = 100*((esum/n)/(dsum/n))^1/2;
%% Plot
df = sampleRateHz/frameSize;
frequencies = -sampleRateHz/2:df:sampleRateHz/2-df;
spec = @(sig) fftshift(10*log10(abs(fft(sig))));
h = plot(frequencies, spec(noisyData(timeIndex)),...
     frequencies, spec(offsetData(timeIndex)),...
     frequencies, spec(corrected(timeIndex)));
 
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
legend('Original','Offset','Corrected','Location','Best');
NumTicks = 5;L = h(1).Parent.XLim;
set(h(1).Parent,'XTick',linspace(L(1),L(2),NumTicks))