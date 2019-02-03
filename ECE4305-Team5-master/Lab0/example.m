%clear, close all
%% Example
% Setup Receiver
rx=sdrrx('Pluto','OutputDataType','double','SamplesPerFrame',2^15, 'CenterFrequency',563.47*10^6);%, 'BasebandSampleRate', 60*10^6);
% Setup Transmitter
tx = sdrtx('Pluto','Gain',-30);
% Transmit sinewave
sine = dsp.SineWave('Frequency',300,...
                    'SampleRate',rx.BasebandSampleRate,...
                    'SamplesPerFrame', 2^12,...
                    'ComplexOutput', true); 
                
%tx.transmitRepeat(sine());
% Transmit continuously
% Setup Scope
samplesPerStep = rx.SamplesPerFrame/rx.BasebandSampleRate;
steps = 3; % ~3 steps to 100ms
%ts = dsp.TimeScope('SampleRate', rx.BasebandSampleRate,...
                   %'TimeSpan', samplesPerStep*steps,...
                   %'BufferLength', rx.SamplesPerFrame*steps);
ts = dsp.SpectrumAnalyzer('SampleRate', 2*10^9,...
                          'FrequencySpan', 'Span and center frequency',...
                          'CenterFrequency', 563.47*10^6);
cd = comm.ConstellationDiagram;
% Receive and view sine
for k=1:steps
    ts(rx());
    cd;
end