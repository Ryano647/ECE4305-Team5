%% Modulation Setup 
amps = [-7, -5, -3, -1, 1, 3, 5, 7];
symb = [000, 001, 010, 011, 111, 110, 100, 101];
msgLength = 1024; %Must be divisible by 3
symMap = containers.Map(symb,amps);

% Creates array of 3-bit symbols 
bits = randi([0 1],msgLength,1);

% symbols = [];
% signal = [];
% 
% for k = 1:((length(bits)/3)-1)
%     inc = k*3;
%     symbols = cat(1,symbols,((bits(inc)*100) + (bits(inc+1)*10) + bits(inc+2)));
% end
% 
% for k = 1:length(symbols)
%     signal = cat(1,signal,symMap(symbols(k)));
% end

%% Receiver and Transmitter Setup 
rx=sdrrx('Pluto','OutputDataType','double','SamplesPerFrame',2^15);
tx = sdrtx('Pluto','Gain',-30);

% Transmit sinewave
sine = dsp.SineWave('Frequency',300,...
'SampleRate',rx.BasebandSampleRate,...
'SamplesPerFrame', 2^12,...
'ComplexOutput', true);

mod = comm.QPSKModulator('BitInput',1);
refC = constellation(mod);
modSig = step(mod,bits);

tx.transmitRepeat(modSig); % Transmit continuously

%% Setup Scope
samplesPerStep = rx.SamplesPerFrame/rx.BasebandSampleRate;
steps = 100;
ts = dsp.TimeScope('SampleRate', rx.BasebandSampleRate,...
'TimeSpan', samplesPerStep*steps,...
'BufferLength', rx.SamplesPerFrame*steps);
cd = comm.ConstellationDiagram('ReferenceConstellation',refC);

% Receive and view sine
for k=1:steps
%ts(rx());
cd(rx());
end
