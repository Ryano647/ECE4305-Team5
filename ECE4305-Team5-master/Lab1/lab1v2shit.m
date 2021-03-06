close all

%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 8;
frameSize = 8;
numFrames = 9;
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
%allData = cat(2,allData,data);
mod = comm.BPSKModulator();
demod = comm.BPSKDemodulator();
modulatedData = mod.step(data);
%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'RolloffFactor', rollOff);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', 1,...
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
downsampledRxData= [];
downsampledRxDatawoff= [];
demodulatedDatawoff= [];
%timeIndex = 0;
w=1; %Downsample increment
for k=1:frameSize:(numSamples)
    w=w+1;

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
    
    % Matched Filter signal
    filteredData = step(RxFlt, offsetData);

    filteredDataRef = step(RxFltRef, noisyData);
    
    %PLL
    %by Sean Brady
   % Error_Signal = real(
    
    %downsample
    %by Sean Brady
    for downsampler_iterator = 1:frameSize
    downsampledRxData= [downsampledRxData,filteredDataRef(...
        (downsampler_iterator-1)*samplesPerSymbol+1)];
    downsampledRxDatawoff= [downsampledRxDatawoff,filteredData(...
        (downsampler_iterator-1)*samplesPerSymbol+1)];
    end
    %demod
    %by Sean Brady
    
    
    demodulatedData = demod.step(downsampledRxData');
    demodulatedDatawoff = demod.step(downsampledRxDatawoff');
    
    if k == 1 
        allFilteredDatawoff= filteredData;
        allFilteredData = filteredDataRef;
        allFilteredTXData = filteredTXData;
        allnoisyData = noisyData;
        alldownsampledRxData = downsampledRxData;
        alldemodulatedData = demodulatedData;
        alldownsampledRxDatawoff = downsampledRxDatawoff;
    end
    if k >1
    allFilteredDatawoff = cat(2,allFilteredDatawoff,filteredData);
    allFilteredData = cat(2,allFilteredData,filteredDataRef);
    allFilteredTXData = cat(2,allFilteredTXData,filteredTXData);
    allnoisyData = cat(2,allnoisyData,noisyData);
    alldownsampledRxData = cat(2,alldownsampledRxData,downsampledRxData);
    alldemodulatedData = cat(2,alldemodulatedData,demodulatedData);
    alldownsampledRxDatawoff = cat(2,alldemodulatedDatawoff,demodulatedDatawoff);
    end
    % Visualize Error
%     step(cdPre,filteredDataRef);
%     step(cdPost,filteredData);pause(0.1); 
    
end


%% Transmit and Receive Plots
%By R O'brian

timeIndex = timeIndex / (1e6);


%Code for creating the correct amountdata points to plot Tx and Rx plots
%Sean Brady

%Constants
% Spaced_out_recieved_data=[];
% Spaced_out_Tx_data = [];
% longTx = [];
% longNoise = [];
% %creating a time variable to match the size we need to plot Tx and Rx
% %filtered data
% tas = 1:(numSamples-1)/((frameSize*samplesPerSymbol*numFrames)-1):numSamples;
% 
% for Frame_iterator= 1:numFrames
%     Symbols_in_frame_iterator = 1;
%     %Getting all of filteredTXData and noisydata in 1 array
%     for All_sample_iterator= 1:(samplesPerSymbol*frameSize)
%         longTx=[longTx, allFilteredTXData(All_sample_iterator,Frame_iterator)];
%         longNoise = [longNoise, allnoisyData(All_sample_iterator,Frame_iterator)];
%     end 
%     %creating a array with properly spaced filteredData points to plot with
%     %filterTX data
%     for Symbols_in_frame_iterator= 1:frameSize
%         hjk = 1;
%         if Symbols_in_frame_iterator ==1
% %             Spaced_out_recieved_data(((Frame_iterator-1)*frameSize*samplesPerSymbol)+1)...
% %                 = allFilteredData(1,Frame_iterator);
%             Spaced_out_Tx_data(((Frame_iterator-1)*frameSize*samplesPerSymbol)+1)...
%                 = modulatedData(1+(samplesPerSymbol*(Frame_iterator-1)));
% 
%             for hjk =1:(samplesPerSymbol-1)
% %                 Spaced_out_recieved_data = [Spaced_out_recieved_data,0];
%                 Spaced_out_Tx_data = [Spaced_out_Tx_data,0];
%             end
%         else
% %             Spaced_out_recieved_data((samplesPerSymbol*(Symbols_in_frame_iterator-1)+...
% %                 ((Frame_iterator-1)*frameSize*samplesPerSymbol)+1))...
% %                 = allFilteredData(Symbols_in_frame_iterator,Frame_iterator);
%             Spaced_out_Tx_data((samplesPerSymbol*(Symbols_in_frame_iterator-1)+...
%                 ((Frame_iterator-1)*frameSize*samplesPerSymbol)+1))...
%                 = modulatedData(Symbols_in_frame_iterator+(samplesPerSymbol*(Frame_iterator-1)));
%             for hjk =1:(samplesPerSymbol-1)
% %                 Spaced_out_recieved_data = [Spaced_out_recieved_data,0];
%                 Spaced_out_Tx_data = [Spaced_out_Tx_data,0];
%             end
%         end 
%     end
% end


% for psd= 1:numFrames
%     stem(tas((1+(psd-1)*samplesPerSymbol*frameSize):(psd*samplesPerSymbol*frameSize)), real(Spaced_out_recieved_data((1+(psd-1)*samplesPerSymbol*frameSize):(psd*samplesPerSymbol*frameSize))));
%     hold on;
%     stem(tas((1+(psd-1)*samplesPerSymbol*frameSize):(psd*samplesPerSymbol*frameSize)), real(Spaced_out_Tx_data((1+(psd-1)*samplesPerSymbol*frameSize):(psd*samplesPerSymbol*frameSize)))); %1:length(noisyData) 
%     hold on;
%     plot(tas((1+(psd-1)*samplesPerSymbol*frameSize):(psd*samplesPerSymbol*frameSize)), real(longTx((1+(psd-1)*samplesPerSymbol*frameSize):(psd*samplesPerSymbol*frameSize)))); %1:length(filteredTXData)
%     hold off;
%     title('Transmit and Receive Plot')
%     xlabel('Time (ms)')
%     ylabel('Amplitude')
%     idc = input('Press return to continue');
% end




