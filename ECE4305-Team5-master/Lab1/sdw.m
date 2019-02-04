    rctFilt = comm.RaisedCosineTransmitFilter('FilterSpanInSymbols', 6);
    fvtool(rctFilt, 'Analysis', 'impulse')
    x = 2*randi([0 1], 96, 1) - 1;
    y = rctFilt(x); 
    plot(y); grid on;
    title('Transmitted Waveform'); 
    xlabel('Signal Index'); ylabel('Amplitude');