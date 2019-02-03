Spaced_out_transmitted_data = [];           
            
            
 for Frame_iterator= 1:numFrames
    Symbol_in_frame_iterator = 1;
    %Getting all of filteredTXData and noisydata in 1 array
    for All_sample_iterator= 1:(samplesPerSymbol*frameSize)
        longTx=[longTx, allFilteredTXData(All_sample_iterator,Frame_iterator)];
        longNoise = [longNoise, allnoisyData(All_sample_iterator,Frame_iterator)];
    end 
    %creating a array with properly spaced filteredData points to plot with
    %filterTX data
    for Symbol_in_frame_iterator= 1:frameSize
        hjk = 1;
        if Symbol_in_frame_iterator ==1
            Spaced_out_transmitted_data(((Frame_iterator-1)*frameSize*samplesPerSymbol)+1)...
                = modulatedData(1+(samplesPerSymbol*(Frame_iterator-1)));
            for hjk =1:(samplesPerSymbol-1)
                Spaced_out_transmitted_data = [Spaced_out_transmitted_data,0];
            end
        else
            Spaced_out_transmitted_data((samplesPerSymbol*(Symbol_in_frame_iterator-1)+...
                ((Frame_iterator-1)*frameSize*samplesPerSymbol)+1))...
                = modulatedData(Symbol_in_frame_iterator+(samplesPerSymbol*(Frame_iterator-1)));
            for hjk =1:(samplesPerSymbol-1)
                Spaced_out_transmitted_data = [Spaced_out_transmitted_data,0];
            end
        end 
    end
end
