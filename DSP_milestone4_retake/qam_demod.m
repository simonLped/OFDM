function [dataOut] = qam_demod(receivedSignal, M)

            k = log2(M);
            dataSymbolsOut = qamdemod(receivedSignal,M,'bin', 'UnitAveragePower', true);
            dataOut = int2bit(dataSymbolsOut,k);
    
end