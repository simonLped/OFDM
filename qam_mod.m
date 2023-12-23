function [dataMod, k] = qam_mod(dataBitsIn, M)
            
        k = log2(M);
        dataSymbolsIn = bit2int(dataBitsIn,k);
        dataMod = qammod(dataSymbolsIn,M, 'bin', 'UnitAveragePower', true);

    

end