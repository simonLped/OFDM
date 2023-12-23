function [errornumber,biterrorrate] = ber(dataIn,dataOut)

            [errornumber,biterrorrate] = biterr(dataIn,dataOut);

end


