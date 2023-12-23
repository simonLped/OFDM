function [data_out, series_msg_length] = ofdm_mod_H_est(packet, P,L)
    %L=prefix length

    
    n = 1* size(packet, 1); % freq-bins = rows of packet*2
    
    
    
    r = size(packet,1); %find rows of ifft_array

    ifft_array = ifft(packet,n,1); %ifft of each column



    cyclic_prefix = ifft_array((r-L+1):r,:); % cyclic prefix
    
    parallel_array = [cyclic_prefix;    % parallel data
                      ifft_array];
    
    
    [r,c] = size(parallel_array);

    len_data_out = r*c;

    data_out = reshape(parallel_array, 1, len_data_out);

    series_msg_length = len_data_out;

end