function [data_out, series_msg_length, mask, zeros_to_append, trainingblock] = ofdm_mod(EQ_mode, dataIN, P,L, H, BW, TB_size, data_size)

    
    %----------------------------------------------------------------------
    %-----------------Bit loading mask from input H estimate---------------
    %----------------------------------------------------------------------
    mask = 69*ones(P+2,1);
    BW_ind = 100 - BW;
    if BW_ind ~= 0
        H_sorted=sort(abs(H));
        H_min = H_sorted(round(length(H)*0.01*BW_ind));
    else
        H_min = 0;
    end

    for i=1:P*2+2
        if abs(H(i)) >= H_min
            mask(i) = 1;
        else
            mask(i) = 0;
        end
    end
    subcarrier_mask = mask(2:P+1);
    rows_sent = sum(subcarrier_mask);


    %----------------------------------------------------------------------
    %-----------------Length of signal-------------------------------------
    %----------------------------------------------------------------------
    length_data = length(dataIN);
    length_data_update = length_data;
    
    while rem(length_data_update, rows_sent) ~= 0
        length_data_update = length_data_update + 1;
    end

    zeros_to_append = length_data_update - length_data;
    dataIN = [dataIN ; zeros(zeros_to_append,1)];
    dataIN_array = reshape(dataIN,rows_sent,[]);

    subcarriers = zeros(P,size(dataIN_array,2));
    k = 1;
    for i=1:P
        if subcarrier_mask(i)
            subcarriers(i,:) = dataIN_array(k,:);
            k=k+1;
        end
    end

    %----------------------------------------------------------------------
    %-----------------Packet of image -------------------------------------
    %----------------------------------------------------------------------
    temp_c = size(subcarriers,2);
    con_subcarriers = flip(conj(subcarriers),1); % conjugated mirrored columns 
    packet = [zeros(1,temp_c);
              subcarriers;
              zeros(1,temp_c);
              con_subcarriers];
    

    if EQ_mode == 0
        %----------------------------------------------------------------------
        %-----------------Training packets & cyclic prefix---------------------
        %----------------------------------------------------------------------
        trainingblock = packet(1:end,1:TB_size);
        number_of_packets = ceil(size(packet,2)/data_size);
        
        n = 1* size(packet, 1); % freq-bins = rows of packet*2
        r = size(packet,1); %find rows of ifft_array
        ifft_array = [];
        packet_sent = [];
        for i=1:number_of_packets
            
            end_pos_data = data_size*i;
            if i == number_of_packets
                end_pos_data = size(packet,2);
            end
            data_block = packet(:, data_size*i-data_size+1:end_pos_data);
            packet_sent = [packet_sent,trainingblock,data_block];
        end
        packet = packet_sent;
        ifft_array = ifft(packet,n,1); %ifft of each column
        cyclic_prefix = ifft_array((r-L+1):r,:); % cyclic prefix
        parallel_array = [cyclic_prefix;    % parallel data
                          ifft_array];
        
        %----------------------------------------------------------------------
        %-----------------Parallel to series-----------------------------------
        %----------------------------------------------------------------------
        [r,c] = size(parallel_array);
        data_out = [];
        for i=1:number_of_packets
            
            size1 = data_size+TB_size;
            end_pos_data = size1*i;
            if i == number_of_packets
                end_pos_data = size(packet,2);
            end
            TB_DB_packet = parallel_array(:, size1*i-size1+1:end_pos_data);
            data_out = [data_out , reshape(TB_DB_packet, 1, [])];
        end
        
            series_msg_length = r*c;
    elseif EQ_mode == 1
        
        trainingblock = packet(1:end,1:TB_size);
        
        packet = [trainingblock, packet];
        n = 1* size(packet, 1); % freq-bins = rows of packet*2
        r = size(packet,1); %find rows of ifft_array

        ifft_array = ifft(packet,n,1); %ifft of each column
        cyclic_prefix = ifft_array((r-L+1):r,:); % cyclic prefix
        parallel_array = [cyclic_prefix;    % parallel data
                          ifft_array];
        
        %----------------------------------------------------------------------
        %-----------------Parallel to series-----------------------------------
        %----------------------------------------------------------------------
        [r,c] = size(parallel_array);
        series_msg_length = r*c;
        data_out = reshape(parallel_array, 1, series_msg_length);
    end

end