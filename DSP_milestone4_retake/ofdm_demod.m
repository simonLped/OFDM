function [bit_stream_out1] = ofdm_demod(EQ_mode, recived_signal, P, L, series_msg_length, mask, zeros_to_append, M, imageSize, bitsPerPixel, colorMap, trainingblock, TB_size, data_size)
    
    %----------------------------------------------------------------------
    %----------------------Cutting serial signal---------------------------
    %----------------------------------------------------------------------
    recived_signal = recived_signal(:, 1:series_msg_length); %cut signal extended by channel to rx signal size
    r1 = ((P*2)+2)*1; %rows of ifft in ofdm_mod, P = subcarriers, *2 = flipped conjugated, +3 = zero padding
    r_with_cp = r1 + L; %rows with cyclic prefix
    number_of_columns = series_msg_length/r_with_cp;

    if EQ_mode == 0
        number_of_packets = ceil(number_of_columns/(data_size+TB_size));
        x = number_of_packets*data_size;
        numer_training_frames = TB_size*number_of_packets;
        number_data_frames = number_of_columns - numer_training_frames;
        size_last_packet = x - number_data_frames; 
        recvd_signal_parallel = [];
    
        %----------------------------------------------------------------------
        %----------------------Series to parallel------------------------------
        %----------------------------------------------------------------------
        for i=1:number_of_packets
            end_pos = r_with_cp*(data_size+TB_size)*i;
            size1 = (TB_size + data_size)*r_with_cp;
            if i == number_of_packets
                end_pos = series_msg_length;
                size1 = (TB_size + data_size-size_last_packet)*r_with_cp;
            end
            series_img_fraction = recived_signal(end_pos-size1+1:end_pos);
            parallel_img_fraction = reshape(series_img_fraction, r_with_cp, []);
            recvd_signal_parallel = [recvd_signal_parallel, parallel_img_fraction ];
        end
        signal_without_cp_parallel = recvd_signal_parallel((L+1):r_with_cp,:); % remove cyclic prefix
        fft_array_with_TB = fft(signal_without_cp_parallel,r1,1).*mask;   %fft of columns
    
        %----------------------------------------------------------------------
        %----------------------Channel estimate & EQ---------------------------
        %----------------------------------------------------------------------
        training_data = [];
        image_data = [];
        image_data_EQ = [];
        H_lsr = zeros(size(fft_array_with_TB,1),number_of_packets);
        for i=1:number_of_packets
            end_pos_data = (data_size+TB_size)*i;
            if i == number_of_packets
                  end_pos_data = size(fft_array_with_TB,2);
            end
            packet_block = fft_array_with_TB(:, end_pos_data - (data_size+TB_size)+1:end_pos_data);
            if i == number_of_packets
                packet_block = fft_array_with_TB(:, end_pos_data - (data_size+TB_size) + size_last_packet +1:end_pos_data);
            end
            trainingblock_rx = packet_block(:,1:TB_size);
            training_data = [training_data, trainingblock_rx];
            image_packet = packet_block(:,TB_size+1:end);
            image_data = [image_data, image_packet];
            for r = 1:size(trainingblock_rx,1)  % subcarriers
                H_lsr(r,i) = conj(lsqr(trainingblock(r,:)',trainingblock_rx(r,:)'));
                clc;
            end
            image_data_EQ = [image_data_EQ, image_packet./H_lsr(:,i)];
        end

    elseif EQ_mode == 1
                
        recvd_signal_parallel = reshape(recived_signal,r_with_cp, []); % parallel of recived_signal
        signal_without_cp_parallel = recvd_signal_parallel((L+1):r_with_cp,:); % remove cyclic prefix
        fft_array_with_TB = fft(signal_without_cp_parallel,r1,1).*mask;   %fft of columns
        
        fft_array_TB = fft_array_with_TB(:,1:TB_size);
        fft_array_image = fft_array_with_TB(:,TB_size+1:end);
        
        
        %----------------------------------------------------------------------
        %----------------------Initial channel estimate------------------------
        %----------------------------------------------------------------------

        H_lsr = zeros(size(fft_array_with_TB,1),1);
        for r = 1:size(fft_array_with_TB,1)  % subcarriers
            H_lsr(r) = conj(lsqr(trainingblock(r,:)',fft_array_TB(r,:)'));
            clc;
        end

        %----------------------------------------------------------------------
        %----------------------LMS channel estimate----------------------------
        %----------------------------------------------------------------------
        W_k_initial_vector = (1+ 0.2)./conj(H_lsr);
        
        H_lms = zeros(P,1);
        iter = size(fft_array_image,2);
        mu = 0.5;
        a = 1;
        error_matrix = zeros(P*2+2,iter);
        H_lms(1) = 0;
        H_lms_matrix = zeros(P,number_of_columns - TB_size);
        for l=2:P+1
            if H_lsr(l) ~= 0
                Y_k = fft_array_image(l,:);
                W_k_initial = W_k_initial_vector(l);
                [error, W_k] = LMS2(Y_k,W_k_initial,iter,mu,a,M);
                error_matrix(l,:) = error;
                H_lms(l-1) = conj(1./W_k(end));
                H_lms_matrix(l-1,:) = conj(1./W_k);
            end
        end

        H_lms = [0; H_lms; 0; flip(conj(H_lms))];
        H_lms_matrix = [zeros(1,number_of_columns - TB_size); 
                        H_lms_matrix; 
                        zeros(1,number_of_columns - TB_size);
                        flip(conj(H_lms_matrix),1)];
%         image_data_EQ = fft_array_image./H_lms;
        
        for i =1:size(H_lms_matrix,2)
            image_data_EQ(:,i) = fft_array_image(:,i)./H_lms_matrix(:,i);
        end


        
    end
    
    %----------------------------------------------------------------------
    %-----------------Mask & visualizing------------------------------------
    %----------------------------------------------------------------------
    [~, temp_c] = size(image_data_EQ); % size of rows
    figure(1)
%     subcarriers_bitloading = zeros(rows_sent,temp_c);
    subcarriers_bitloading = zeros(sum(mask)/2,temp_c);
    [r2, c1] = size(subcarriers_bitloading);

    k = 1;
    z=1;
    prev_m=1;
    for m=1:c1
        for i=2:P+1

            if (mod(m,data_size) == 0) && (m~=prev_m) && (EQ_mode==0) && (z~=number_of_packets)
                z=z+1; 
            end
            prev_m = m;

            % Mask for ON-OFF-bitloading
            if mask(i)
                subcarriers_bitloading(k,m) = image_data_EQ(i,m);
                k = k +1;
            end
        end
        k=1;

        % parallel to series
        recived_serial = conj(reshape(subcarriers_bitloading, 1,(c1*r2))');
        recived_serial = recived_serial(1:end-zeros_to_append);

        % QAM-demod
        [bit_stream_out1] = qam_demod(recived_serial, M);

        % Visualize
        imageRx = bitstreamtoimage(bit_stream_out1, imageSize, bitsPerPixel);
        subplot(2,2,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
        if EQ_mode == 0
            subplot(2,2,1); plot(abs(H_lsr(:,z))); drawnow;
        else
            subplot(2,2,1); plot(abs(H_lms_matrix(:,m))); drawnow;
        end
        % subplot(2,2,3); plot(abs(ifft(H_lsr(:,z)))); drawnow;
        
    end
    
end