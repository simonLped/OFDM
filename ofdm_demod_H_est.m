function [H_est] = ofdm_demod_H_est(recived_signal, P, L, series_msg_length,trainingblock)


    recived_signal = recived_signal(:, 1:series_msg_length); %cut signal extended by channel to rx signal size
    r1 = ((P*2)+2)*1; %rows of ifft in ofdm_mod, P = subcarriers, *2 = flipped conjugated, +3 = zero padding
    r_with_cp = r1 + L; %rows with cyclic prefix
    recvd_signal_parallel = reshape(recived_signal,r_with_cp, []); % parallel of recived_signal
    
    fft_array = fft(recvd_signal_parallel((L+1):r_with_cp,:)); % remove cyclic prefix

    % fft_array_with_TB = fft(signal_without_cp_parallel,r1,1).*mask;   %fft of columns
    
    H_est = zeros(size(fft_array,1),1);

    % for i = 1:size(fft_array,1)  % subcarriers
    %     H_est(i) = mean(conj(fft_array(i,:)'./trainingblock(i,:)'));
    % end

    for r = 1:size(fft_array,1)  % subcarriers
        H_est(r) = mean(conj(fft_array(r,:)'./trainingblock(r,:)'));
    end
end