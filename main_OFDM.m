% OFDM LS packet based & LMS adaptive
load("channel.mat") % simulate fast varying channel
load("channel_constant_time.mat") % simulate constant channel
%--------------------------------------------------------------------------
%----------------------------READ ME---------------------------------------
%--------------------------------------------------------------------------

% 1) Chose adaptive LMS or packet based LS with EQ_mode=1 or EQ_mode=0

% 2) Chose channel: simulated fast-varying: channel real_channel_on=0
%                   over speaker: real_channel_on=1
%                   simulated packet-varying: channel real_channel_on=2
%             simulated constant channel: & align-pulse: real_channel_on= 3

% 3) Chose P=subcarriers, L=cyclic prefix, BW=Bandwidth (on-off-bitloading)
%          TB_size = training frames per packet (used if EQ_mode=1)
%          Else TB_size = LS starting values for LMS
%          TB_size_H_est = training frames used for on-off-bit-loading
%          snr = awgn noise added to signal
%          margin = margin for aligned signal
%          M = QAM size

%--------------------------------------------------------------------------
%----------------- Setup and parameters------------------------------------
%--------------------------------------------------------------------------

[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
M = 16; % QAM size
[qam_symbols, k] = qam_mod(bitStream, M);

P = 512; % subcarriers
L = 300; % cyclic prefix > length of channel
BW = 50; % Bit loading 
TB_size = 10; % training frames
data_size = 10; % Data frames per packet
TB_size_H_est = 10; % number of frames to estimate channel for bit loading
fs = 16000; % sample rate (fixed)
snr = 20;
mu = 0.5;
alpha = 1;
margin = 71; % correlation margin
EQ_mode = 1; % packet-based-LS = 0, adaptive LMS = 1, 
real_channel_on = 0; % simulated_channel_milestone = 0, speaker = 1, simulated_channel = 2, simulated constant = 3
plot_error = 0; % plot error for LMS (slow, but cool), can be turned off (0)

%--------------------------------------------------------------------------
%------------------------------------Training data H_est-------------------
%--------------------------------------------------------------------------

H_est = ones(P*2+2,1); % Needed as input for next line so that BW = 100
[~, ~, ~, ~, trainingblock] = ofdm_mod(EQ_mode, qam_symbols, P,L, H_est, 100,TB_size_H_est,5);
trainingblock_H_est = trainingblock;

%--------------------------------------------------------------------------
%-------------------------------Align pulse--------------------------------
%--------------------------------------------------------------------------

% Parameters for the sine wave align signal
frequency = 0.01; % Frequency in Hz
duration = 10000; % Duration in seconds
t = 0:1:duration;
zM_Nsubc = 0.3*sin(2 * pi * frequency * t); % Sin wave
zM_Nsubc = [zM_Nsubc, zeros(1,1000)]; % Needs to have zeros longer than channel

%--------------------------------------------------------------------------
%-------------------Estimate Channel For Bit Loading-----------------------
%--------------------------------------------------------------------------
% Modified OFDM to make training data estimate channel used for bitloading
[ofdm_data_H_est, series_msg_length_H_est] = ofdm_mod_H_est(trainingblock_H_est, P,L);

if real_channel_on == 1
    ofdm_data_H_est = [zM_Nsubc, ofdm_data_H_est];
    ofdm_data_H_est = awgn(ofdm_data_H_est,snr, 'measured');
    
    [simin, nbsecs, fs] = initparams(ofdm_data_H_est, fs);
    
    sim('recplay'); % simulation line
    out = simout.signals.values; % recieved signal
    % % % soundsc(out,fs); % play the sound back
    signal_over_channel_h_est = out';
    [C21,lag21] = xcorr(signal_over_channel_h_est, zM_Nsubc);
    C21 = C21/max(C21);
    [M21,I21] = max(C21);
    t21 = lag21(I21);
    signal_over_channel_h_est = signal_over_channel_h_est(t21+length(zM_Nsubc)-margin:end);

elseif real_channel_on == 0
    signal_over_channel_h_est = simulate_channel_milestone(ofdm_data_H_est, TB_size_H_est, 0, P*2+2, L, snr);
elseif real_channel_on == 2
    signal_over_channel_h_est = simulate_channel(ofdm_data_H_est, TB_size_H_est, 0, P*2+2, L);
else
    ofdm_data_H_est = [zM_Nsubc, ofdm_data_H_est];
    ofdm_data_H_est = awgn(ofdm_data_H_est,snr, 'measured');
    snr_H_est = 10000;
    signal_with_noise_H_est = awgn(ofdm_data_H_est,snr, 'measured');
    signal_over_channel_h_est = conv(signal_with_noise_H_est,h);
    [C21,lag21] = xcorr(signal_over_channel_h_est, zM_Nsubc);
    C21 = C21/max(C21);
    [M21,I21] = max(C21);
    t21 = lag21(I21);
    signal_over_channel_h_est = signal_over_channel_h_est(t21+length(zM_Nsubc)-margin:end);
end
[H_est] = ofdm_demod_H_est(signal_over_channel_h_est, P, L, series_msg_length_H_est, trainingblock_H_est);

%--------------------------------------------------------------------------
%----------------------Modulating And Sending Data-------------------------
%--------------------------------------------------------------------------

[ofdm_data, series_msg_length, mask, zeros_to_append, trainingblock] = ofdm_mod(EQ_mode, qam_symbols, P,L, H_est, BW,TB_size,data_size);

if real_channel_on == 1
    ofdm_data = [zM_Nsubc, ofdm_data];
    ofdm_data_Tx = awgn(ofdm_data,snr, 'measured');
    [simin, nbsecs, fs] = initparams(ofdm_data_Tx, fs);
    sim('recplay'); % simulation line
    out = simout.signals.values; % recieved signal
    % % soundsc(out,fs); % play the sound back
    signal_over_channel = out';

    [C21,lag21] = xcorr(signal_over_channel, zM_Nsubc);
    C21 = C21/max(C21);
    [M21,I21] = max(C21);
    t21 = lag21(I21);
    
    signal_over_channel = signal_over_channel(t21+length(zM_Nsubc)-margin:end);

elseif real_channel_on == 0
    signal_over_channel = simulate_channel_milestone(ofdm_data, TB_size, data_size, P*2+2, L,snr);
elseif real_channel_on == 2
    signal_over_channel = simulate_channel(ofdm_data, TB_size, data_size, P*2+2, L);
else
    ofdm_data = [zM_Nsubc, ofdm_data];
    signal_with_noise = awgn(ofdm_data,snr, 'measured');
    signal_over_channel = conv(signal_with_noise,h);
    [C21,lag21] = xcorr(signal_over_channel, zM_Nsubc);
    C21 = C21/max(C21);
    [M21,I21] = max(C21);
    t21 = lag21(I21);
    signal_over_channel = signal_over_channel(t21+length(zM_Nsubc)-margin:end);
end

%--------------------------------------------------------------------------
%----------------------Demodulating Recieved Data--------------------------
%--------------------------------------------------------------------------


if EQ_mode == 1
    figure(2)
    colormap(colorMap); image(imageData); axis image; title('Transmitted image'); drawnow;
elseif EQ_mode == 0
    subplot(2,2,4); colormap(colorMap); image(imageData); axis image; title('Transmitted image'); drawnow;
end


[recived_ofdm] = ofdm_demod(EQ_mode, signal_over_channel,P, L, series_msg_length, mask, zeros_to_append,M, imageSize, bitsPerPixel, colorMap, trainingblock, TB_size, data_size,mu,alpha,plot_error);
[er,bit1] = ber(bitStream,recived_ofdm); 
disp(bit1)




