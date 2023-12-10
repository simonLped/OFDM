% Construct received ofdm_stream using a different channel for each
% OFDM frame
% Input arguments
%   Tx: signal to be transmitted over channel
%   N: used DFT-size
%   L: used cylcic prefix length

% Output arguments
%   Rx: received ofdm_stream over changing channel


function Rx = simulate_channel_milestone(Tx, Lt, Ld, N, L, SNR)
TX_length = length(Tx);
Rx = zeros(size(Tx));

P = ceil(length(Tx)/(N+L)); % number of OFDM frames
step = N + L;

padLength = step*P - TX_length;
Tx = [Tx; zeros(padLength,1)];

load("channel.mat", "H");
for p=1:P
    h = H(p,:);
    sig = fftfilt(h, Tx( (p-1)*step+1 : p*step ));
    Rx( (p-1)*step+1 : p*step ) = awgn(sig,SNR, "measured");
end
Rx = Rx(1:TX_length);

end