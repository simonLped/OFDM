% Construct received ofdm_stream using a different channel for each
% training/data block
% Input arguments
%   Tx: signal to be transmitted over channel
%   Lt: number of training packets per training/data block
%   Ld: number of data packets per training/data block
%   N: used DFT-size
%   L: used cylcic prefix length
% Output arguments
%   Rx: received ofdm_stream over changing channel

function Rx = simulate_channel(Tx, Lt, Ld, N, L)

Rx = zeros(size(Tx));

P = ceil(length(Tx)/((Lt+Ld)*(N+L)));
step = (Lt+Ld)*(N+L);

for p=1:P-1
    h = randn(L,1);
    Rx( (p-1)*step+1 : p*step ) = fftfilt(h, Tx( (p-1)*step+1 : p*step ));
end

h = randn(L,1);
Rx( (P-1)*step+1 : end ) = fftfilt(h, Tx( (P-1)*step+1 : end ));

end
