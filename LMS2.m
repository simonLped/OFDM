function [error, W_k] = LMS2(Y_k,W_k_initial,iter,mu,a,M)
% Y_k = one row 150x1 recived data one subcarrier
% W_k_initial = initial weights 1026x1 
% iter = length of Y_k = 150



W_k = zeros(1,iter);
W_k(1) = W_k_initial;
error = zeros(1,iter);



for l=1:iter-1

    X_k_hat = conj(W_k(l))*Y_k;
    estimated_bits = qam_demod(X_k_hat, M);
    X_k_hat_decision = qam_mod(estimated_bits, M);
    
    

%     W_k(l+1) = W_k(l) + (mu/(a + abs(Y_k(l+1)))) * Y_k(l+1)* conj(X_k_hat_decision(l+1) - conj(W_k(l))*Y_k(l+1));
    W_k(l+1) = W_k(l) + (mu/(a + ctranspose(Y_k(l+1)) * Y_k(l+1))) * Y_k(l+1) * conj(X_k_hat_decision(l+1) - ctranspose(W_k(l)) * Y_k(l+1));
    
    error(l) = abs(X_k_hat_decision(l+1) - conj(W_k(l))*Y_k(l+1));


end




















