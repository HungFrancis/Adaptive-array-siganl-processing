function [x_hat_new] = ISTA(X, A, lambda, beta, eps)
% X is the (non-zero-filled) undersampled Fourier measurement. 
% A is the subsampling matrix
i = 0;
x_hat_prev = 4.*ifft(A'*X);
x_hat_new = x_hat_prev;
FH = conj(dftmtx(size(X,1)))./size(X,1);
% FH = dftmtx(size(X,1));
    while (i == 0 || norm(x_hat_new - x_hat_prev)>eps)
        x_hat_prev = x_hat_new;
         % Gradient update step
         x_hat_new = x_hat_prev-beta*FH*(A'*(A*fft(x_hat_prev)-X));
%       x_hat_new = x_hat_prev-beta*ifft(A'*(A*fft(x_hat_prev)-X))*size(A,2);
        % Soft-thresholding
        x_hat_new = soft_thresholding(x_hat_new,lambda);
        i = i + 1;
        % Use this plot for debugging, but comment it out if you want quicker
        % running.
%         norm(x_hat_new - x_hat_prev)
%         subplot(2,1,2)
%         stem(x_hat_new)    
%         title('x_{hat}')
%         pause(0.05)
    end
%         subplot(2,1,2)
%         stem(x_hat_new)    
%         title('x_{hat}')
%         pause(0.05)
end
