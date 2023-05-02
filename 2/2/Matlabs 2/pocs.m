function [x_hat_new] = pocs(X, A, lambda, eps)
% X is the (non-zero-filled) undersampled Fourier measurement. 
% A is the subsampling matrix
% Initialize
i = 0;
x_hat_new = 0;
x_hat_prev = 0;
X0 = A'*X;
Fx_hat_prev = X0;
    while (i == 0 || norm(x_hat_new-x_hat_prev)>eps)
        x_hat_prev = x_hat_new;
        % Soft-thresholding
        x_hat_new = soft_thresholding(ifft(Fx_hat_prev),lambda);
        Fx_hat_new = fft(x_hat_new);
        % Data consistency step
        Fx_hat_new(X0 ~= 0) = X0(X0 ~= 0);
        % 
        Fx_hat_prev = Fx_hat_new;
        i = i+1
        % Use this plot for debugging, but comment it out if you want quicker
        % running.
%         subplot(2,1,2)
%         stem(abs(x_hat_new))
%         title('x_{hat}')
%         pause(.1)
    end
end

% function [x_hat_new] = pocs(X, A, lambda, eps)
% i=0; x_hat_new=0; x_hat_prev=0; 
% X0 = A'*X; 
% X_hat = X0;
% while (i == 0 || norm(x_hat_new-x_hat_prev)>eps)
% x_hat_prev = x_hat_new;
% x_hat_new = soft_thresholding(ifft(X_hat),lambda);
% X_hat = fft(x_hat_new);
% X_hat(X0~=0) = X0(X0~=0);
% i = i+1;
% end
% end