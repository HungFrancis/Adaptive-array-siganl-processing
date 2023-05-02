clc;clear;close all
data = load('brain.mat');
im = imresize(abs(data.im), [64,64]);
W = dbwavf('db32');
mask_vardens = fspecial('gaussian',64,50);
DATA = fft2(im).*mask_vardens;
im_cs = ifft2(DATA); % initial value
figure;
for iter=1:20
im_cs = W'*(soft_thresholding(W*im_cs,0.025));
im_cs = ifft2(fft2(im_cs).*(1-mask_vardens) + DATA); % Data consistency
end
imshow(abs(im_cs),[]);

im_reconstruct = pocs2d(mask_vardens*fft2(im).',mask_vardens,W,0.025,1e-10);
figure
imshow(abs(im_reconstruct),[])

function [x_hat_new] = pocs2d(X, A, w,lambda, eps)
x_hat_new = 0;
x_hat_prev = 0;
Fx_hat_prev = A'*(X);
i = 0;
X0 = Fx_hat_prev;
while(i==0 || norm(x_hat_prev-x_hat_new,2)>eps)
    x_hat_prev = x_hat_new;
    x_hat_new = ifft2(Fx_hat_prev);
    x_hat_new = w'*soft_thresholding(w*x_hat_new,lambda);
    Fx_hat_new = fft2(x_hat_new);
    Fx_hat_new(X0 ~= 0) = X0(X0 ~= 0);
    i=i+1;
end
end