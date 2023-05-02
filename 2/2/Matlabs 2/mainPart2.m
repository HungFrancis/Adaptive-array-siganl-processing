close all;
clear all;
clc;
%% Question 17: Create sampling matrix
M = 1024;
tic
% In each row, select a random column
idx = (0:1023)*4096 + randperm(64*64,M);
A_1d = zeros(1,M*64*64);
A_1d(idx) = 1;
% Reshape by 4096*M, and transpose
A = reshape(A_1d,[64*64,M])';
toc
%%%%%%%%%% TRY %%%%%%%%%%%%%%%
%%%%%%%%%% TRY %%%%%%%%%%%%%%%
% Uncomment to sanity check
unique(sum(A,1));
unique(sum(A,2));
figure
imshow(A,[])
%%
M = 1024;
tic
A = [eye(M) zeros(M,64^2-M)]; % create matrix A
A = A(:,randperm(64^2)); % randomize matrix A
toc
%% Question 18: IFFT reconstruction
% Subsample the image
im = phantom('Modified Shepp-Logan', 64);
% im_reshape = reshape(im,1,[]);
im_reshape = im(:);
y = (A'*A)*(fft(im_reshape));
% Reconstruction via inverse fourier transform
invReconstruction = 64^2/M*(ifft2(reshape(y,64,64)'))/sqrt(M);
figure()
subplot(1,3,1)
imshow(im, [])
title('Original image')
subplot(1,3,2)
imshow(abs(invReconstruction), [])
title(sprintf('IFFT of zero-filled \nsubsampled k-space'))
%%
psi = dftmtx(4096);
figure
imshow((A*psi));

%% Question 19: TV based reconstruction 
N = 64; % x and y dimension of the image
F = dftmtx(N);
lambda = 0.1;
% Matrix with diagonal:-1, 1 element above diagonal is 1
% Dx(M,N) = im(M,N) - im(M-1,N)
% Dy(M,N) = im(M,N) - im(M,N-1)
v1 = -1*ones(N,1);v2 = ones(N-1,1);
diff_ops = diag(v1,0)+diag(v2,1);
% Measurment data y
y = A*(reshape(fft2(im),[],1));
% Optimizer
cvx_begin
    cvx_precision low
    variable imhat(N,N)
      Dx_imhat = imhat*diff_ops;
      Dy_imhat = diff_ops* imhat;
    minimize(norm(A*reshape(F*imhat*F,[],1)-y) + lambda * (norm(Dx_imhat(:),1) + norm(Dy_imhat(:),1)) );
cvx_end
TVreconPhantom = imhat;
% figure()
% subplot(1,3,1)
% imshow(im, [])
% title('Original image')
subplot(1,3,3)
imshow(abs(TVreconPhantom), [])
title(sprintf('TV-domain reconstruction'))
%% Question 20: IFFT vs TV-based MRI reconstruction
% Brain data
data = load('brain.mat');
im = imresize(abs(data.im), [64,64]);
% IFFT reconstruction
% y =  (A'*A)*(fft(reshape(im,[],1)));
y =  (A'*A)*(fft(im(:)));
invReconstructionMRI = 64^2/M*ifft2(reshape(y,64,64)');
F = dftmtx(64);N=64;
%%
figure
imshow(abs(invReconstructionMRI),[])
%% TV-based reconstruction
lambda = 0.1;
F_X_hat_F = reshape(fft2(im),[],1);
y = A*(F_X_hat_F);
cvx_begin
    cvx_precision low
    variable imhat(N,N)
    Dx_imhat =  imhat*diff_ops;;
    Dy_imhat =  diff_ops* imhat;
    F_imhat = F*imhat*F;
    minimize(norm(A*F_imhat(:)-y,2) + lambda * (norm(Dx_imhat(:),1) + norm(Dy_imhat(:),1)) );
cvx_end
TVreconMRI = imhat;

%% Plot the original image and the two reconstructions
figure()
subplot(1,3,1)
imshow(im,[])
title('Original image')
subplot(1,3,2)
imshow(invReconstructionMRI,[])
title(sprintf('IFFT of zero-filled \nsubsampled k-space'))
subplot(1,3,3)
imshow((TVreconMRI),[])
title('TV reconstruction')
% 
%%
figure()
imshow((diff_ops*imhat))

%%
im = phantom('Modified Shepp-Logan', 64); 
data = load('brain.mat');
M = 1024;
tic
A = [eye(M) zeros(M,64^2-M)]; % create matrix A
A = A(:,randperm(64^2)); % randomize matrix A
im = imresize(abs(data.im), [64,64]);
im_k=fft2(im);
y = A*im_k(:); % sparse measurement in k-space
eps = 1e-10; lambda = 0.01; % POCS variables
x_hat_pocs = pocs2d(y, A, lambda, eps); % see function below
pocsReconstruction = reshape(x_hat_pocs,64,64); % reshape to image
figure()
imshow(pocsReconstruction,[])
%% 
data = load('brain.mat');
im = imresize(abs(data.im), [64,64]);
w = dbwavf('db32');
gaussMatrix = random('normal',0,1,size(im))
im_Reconstruct_POCS2d = pocs2d(gaussMatrix*fft2(im),gaussMatrix,w, 0.025,1e-10);
figure
imshow(im_Reconstruct_POCS2d,[])
function [x_hat_new] = pocs2d(X, A, w,lambda, eps)
i=0; 
x_hat_new=zeros(64); 
x_hat_prev=zeros(64); 
X0 = reshape(A'*X,64,64); 
X_hat = X0;
while (i == 0 || norm(x_hat_new-x_hat_prev)>eps)
x_hat_prev = x_hat_new;
x_hat_temp = ifft2(X_hat);
Fx_hat_new = w'*soft_thresholding(w*x_hat_temp,lambda);
x_hat_new = fft2(reshape(Fx_hat_new,64,64));
x_hat_new(X0~=0) = X0(X0~=0);
i = i+1;
norm(x_hat_new(:)-x_hat_prev(:))
end
end
