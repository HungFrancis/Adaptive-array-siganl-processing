clear all;
close all;
clc;
set(0,'defaultfigurecolor','w')
%% Question 3 
% Generate signal
x = [[1:5]/5 zeros(1,128-5)];
x = x(randperm(128));
y = x + 0.05*randn(1,128);
lambdas = [0.1,1,3];
figure
for i = 1:length(lambdas)
    subplot(1,3,i)
    stem(y,'filled');
    hold on
    x_hat = 1/(1+lambdas(i))*y;
    stem(x_hat, 'filled');
    legend('$y$', '$\hat{x}$','Interpreter','latex')
    grid on;
    title(sprintf('Estimate of x for \\lambda = %.1d',lambdas(i)))
end
%% Question 5: Soft-thresholding function

% Applying the function on real values
lambda = 2;
x_thr = soft_thresholding(-10:10,lambda);
figure()
hold on;
stem(-10:10, -10:10,'filled');
stem(-10:10,x_thr,'filled');
grid on
title('Soft-thresholding real values for lambda = 2')
legend('x', 'soft\_thr(x)', 'Location', 'northwest')

% Applying the function on complex values

% A list of complex numbers of which the magnitude ranges between
% 0 and 10. Note; we use a random phase of 1 here.
x_complex = (0:10).*exp(1i);
x_thr = soft_thresholding(x_complex,lambda);
figure()
hold on;
stem(0:10, 0:10,'filled');
stem(0:10,abs(x_thr),'filled');
grid on
title('Soft-thresholding complex values for \lambda = 2')
legend('x', 'soft\_thr(x)', 'Location', 'northwest')

%% Question 6: Soft-thresholding on signal y
figure
lambdas = [0.01, 0.05, 0.1, 0.2]
for i = 1:length(lambdas)
    subplot(2,2,i)
    stem(y, 'filled');
    hold on
    x_hat = soft_thresholding(y,lambdas(i));
    stem((x_hat),'filled');
    legend('$y$', '$\hat{x}$','Interpreter','latex') %This latex syntax requires Matlab R2018b or newer.
    grid on;
    title(sprintf('Estimate of x for \\lambda = %.1d',lambdas(i)))
    axis tight
end

%% Question 10: Uniform subsampling
y = x ;
X = fft(x);
% figure
% plot(abs(X))
% Uniform subsampling matrix 
Au = zeros(size(X,2));
space = 4;
for i = 1:floor(size(X,2)/space)
    Au(i*space,i*space) = 1;
end
% for i = 1:64
%     Au(i,i)=1;
% end
% Au = eye(128)
% Use .' to transpose X if needed for matrix multiplication.
% Since X contains complex numbers, if you'd just apply ' (non-element-wise
% transpose), you calculate the complex conjugate. 
% Zero-filled undersampled X
Xu = 4*(Au'*Au)*(X.');
% Plot here your Xu so check that it takes the form as you expect!
figure
stem(Xu)
% Fourier inverse reconstruction
x_hat_u = (ifft(Xu));

figure()
stem(x, 'fill');
hold on;
stem(abs(x_hat_u));
hold off
title('Reconstruction from uniform undersampling')
legend('original signal x', 'reconstructed signal $\hat{x_u}$', 'Interpreter', 'Latex')
ylabel('Amplitude');xlabel('Number of sample')
axis tight;ylim([-0.1, 1.3]);grid on;

%% Question 11: Random subsampling
% Find random sampling matrix Ar
Ar = zeros(size(X,2));
idx = randperm(128,32)
for i = 1:length(idx) 
    Ar(idx(i),idx(i)) = 1;
end
% figure
% imshow(Ar)
% Zero-filled subsampled X
Xr = 4.*(Ar*Ar')*(X.'); 
figure
stem(Xr);
figure
stem(4.*(Ar')*(X.'))
% Plot here your Xr so check that it takes the form as you expect!

% Reconstruction via inverse Fourier transform
x_hat_r =(ifft(Xr));
% figure()
% stem(x, 'fill');
% hold on;
% stem(abs(x_hat_r));
% title('Reconstruction from random undersampling')
% legend('original signal x', 'undersampled signal $\hat{x_r}$', 'Interpreter', 'Latex')
% axis tight; ylim([-0.1, 1.3])
% ylabel('Amplitude');xlabel('Number of sample')
% grid on;

%% Question 12: POCS
eps = 1e-10;
x_hat_pocs_res = zeros(3,length(x));
lambdas = [0.001, 0.01, 0.1];

for i = 1:length(lambdas)
%Use the plotting for debugging purposes
    figure()
    subplot(2,1,1)
    stem(x)
    title('original x')
    x_hat_pocs = pocs(Ar*fft(x).',Ar,lambdas(i),eps);
    % Save the result for this lambda
    x_hat_pocs_res(i,:) = x_hat_pocs;
end

% Plot the results for all lambdas
figure()
for i = 1:length(lambdas)
    subplot(1,3,i)
    stem(x);
    hold on;
    stem(abs(x_hat_pocs_res(i,:)));
    hold off;
    title(sprintf('\\lambda = %.3d',lambdas(i)))
    legend('original signal x', 'Reconstructed signal')
    axis tight;
    ylim([-0.1, 1.3])
    grid on;    
end
sgtitle('POCS random undersampling')

%% Question 13: ISTA
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
beta = 1e-2;
eps = 1e-10;
x_hat_ista_res = zeros(3,length(x));
lambdas = [1e-5, 1e-4, 1e-3];
% lambdas = [0.001, 0.01, 0.1];
for i = 1:length(lambdas)
    %Use the plotting for debugging purposes
    figure()
    subplot(2,1,1)
    stem(x)
    title('original x')
    x_hat_ista = ISTA(Ar*fft(x).',Ar,lambdas(i),beta,eps);
    
    % Save the result for this lambda
    x_hat_ista_res(i,:) = x_hat_ista;

end

% Plot the results for all lambdas
figure()
for i = 1:length(lambdas)
    subplot(1,3,i)
    stem(x);
    hold on;
    stem(abs(x_hat_ista_res(i,:)));
    hold off;
    title(sprintf('\\lambda = %.3d',lambdas(i)))
    legend('original signal x', 'Reconstructed signal')
    axis tight;ylim([-0.1, 1.3])
    grid on;
end
sgtitle('ISTA random undersampling')


