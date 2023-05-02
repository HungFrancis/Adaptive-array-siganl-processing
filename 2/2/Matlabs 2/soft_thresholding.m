function [x_thr] = soft_thresholding(y, lambda)
% Write one implementation that works both for y being real and complex.
% x_thr =  wthresh(y,'h',lambda);
x_thr = zeros(length(y),1);
x_thr = (abs(y)-lambda) .* exp(1i*angle(y)) .*(abs(y)>lambda);
% for i = 1:length(y)
%     if abs(y(i)) <= lambda
%         x_thr(i) = 0;
%     else
%         x_thr(i) = ((abs(y(i))-lambda)/abs(y(i)))*y(i);
%     end
end
