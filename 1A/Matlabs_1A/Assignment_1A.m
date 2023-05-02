clear all;close all;
set(0,'defaultfigurecolor','w') 
%% generate x1 and x2
nsample=200;
[x e]=generate_input(nsample);

%% set up filter
filterA1=adaptive_filter(2,'SGD',0.01);
% filterA1 = adaptive_filter(2,'Newton',0.3);
% filterA1 = adaptive_filter(300,'LMS',0.01);
% filterA1 = adaptive_filter(300,'NLMS',0.01);
% filterA1 = adaptive_filter(300,'RLS',0.03);
filterA1 = adaptive_filter(300,'FDAF',0.9);
%% perform nsample iterations
tic
for sample=1:nsample
   filterA1=filterA1.filter(x(sample),e(sample)); 
end
toc
% sprintf("W1 last")
% filterA1.w_history(200,1)
% sprintf("W2 last")
% (filterA1.w_history(200,2))
%% plot filter coefficients
hold on
plot(filterA1.w_history(:,1),filterA1.w_history(:,2),'x-');
title(strcat('filter algorithm: ',filterA1.type,',adaptation constant: ',num2str(filterA1.adaptation_constant)))
hold off