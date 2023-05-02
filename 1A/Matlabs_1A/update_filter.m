function [ filter ] = update_filter( filter,e )
x=filter.x_delayed;
w_old = filter.w;
r=filter.r;
filter_type=filter.type;
rx_inv_est_old = filter.inv_Rx_RLS;
rex_est_old = filter.rex_corr;
filter.w=w_old; %default output

%% A1 scenario 1:f
if strcmpi(filter_type,'SGD')
    %implement the SGD update rule here
    alpha=filter.adaptation_constant;
    Rx = [2,-1;-1,2];
    rex = [0;3];
    filter.w= w_old+ 2 * alpha *(rex-Rx*w_old);
end

%% A1 scenario 1:i
if strcmpi(filter_type,'Newton')
    %implement the Newton update rule here
    alpha=filter.adaptation_constant;
    Rx = [2,-1;-1,2];
    rex = [0;3];
    inv_Rx = inv(Rx);
    filter.w= w_old + 2* alpha*inv_Rx*(rex-Rx*w_old)
end

%% A1 scenario 2:a
if strcmpi(filter_type,'LMS')
    alpha=filter.adaptation_constant;
%     tic
    %implement the LMS update rule here
    filter.w= w_old + 2*alpha*x*r % O(n^2)
%     toc
end

%% A1 scenario 2:b
if strcmpi(filter_type,'NLMS')
    %implement the NLMS update rule here
    beta = 0.8;
%     tic
    alpha=filter.adaptation_constant;
    filter.var_est = beta*filter.var_est+(1-beta)*((x'*x)/length(x));%(o(n))
    filter.w=w_old+2*(alpha/filter.var_est)*x*r;% O(n^2)    
%     toc
end

%% A1 scenario 2:d
if strcmpi(filter_type,'RLS')
     %implement the RLS update rule here
    gamma = filter.adaptation_constant;
%     tic
    G = (rx_inv_est_old*x)/(gamma^2+x'*rx_inv_est_old*x);%(n*1)*(n*1)/ (1+(1*n)*(n*n)*(n*1))
    filter.inv_Rx_RLS = gamma^-2*(rx_inv_est_old-G*x'*rx_inv_est_old); % (N*1)*(1*N)*(N*N) 
    filter.rex_corr = (gamma^2*rex_est_old)+x*e; % (N*1) * (1*1)
    filter.w = filter.inv_Rx_RLS*filter.rex_corr;% (N*N) * (N*1)
%     toc
end % 3 O(n^2)

%% A1 scenario 2:e
if strcmpi(filter_type,'FDAF')
    alpha=filter.adaptation_constant;
    beta = 0.5;
%     tic
    X = filter.F*x; % (N*N) * (N*1) = O(n^2) multiply operation, O(n^2) multply operation
    filter.P_est=beta*filter.P_est+(1-beta)*X*conj(X)'/length(x); % (N*1) * (1*N) = N*N
    W_old = filter.F_inv*w_old; % (N*N) * (N*1) = O(n^2) multiply operation, O(n^2) multply operation
    W_new = W_old+2*alpha*1./diag(filter.P_est).*conj(X)*r;  % (N*N) * (N*1) *1
    filter.w=filter.F*W_new;% (N*N) * (N*1)
%     toc
end

end

