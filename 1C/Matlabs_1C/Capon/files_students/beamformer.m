function beamformed_data = beamformer(dataset, scan, pw_indices, adaptive)

    time = (0:dataset.N_ax-1).'/dataset.fs;
       
    [X, Z] = meshgrid(scan.x, scan.z);

    beamformed_channel_data = zeros(dataset.N_ax*dataset.N_el, length(pw_indices));

    tic 
    % Complexity: (num_angle * number of apodization * num_of_smaple -
    % num_of subarray)
    for pw=1:length(pw_indices) %for every selected plane wave     
        %% Time-Of-Flight correction
        tof_data = zeros(dataset.N_ax*dataset.N_el,dataset.N_el);
        transmit_delay = Z*cos(dataset.angles(pw_indices(pw)))+X*sin(dataset.angles(pw_indices(pw)));
        
        for nrx=1:dataset.N_el
            receive_delay = sqrt((dataset.probe_geometry(nrx)-X).^2+(0-Z).^2);
            delay = (transmit_delay+receive_delay)/dataset.c0;
            tof_data(:,nrx) = interp1(time, dataset.data(:,nrx,pw_indices(pw)), reshape(delay,[],1),'spline',0);
        end 
        %% Apodization
        for i = 1:(2048*128)
            warning('off','MATLAB:singularMatrix');
            y = tof_data(i,:);
            N = length(y);
            L = max(round(N/3),1);            
            a = ones(L,1);
            if adaptive % set adaptive = 1
				% Code your adaptive capon estimator for the weights w here
                Rx_est = zeros(L,L);
				%calculate the R_x_est
                for j = 1:(N-L+1)
                    Rx_est = Rx_est + y(j:j+L-1)'*y(j:j+L-1);
                end
                % Add delta= L/100, to ensure well conditied cov matrix
                epsilon = (1/(L))*trace(Rx_est./(N-L+1));
                Rx_est = (Rx_est./(N-L+1))+epsilon.*eye(L);
                
%                 Rx_est = (Rx_est./(N-L+1));
                Rx_inv = (Rx_est)^(-1);
%                 Rx_inv = (Rx_est)^(-1);
                % Compute the capon filter
                w = Rx_inv*a*(a'*Rx_inv*a)^(-1);
                % Estimated amplitude
                z = 0;
                for j = 1:(N-L+1)
                   z = z + w'*y(j:j+L-1)';    
                end
                beamformed_channel_data(i, pw) = z./(N-L+1);

            else % set adaptive = 0, delay an sum beamformmer
                w = ones(N,1);
                beamformed_channel_data(i, pw) = w'*y';
                
            end               
        end
          
    end
    
    time = toc; 
    
    fprintf('Reconstruction time:  %f s \n',time)
    
    % Sum the contributions of different PW's. 
    % Alternatively, (adaptive) apodization could also be applied here. 
    beamformed_data = sum(beamformed_channel_data, 2); 
    beamformed_data(isnan(beamformed_data)) = 0;
    
    beamformed_data = reshape(beamformed_data, 2048, 128);
 
end
