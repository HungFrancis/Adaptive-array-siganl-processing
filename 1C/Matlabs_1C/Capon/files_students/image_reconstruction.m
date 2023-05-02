clear all
close all

%% Load data
dataset_name = 'point_scat';
% dataset_name = 'cyst';

load(['rf_data/', dataset_name, '.mat']);

dataset.data = data;
dataset.c0 = 1540;
dataset.f0 = 5.208e6;
dataset.fs = dataset.f0*4;

sz = size(data);
dataset.N_ax = sz(1); % Number of axial samples per RF line
dataset.N_el = sz(2); % Number of receiving elements
dataset.N_aq = sz(3); % Number of plane wave acquisitions

dataset.angles = linspace(-16*pi/180,16*pi/180, dataset.N_aq); % Define 75 equally spaced transmit angles

% Probe geometry
pitch = 0.3;
dataset.probe_geometry = 1e-3*(0:pitch:127*pitch);
dataset.probe_geometry = dataset.probe_geometry-mean(dataset.probe_geometry);

% Define the scan region which determines the grid of focussing points. For
% simplicity this grid currently containts the same amount of axial and
% lateral samples as the raw RF data, however this is not required (because
% of interpolation during TOF correction).
scan_size = [-0.019, 0.019; 0 ,0.063]; % Size of imaging region (xmin,xmax;zmin,zmax)
scan.x = linspace(scan_size(1,1),scan_size(1,2), dataset.N_el)';
scan.z = linspace(scan_size(2,1),scan_size(2,2), dataset.N_ax)';

%% Beamforming
pw_indices = 38; % 1:75 to use all PW's 
% pw_indices = 1:75;
% % 
adaptive = 0;
beamformed_data = beamformer(dataset, scan, pw_indices, adaptive);

adaptive = 1;
beamformed_data_mv = beamformer(dataset, scan, pw_indices, adaptive);

%% Envelope detection and logarithmic compression
env_data = abs(hilbert(beamformed_data));
image = 20*log10(env_data./max(env_data(:)));

env_data_mv = abs(hilbert(beamformed_data_mv));
image_mv = 20*log10(env_data_mv./max(env_data_mv(:)));
%% Visualize data
figure();
subplot(1,2,1)
imagesc((scan.x)*1e3,(scan.z)*1e3,image); 
colormap gray; caxis([-60 0]); colorbar; 
axis equal tight;
xlabel('x [mm]');
ylabel('z [mm]');
title('Beamformed image');
subplot(1,2,2)
imagesc((scan.x)*1e3,(scan.z)*1e3,image_mv); 
colormap gray; caxis([-60 0]); colorbar; 
axis equal tight;
xlabel('x [mm]');
ylabel('z [mm]');
title('Adaptively beamformed image');

%%%%%%%%%%%%% Got the image from otehr angles %%%%%%%%%%%%%%%%%%%%%%%%%%
Compare_angle_FLAG = false;
if Compare_angle_FLAG

    pw_indices_middle = 38;
    pw_indices_left = 1;
    pw_indices_right = 75;
    %
    adaptive = 0;
    beamformed_data_left = beamformer(dataset, scan, pw_indices_left, adaptive);
    adaptive = 1;
    beamformed_data_mv_left = beamformer(dataset, scan, pw_indices_left, adaptive);
    % Use the previous reuslt to save time
%     adaptive = 0;
%     beamformed_data_m = beamformer(dataset, scan, pw_indices_middle, adaptive);
%     adaptive = 1;
%     beamformed_data_mv_m = beamformer(dataset, scan, pw_indices_middle, adaptive);
    %
    adaptive = 0;
    beamformed_data_right = beamformer(dataset, scan, pw_indices_right, adaptive);
    adaptive = 1;
    beamformed_data_mv_right = beamformer(dataset, scan, pw_indices_right, adaptive);
     %
     % Envelope detection and logarithmic compression
    env_data_left = abs(hilbert(beamformed_data_left));
    env_data_right = abs(hilbert(beamformed_data_right));
    image_left = 20*log10(env_data_left./max(env_data_left(:)));
    image_right = 20*log10(env_data_right./max(env_data_right(:)));
    env_data_mv_left = abs(hilbert(beamformed_data_mv_left));
    env_data_mv_right = abs(hilbert(beamformed_data_mv_right));
    image_mv_left = 20*log10(env_data_mv_left./max(env_data_mv_left(:)));
    image_mv_right = 20*log10(env_data_mv_right./max(env_data_mv_right(:)));
    %%
    % Visualize data
    figure();
    s=sprintf('Beeamformed Image from different angles \n (First row: Non-adaptive, Second row: Adaptive)');
    sgtitle(s)
    set(gca, 'LooseInset', get(gca,'TightInset'))
    % Non-adaptive
    subplot(2,3,1)
    imagesc((scan.x)*1e3,(scan.z)*1e3,image_left); 
    colormap gray; caxis([-60 0]); colorbar; 
    axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Non-adaptive: -16 deg.');
    subplot(2,3,2)
    imagesc((scan.x)*1e3,(scan.z)*1e3,image); 
    colormap gray; caxis([-60 0]); colorbar; 
    axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Non-adaptive: 0 deg.');
    subplot(2,3,3)
    imagesc((scan.x)*1e3,(scan.z)*1e3,image_right); 
    colormap gray; caxis([-60 0]); colorbar; 
    axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Non-adaptive: 16 deg.');
    % Adaptive
    subplot(2,3,4)
    imagesc((scan.x)*1e3,(scan.z)*1e3,image_mv_left); 
    colormap gray; caxis([-60 0]); colorbar; 
    axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Adaptive: -16 deg.')
    subplot(2,3,5)
    imagesc((scan.x)*1e3,(scan.z)*1e3,image_mv); 
    colormap gray; caxis([-60 0]); colorbar; 
    axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Adaptive: 0 deg.');
    subplot(2,3,6)
    imagesc((scan.x)*1e3,(scan.z)*1e3,image_mv_right); 
    colormap gray; caxis([-60 0]); colorbar; 
    axis equal tight;
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('Adaptive: 16 deg.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%