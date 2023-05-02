%% Assignment 1C Scenario 1: Maximized steered response DOA estimation
%
%   Implement DOA estimation
%   Requires tools from A2
%

clear all;
close all;
clear classes;
clc;
load Observations_1C
set(0,'defaultfigurecolor','w') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1Cs1 b and c)
% 1. Setup the ULA settings

J = 6;                  % Number of sensors
dy = 0;                 % meters of element spacing in y-direction
dx = 3.4e-2;                 % meters of element spacing in x-direction
nb_f = 2.5e3;               % narrowband (nb) frequency in Hz

% Setup an ULA array from the settings and plot the array configuration
my_array = arrays.ULA(J,dx,dy);

% Create a beamformer object and put settings in the beamformer object.
b = beamformer;
set(b, 'array',         my_array);
set(b, 'angles',        -180:0.1:180);
set(b, 'nb_frequency',  nb_f);
%%%%%%%%%%%%%%%% Write code here %%%%%%%%%%%%%
% calculate the steered response power here
R_x = zeros(size(observations,1),size(observations,1),size(observations,2));
for i = 1:size(observations,2)
    x_k = observations(:,i);
    R_x(:,:,i) = x_k*x_k';
end
R_x_avg = sum(R_x,3)/size(observations,2);
P_srp = zeros(size(b.angles));
% 1.c, Only use 4 microphone that closest to the origin
UseFourSensor_FLAG = false;
for i = 1:size(b.angles,2)
    a_theta = b.array_response_vector(b.angles(i),b.nb_frequency);
    if(UseFourSensor_FLAG == true)
        a_theta(1) = 0; a_theta(6)=0;
        J = 4;
    end
    P_srp(i)=abs(a_theta'*R_x_avg*a_theta./J);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);plot(b.angles,10*log10(P_srp),'LineWidth',1.5);xlabel('angle [degrees]');ylabel('steered response power [dB]');grid on;
if(UseFourSensor_FLAG == true)
    title(strcat('Steered response spatial power spectrum, 4 sensor ULA'));
else
    title(strcat('Steered response spatial power spectrum,',num2str(b.array.number_of_sensors), ' sensor ULA'));
end
axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare result in 1.b and 1.c
a_theta = b.array_response_vector(b.angles,b.nb_frequency);
P_srp_1c = zeros(size(b.angles));
UseFourSensor_FLAG = false;
if(UseFourSensor_FLAG == true)
    a_theta(1,:) = int8(0); a_theta(6,:)=int8(0);
    J = 4;
end
for i = 1:size(b.angles,2)
    P_srp_1c(i)=abs(a_theta(:,i)'*R_x_avg*a_theta(:,i)./J);
end
localmax_Psrp_1b = islocalmax(P_srp);
localmax_Psrp_1c = islocalmax(P_srp_1c);
figure; grid on;hold on;
plot(b.angles,10*log10(P_srp),'LineWidth',1.5);
plot(b.angles,10*log10(P_srp_1c),'LineWidth',1.5);
xline(b.angles(localmax_Psrp_1b),'b--')
xline(b.angles(localmax_Psrp_1c),'r--')
xlabel('angle [degrees]');ylabel('steered response power [dB]');grid on;
title('Compare steered spatial power spectrum between 4 and 6 sensors ULA')
legend('6 Sensor ULA','4 Sensors ULA')
axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1Cs2 a)
% return %%remove this to continue
% calculate the spatial spectrum using MUSIC here
% 1. Estimate the Rx
% Rank the eigen value, 
Num_of_source = 2;J = 6;
for i = 1:size(observations,2)
    x_k = observations(:,i);
    R_x(:,:,i) = x_k*x_k';
end
R_x_est = sum(R_x,3)/size(observations,2);
% 2. Compute the EVD of Rx
[Us,Lambda_s] = eigs(R_x_est,Num_of_source);
[Un,Lambda_n] = eigs(R_x_est,(J-Num_of_source),'sm');
P_n = Un*Un';
P_music = zeros(size(b.angles));
% Compute projection matrix
for i = 1:size(b.angles,2)
    a_theta = b.array_response_vector(b.angles(i),b.nb_frequency);
    P_music(i) =abs(J/(a_theta'*P_n*a_theta));
end
figure;plot(b.angles,10*log10(P_music),'LineWidth',1.5);xlabel('angle [degrees]');ylabel('Spatial spectrum [dB]');title(strcat('MUSIC spatial pseudo spectrum ',num2str(b.array.number_of_sensors), ' sensor ULA'));
grid on; axis tight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1Cs2 c)
% return %%remove this to continue
% using the number of sources wrongly and apply the spectral-MUSIC
% algorithm
P_music_3sources = spectrum_music(3,J,observations,b);
P_music_1sources = spectrum_music(1,J,observations,b);
P_music_4sources = spectrum_music(4,J,observations,b);
P_music_5sources = spectrum_music(5,J,observations,b);
P_music_2sources = spectrum_music(2,J,observations,b);
localmax_Psrp_2s = islocalmax(P_music_2sources);
localmax_Psrp_1s = islocalmax(P_music_1sources);
localmax_Psrp_3s = islocalmax(P_music_3sources);
localmax_Psrp_4s = islocalmax(P_music_4sources);
localmax_Psrp_5s = islocalmax(P_music_5sources);
figure; hold on;
plot(b.angles,10*log10(P_music_2sources),'b','LineWidth',1.5);
plot(b.angles,10*log10(P_music_3sources),'r','LineWidth',1.5);
plot(b.angles,10*log10(P_music_1sources),'k','LineWidth',1.5);
xline(b.angles(localmax_Psrp_2s),'b--','LineWidth',1.25)
xline(b.angles(localmax_Psrp_1s),'k--','LineWidth',1.25)
xline(b.angles(localmax_Psrp_3s),'r--','LineWidth',1.25)
xline(b.angles(localmax_Psrp_4s),'LineWidth',1.25)
xline(b.angles(localmax_Psrp_5s),'LineWidth',1.25)
plot(b.angles,10*log10(P_srp),'LineWidth',1.5)
(b.angles(localmax_Psrp_2s))
(b.angles(localmax_Psrp_1s))
(b.angles(localmax_Psrp_3s))
(b.angles(localmax_Psrp_4s))
(b.angles(localmax_Psrp_5s))
% plot(b.angles,10*log10(P_music_4sources),'LineWidth',1.5);
% plot(b.angles,10*log10(P_music_5sources),'LineWidth',1.5);
% plot(b.angles,(P_music),'LineWidth',1.5);
% plot(b.angles,(P_music_3sources),'LineWidth',1.5);
% plot(b.angles,(P_music_1sources),'LineWidth',1.5);
% plot(b.angles,(P_music_4sources),'LineWidth',1.5);
% plot(b.angles,(P_music_5sources),'LineWidth',1.5);
xlabel('angle [degrees]');ylabel('Spatial spectrum [dB]');
title(strcat('MUSIC spatial pseudo spectrum ',num2str(b.array.number_of_sensors), ' sensor ULA'));
grid on; axis tight;legend('2 sources(correct)','3 sources (wrong)','1 source (wrong)')
xlim([-90 90])
function P_music = spectrum_music(Num_of_source,J,observations,b)
    %1. Compute the covariance
    for i = 1:size(observations,2)
        x_k = observations(:,i);
        R_x(:,:,i) = x_k*x_k';
    end
    R_x_est = sum(R_x,3)/size(observations,2);
    % 2. Compute the EVD of Rx
    [Us,Lambda_s] = eigs(R_x_est,Num_of_source);
    [Un,Lambda_n] = eigs(R_x_est,(J-Num_of_source),'sm');
    P_n = Un*Un';
    P_music = zeros(size(b.angles));
    % Compute projection matrix
    for i = 1:size(b.angles,2)
        a_theta = b.array_response_vector(b.angles(i),b.nb_frequency);
        P_music(i) =abs(J/(a_theta'*P_n*a_theta));
    end
end
% Observation:
% More sensor: P_n smaller -> PSM larger
% Less sensor: P_n larger -> PSM smaller

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
