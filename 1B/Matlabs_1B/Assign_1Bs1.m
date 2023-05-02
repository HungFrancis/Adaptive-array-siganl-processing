%% Assignment 1B Scenario 1: Narrowband beamformers
%
%   Implement tools for narrowband beamformer design and evaluation
%   Result plots: sensor locations, beampatterns in plot format.
%   Note:   Each step uses different functions of the beamformermer class.
%           These functions are denoted with each step.
%

clear all;
close all;
clear classes;
clc;
set(0,'defaultfigurecolor','w') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1BS1 c)
% 1. Setup the ULA settings as in figure 1 of the assignment
J = 4;                  % Number of sensors
dy = 0;                 % meters of element spacing in y-direction
dx = 0.034;                 % meters of element spacing in x-direction
nb_f = 2500;               % narrowband (nb) frequency in Hz

% Setup an ULA array from the settings and plot the array configuration
my_array = arrays.ULA(J,dx,dy);
%%%%% For 1.d
% rotated_mat = [cosd(90) sind(90);-sind(90) cosd(90)];
% my_array.sensor_positions = my_array.sensor_positions * rotated_mat;
%%%%%%
% Use the plot function that belongs to the array class (in @array
% directory)
figure;
my_array.plot();

% Create a beamformer object and put settings in the beamformer object.
b = beamformer;
set(b, 'array',         my_array);
set(b, 'angles',        -180:1:180);
set(b, 'nb_frequency',  nb_f);

% Display all properties of the beamformer b:
b

% 2. Implement the array_response_vector.m method that is located in the
%    @beamformer folder.
% 3. Implement the calc_nb_beampattern.m method that is located in the
%    @beamformer folder
% 4. Verify the result of the matlab function with your answer at a)

% Set the beamformer weights to 1
b.nb_weights = ones(J,1);
b.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
subplot(1,2,1);b.plot_nb;grid on;
subplot(1,2,2); %plot your equation from a) here
v = wave_vector(30);
c = 340; fd = 2.5e3;
p = b.array.sensor_positions;
w = ones(1,J);
angle = linspace(-180,180);
array_response = @(theta) abs((sin(pi.*sind(theta)))./(sin((pi/4).*sind(theta))));
b_theta = (array_response(angle).^2)/(J^2);
linspec = {'b-','LineWidth',2};
plot(angle,10*log10(b_theta),linspec{:});
xlabel('Angle in [degrees]')
ylabel('Beamformer gain in [dB]')
axis tight
title('Beampattern for narrowband beamformer (reference)')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1BS1 d)
% 1. Implement the beam_steering_nb.m method that is located in the
%    @beamformer folder.
% 2. Verify the result of the matlab function by visual inspection.

% Remove this return to continue with the assignment
% return;

theta = [30]; % row vector containing angles for which constraints hold
target = [1]; % row vector containing target values for the beampattern
b.beam_steering(theta, target);
b.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b.plot_nb([],theta, {'k-.','LineWidth',2});
title("Filter weights w that give the array unit response at θ_d = 30 deg.")
ylim([-60 0]);
grid on ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1BS1 f)
% 1. Add the undesired source direction and make sure that the beamformer
% has unity response at 30 degrees and a zero response at -60 degrees.

% Remove this return to continue with the assignment
% return;

theta = [30 29]; % row vector containing angles for which constraints hold
target = [1;0]; % row vector containing target values for the beampattern
b.beam_steering(theta, target);
b.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b.plot_nb([],theta, {'k-.','LineWidth',2});
% ylim([-45 20]);
grid on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1BS1 g)
% 1. Implement the rectangular_array.m function that is located in the
%    + array folder.
J = 4;                  % Number of sensors
dy = 1.7e-2;                 % meters of element spacing in y-direction
dx = 3.4e-2;                 % meters of element spacing in x-direction
nb_f = 2500;               % narrowband (nb) frequency in Hz
my_array = arrays.rectangular_array(J,dx,dy);
b2 = beamformer;
set(b2, 'array',         my_array);
set(b2, 'angles',        -180:1:180);
set(b2, 'nb_frequency',  nb_f);
figure;
my_array.plot();

theta = [30 -60]; % row vector containing angles for which constraints hold
target = [1;0]; % row vector containing target values for the beampattern
b2.beam_steering(theta, target);
b2.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b2.plot_nb([],theta, {'k-.','LineWidth',2});
ylim([-50 -10]);
grid on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign 1BS1 h)
% 1. Compare the square array results for different frequencies
J = 9;                  % Number of sensors
dy = 1.7e-2;                 % meters of element spacing in y-direction
dx = 3.4e-2;                 % meters of element spacing in x-direction
nb_f = 4000;               % narrowband (nb) frequency in Hz
my_array = arrays.rectangular_array(9,dx,dy);
b2 = beamformer;
set(b2, 'array',         my_array);
set(b2, 'angles',        -180:1:180);
set(b2, 'nb_frequency',  nb_f);
figure;
my_array.plot();

theta = [30 -60]; % row vector containing angles for which constraints hold
target = [1;0]; % row vector containing target values for the beampattern
b2.beam_steering(theta, target);
b2.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b2.plot_nb([],theta, {'k-.','LineWidth',2});
ylim([-60 -10]);
grid on 












