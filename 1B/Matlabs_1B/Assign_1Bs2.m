clear all;clc;close all
load Assign_1B_testsignals.mat
set(0,'defaultfigurecolor','w') 
DOA=30;
fs=25e3;

J = 4;                  % Number of sensors
dy = 0;                 % meters of element spacing in y-direction
dx = 3.4e-2;                 % meters of element spacing in x-direction
ULA_array = arrays.ULA(J,dx,dy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1B scenario 2 d)
%  1.Calculate the required delay for each sensor signal that will
%    time-align the input signals
% d={0,1, ..., dmax}
% fd={0,1, ... L-1}
% L={1,2,3,...}
% dmax={0,1,2,...}
% T = linspace(0,50,1000);
NSample = 20;
% t = linspace(0,50,NSample+1);
t = 0:0.01:NSample-1;
d0 = sinc((3-t));
d1 = sinc((5.15-t));
figure
% stem(T,d0(T),'b-')
plot(t,d0,'b-','LineWidth',2)
title("d_0(t)"); ylabel("Amplitude");xlabel("sample");grid on
figure
% stem(T,d1(T),'b-')
plot(t,d1,'b-','LineWidth',2)
grid on;ylabel("Amplitude");xlabel("sample");title("d_1(t)")
%%
%  2.Use delay.m to generate impulse responses that will perform the required delay. Take
tau0 = 3;
tau1 = 5.15;
FracDelay = tau0; % Set the number of delay sample
FracDelayInt = floor(FracDelay); % Largest integer that is less than the fractional delay
FracDelayRem = FracDelay - FracDelayInt; %Fractional part
L = 100; 
h = delay(FracDelayInt, round(FracDelayRem * L), L, ceil(FracDelay)); %design the filter
ProcessingDelay = length(h) - ceil(FracDelay); %The processing delay is determined by the length of the filter and dmax
% d0_hat = h(floor(ProcessingDelay / 2) + 1: floor(ProcessingDelay / 2) + NSample).*[zeros(1, ceil(FracDelay)-1), ones(1, NSample - ceil(FracDelay)+1)]; 
d0_hat = h(floor(ProcessingDelay / 2) + 1: floor(ProcessingDelay / 2) + NSample); 
%
FracDelay = tau1; % delay 3.75 samples, when the FracDelay is integer, it would be a delta funtion 
FracDelayInt = floor(FracDelay); % Largest integer that is less than the fractional delay
FracDelayRem = FracDelay - FracDelayInt; %Fractional part
L = 100; 
h = delay(FracDelayInt, round(FracDelayRem * L), L, ceil(FracDelay)); %design the filter
ProcessingDelay = length(h) - ceil(FracDelay); %The processing delay is determined by the length of the filter and dmax
d1_hat = h(floor(ProcessingDelay / 2) + 1: floor(ProcessingDelay / 2) + NSample); 
figure;hold on;box on;grid on;xlabel('sample');ylabel('Amplitude');
% Plot 
plot(t,d0,'LineWidth',1.25)
stem(0:NSample-1,d0_hat,'LineWidth',1.25)
legend('$$d_0(t)$$','$\hat{d}_0(t)$','Interpreter','Latex')
figure;hold on;box on;grid on;
plot(t,d1,'LineWidth',1.25)
stem(0:NSample-1,d1_hat,'LineWidth',1.25)
xlabel('sample');ylabel('Amplitude')
legend('$$d_1(t)$$','$\hat{d}_1(t)$','Interpreter','Latex')
%%
%  3.Filter the input signals with their corresponding filters. Verify
%    whether the signals have been time-aligned.
% Use the delay filter on the signal
% X show the 4 sensors te 1000 samples
% when the FracDelay is integer, it would be a delta funtion 
FilteredSignal = zeros(J,length(x));
for i = 1:J
    % The sensor that first receive signal -> delay to align with the farest sensor
    % The farest sensor-> no delay
    FracDelay = (dx*(J-i+1)*sind(DOA)/340)*fs;
    FracDelayInt = floor(FracDelay);
    FracDelayRem = FracDelay - FracDelayInt; %Fractional part
    DelayFilter = delay(FracDelayInt, round(FracDelayRem * 100),100, ceil(FracDelay));
     % Apply the delay filter to the corresponding sensor
    SignalFracDelay = conv(x(i,:),DelayFilter);
    ProcessingDelay = length(DelayFilter) - ceil(FracDelay);
    ResultSignal = SignalFracDelay(floor(ProcessingDelay / 2) + 1: floor(ProcessingDelay / 2) + length(x));
    FilteredSignal(i,:) =ResultSignal;
%     % To cope with casuality
%     CasualFactor =  ceil(FracDelay);
%     FilteredSignal(i,:) = ResultSignal.*...
%     [zeros(1,CasualFactor), ones(1, length(x)-CasualFactor)]; 
end
figure;
plot(x(1,:))
hold on;grid on;box on;
plot(x(2,:))
plot(x(3,:))
plot(x(4,:))
title("Before time aligning")
legend("Sensor 1","Sensor 2","Sensor 3","Sensor 4")
xlabel('sample');ylabel('Amplitude')
figure;
hold on;grid on;box on;
plot(FilteredSignal(1,:))
legend
plot(FilteredSignal(2,:))
plot(FilteredSignal(3,:))
plot(FilteredSignal(4,:))
title("After time aligning")
legend("Sensor 1","Sensor 2","Sensor 3","Sensor 4")
xlabel('sample');ylabel('Amplitude')
% cope with casuality