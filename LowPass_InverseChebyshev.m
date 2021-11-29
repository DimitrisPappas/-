clear all
clc

fp = 4400;
fs = 7480;
amax = 0.25;
amin = 28;

% Poles
wo1 = 54291.88089;
wo2 = 44916.20468;
wo3 = 36545.35063;
% Q
q1 = 0.5;
q2 = 0.74704;
q3 = 2.40383;
% Zeros
wz2 = 49416.754815;
wz3 = 79758.082064;
% gain
k1 = 1;         % gain monada 1
k2 = 0.2614;    % gain monada 2
k3 = 0.7245;    % gain monada 3
K = 0.916;      % gain gia na ftaso to epithumito kerdos (0 dB)

% Transfer Function
sys1 = tf([0,0,wo1],[0,1,wo1]);             % (wo1)
sys2 = tf([1,0,wz3^2],[1,wo2/q2,wo2^2]);    % (wo2,wz3)
sys3 = tf([1,0,wz2^2],[1,wo3/q3,wo3^2]);    % (wo3,wz2)

sys1 = k1*sys1;
sys2 = k2*sys2;
sys3 = k3*sys3;

sys = series(sys1,sys2);
sys = series(sys,sys3);

sys = K*sys;

a = inv(sys);   % to a thelo na exei monadiaio kerdos (0 dB)

% bode
% plot_transfer_function(sys1,wo1/(2*pi));
% plot_transfer_function(sys2,[10 wo2/(2*pi)]);
% plot_transfer_function(sys3,[10 wo3/(2*pi)]);
plot_transfer_function(sys,[1 fp fs]);
% plot_transfer_function(a,[1 fp fs]);
% ltiview('bodemag', sys);


% Square signal
Fs=100000;
T=1/Fs;
t= 0 : T : 0.02;
L=length(t);
f=Fs*(0:(L/2))/L;

% vin
ft=2e+3;
sqr = square(2*pi*ft*t, 40); % 40% positive
sqr = (sqr+1)/2;
vin = sqr;
vin = transpose(vin);
% dutycycle(vin, t);
%ipologiso fft vin
Vin = fft(vin);
Vin=abs(Vin/L);
Vin=Vin(1:round(L/2)+1);
Vin(2:end-1)=2*Vin(2:end-1);

% vout
vout=lsim(-sys,vin,t);
% dutycycle(vout, t);
Vout=fft(vout);
Vout=abs(Vout/L);
Vout=Vout(1:round(L/2)+1);
Vout(2:end-1)=2*Vout(2:end-1);

figure()
subplot(2,1,1);
plot(t,vin)
xlabel('sec')
ylabel('vin')
xlim([0.010 0.012])

subplot(2,1,2);
plot(t,vout)
xlabel('sec')
ylabel('vout')
xlim([0.010 0.014])


figure()
subplot(2,1,1);
plot(f,Vin)
xlabel('f(Hz)')
ylabel('|Vin|')

subplot(2,1,2);
plot(f,Vout)
xlabel('f(Hz)')
ylabel('|Vout|')


