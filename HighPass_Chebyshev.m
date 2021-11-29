clear all
clc

fp = 4000;
fs = 2222.2;
amax = 0.5278;
amin = 29;

% Poles
wo1 = 70459.0446;
wo2 = 36551.3979;
wo3 = 24741.8204;
% Q
q1 = 0.5;
q2 = 1.1913;
q3 = 4.6089;

% gain
k1 = 1;         % gain monada 1
k2 = 1;         % gain monada 2
k3 = 1;         % gain monada 3
K = 1.77828;    % gain gia na ftaso to epithumito kerdos (5 dB)

% Transfer Function
sys1 = tf([0,1,0],[0,1,wo1]);             % (wo1)
sys2 = tf([1,0,0],[1,wo2/q2,wo2^2]);    % (wo2)
sys3 = tf([1,0,0],[1,wo3/q3,wo3^2]);    % (wo3)

sys1 = k1*sys1;
sys2 = k2*sys2;
sys3 = k3*sys3;

sys = series(sys1,sys2);
sys = series(sys,sys3);

sys = K*sys;

a = inv(sys);

% bode
% plot_transfer_function(sys1,wo1/(2*pi));
% plot_transfer_function(sys2,wo2/(2*pi));
% plot_transfer_function(sys3,wo3/(2*pi));
plot_transfer_function(sys,[1e+6 fp fs]);
% plot_transfer_function(a,[1e+6 fp fs]);


% Signal
Fs=1000000;
T=1/Fs;
t=0:T:0.02;
L=length(t);
f=Fs*(0:(L/2))/L;

vin = cos(2*pi*0.4*fs*t) + 0.5*cos(2*pi*0.9*fs*t) + cos(2*pi*1.4*fp*t) + 0.7*cos(2*pi*2.4*fp*t) + 0.5*cos(2*pi*4.5*fp*t);


% fft vin
Vin = fft(vin);
Vin=abs(Vin/L);
Vin=Vin(1:round(L/2)+1);
Vin(2:end-1)=2*Vin(2:end-1);

% fft vout
vout=lsim(-sys,vin,t);
Vout=fft(vout);
Vout=abs(Vout/L);
Vout=Vout(1:round(L/2)+1);
Vout(2:end-1)=2*Vout(2:end-1);

figure()
subplot(2,1,1);
plot(t,vin)
xlabel('sec')
ylabel('vin')
xlim([0.010 0.014])

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




