clear all
clc
close all

aem = [8 3 9 1];

fo = 1000;
f1 = 650+25*aem(4);
f2 = fo^2/f1;
d = 2.1*(fo^2-f1^2)/f1;
f3 = (-d+sqrt(d^2+4*fo^2))/2;
f4 = fo^2/f3;

amin = 35-aem(3);
amax = 0.4+aem(4)/36;

wo = 2*pi*fo;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

bw = w2-w1;

Ws = (w4-w3)/(w2-w1);
Wp = 1;
Wp = Wp/Ws;
Ws = Ws/Ws;

n = acosh(((10^(amin/10)-1)/(10^(amax/10)-1))^(1/2))/acosh(1/Wp);
n = ceil(n);

e = (10^(amin/10)-1)^(-1/2);
a = asinh(1/e)/n;

% poles
% Low Pass Chebyshev
y = [22.5 67.5];        % degrees
s = sinh(a)*cosd(y);
w = cosh(a)*sind(y);
Wo = sqrt(s.^2+w.^2);
Q = Wo./(s.*2);

% Low Pass Inverse Chebyshev
Wo = Wo.^(-1);
% klimakopoioume gia na pame se Wp=1 Ws=2.1
Ws = Ws/Wp;
Wp = Wp/Wp;
Wo = Wo.*Ws;
S = Wo./Q./2;
W = sqrt(Wo.^2-S.^2);
Wo = sqrt(S.^2+W.^2);

% zeros
kappa = [1 3];
Wz = sec((kappa*pi)/(2*n));
% klimakopoioume gia na pame se Wp=1 Ws=2.1
Wz = Wz.*Ws;

% Band Pass Inverse Chebyshev Poles
qc = wo/bw;
% Geffe
C = S.^2 + W.^2;
D = S.*2/qc;
E = 4+C./qc^2;
G = sqrt(E.^2-4*D.^2);
Q = sqrt((E+G)/2)./D;
k_poles = S.*Q/qc;
Wgeffe = k_poles+sqrt(k_poles.^2-1);
% BP IC Poles
wo_poles = [0 0 0 0];
wo_poles(1) = wo/Wgeffe(1);
wo_poles(2) = wo*Wgeffe(1);
wo_poles(3) = wo/Wgeffe(2);
wo_poles(4) = wo*Wgeffe(2);

% Band Pass Inverse Chebyshev Zeros
k_zeros = 2+Wz.^2/qc^2;
x = (k_zeros+sqrt(k_zeros.^2-4))/2;
% BP IC Zeros
wz_zeros = [0 0 0 0];
wz_zeros(1) = wo*sqrt(x(1));
wz_zeros(2) = wo/sqrt(x(1));
wz_zeros(3) = wo*sqrt(x(2));
wz_zeros(4) = wo/sqrt(x(2));

% Q
q = [0 0 0 0];
q(1) = Q(1);
q(2) = Q(1);
q(3) = Q(2);
q(4) = Q(2);

% Transfer Function
sys1 = tf([1,0,wz_zeros(2)^2],[1,wo_poles(1)/q(1),wo_poles(1)^2]);    % (wo1,q1,wz2)
sys2 = tf([1,0,wz_zeros(1)^2],[1,wo_poles(2)/q(2),wo_poles(2)^2]);    % (wo2,q2,wz1)
sys3 = tf([1,0,wz_zeros(4)^2],[1,wo_poles(3)/q(3),wo_poles(3)^2]);    % (wo3,q3,wz4)
sys4 = tf([1,0,wz_zeros(3)^2],[1,wo_poles(4)/q(4),wo_poles(4)^2]);    % (wo4,q4,wz3)

% MONADES
% Monada 1
% Boctor      
BoctorHighPass( wz_zeros(2), wo_poles(1), q(1), 10^3, 10^(-8) );        % (wo1,q1,wz2)
k_m1 = 2;       % H = 2

% Monada 3
% BoctorHighPass( wz_zeros(4), wo_poles(3), q(3), 10^2, 10^(-8) );        % (wo3,q3,wz4)
% High Pass Notch 7.21
wo_m3 = 1;
wz_m3 = wz_zeros(4)/wo_poles(3);
k1_m3 = wo_m3^2/wz_m3^2 - 1;
k2_m3 = ((2+k1_m3)*q(3)^2)/((2+k1_m3)*q(3)^2+1);
k_m3 = k2_m3*(wo_m3^2/wz_m3^2);
R1_m3 = 1;
R2_m3 = q(3)^2*(k1_m3+2)^2;
R3_m3 = 1;
R4_m3 = q(3)^2*(k1_m3+2);
C_m3 = 1/(q(3)*(2+k1_m3));
C1_m3 = k1_m3*C_m3;
% klimakopoihsh
kf_m3 = wo_poles(3);
Creal_m3 = 10^(-8);
km_m3 = C_m3/(Creal_m3*kf_m3);
% pragmatikes times
R31 = km_m3*R1_m3;
R32 = km_m3*R2_m3;
R33 = km_m3*R3_m3;
R34 = km_m3*R4_m3;
C31 = Creal_m3;
C32 = C31;
C33 = C1_m3/(km_m3*kf_m3);
k_m3;

% Monada 2
% Boctor LPN 7.24(a)
wo_m2 = 1;
wz_m2 = wz_zeros(1)/wo_poles(2);
k1_m2 = ((wo_m2^2/wz_m2^2)+1)/2;
R1_m2 = 2/(k1_m2*wz_m2^2-1);
R2_m2 = 1/(1-k1_m2);
R3_m2 = (k1_m2/q(2)^2+k1_m2*wz_m2^2-1)/2;
R4_m2 = 1/k1_m2;
R5_m2 = 1;
R6_m2 = 1;
C1_m2 = k1_m2/(2*q(2));
C2_m2 = 2*q(2);
k_m2 = 1/((k1_m2/q(2)^2+k1_m2*wz_m2^2+1)/2);
% klimakopoihsh
kf_m2 = wo_poles(2);
C1real_m2 = 10^(-8);
km_m2 = C1_m2/(C1real_m2*kf_m2);
% pragmatikes times
R21 = R1_m2*km_m2;
R22 = R2_m2*km_m2;
R23 = R3_m2*km_m2;
R24 = R4_m2*km_m2;
R25 = R5_m2*km_m2;
R26 = R6_m2*km_m2;
C21 = C1_m2/(km_m2*kf_m2);
C22 = C2_m2/(km_m2*kf_m2);
k_m2;

% Monada 4
% Boctor LPN 7.24(a)
wo_m4 = 1;
wz_m4 = wz_zeros(3)/wo_poles(4);
k1_m4 = ((wo_m4^2/wz_m4^2)+1)/2;
R1_m4 = 2/(k1_m4*wz_m4^2-1);
R2_m4 = 1/(1-k1_m4);
R3_m4 = (k1_m4/q(4)^2+k1_m4*wz_m4^2-1)/2;
R4_m4 = 1/k1_m4;
R5_m4 = 1;
R6_m4 = 1;
C1_m4 = k1_m4/(2*q(4));
C2_m4 = 2*q(4);
k_m4 = 1/((k1_m4/q(4)^2+k1_m4*wz_m4^2+1)/2);
% klimakopoihsh
kf_m4 = wo_poles(4);
C1real_m4 = 10^(-8);
km_m4 = C1_m4/(C1real_m4*kf_m4);
% pragmatikes times
R41 = R1_m4*km_m4;
R42 = R2_m4*km_m4;
R43 = R3_m4*km_m4;
R44 = R4_m4*km_m4;
R45 = R5_m4*km_m4;
R46 = R6_m4*km_m4;
C41 = C1_m4/(km_m4*kf_m4);
C42 = C2_m4/(km_m4*kf_m4);
k_m4;

% Gain
% gain @jwo
ss = 1i*2*pi*fo;
k1_jwo = abs(k_m1*(ss^2+wz_zeros(2)^2)/(ss^2+ss*wo_poles(1)/q(1)+wo_poles(1)^2));
k2_jwo = abs(k_m2*(ss^2+wz_zeros(1)^2)/(ss^2+ss*wo_poles(2)/q(2)+wo_poles(2)^2));
k3_jwo = abs(k_m3*(ss^2+wz_zeros(4)^2)/(ss^2+ss*wo_poles(3)/q(3)+wo_poles(3)^2));
k4_jwo = abs(k_m4*(ss^2+wz_zeros(3)^2)/(ss^2+ss*wo_poles(4)/q(4)+wo_poles(4)^2));
A = k1_jwo*k2_jwo*k3_jwo*k4_jwo;
% A = k_m1*k_m2*k_m3*k_m4
% gia na paroume kerdos 0 dB prosthetoume kerdos K sto telos
K = 1/A;

sys1 = sys1*k_m1;
sys2 = sys2*k_m2;
sys3 = sys3*k_m3;
sys4 = sys4*k_m4;

sys = series(sys1,sys2);
sys = series(sys,sys3);
sys = series(sys,sys4);

sys = sys*K;
a = inv(sys);

% Plots
% plot_transfer_function(sys1,wo_poles(1)/(2*pi))
% plot_transfer_function(sys2,wo_poles(2)/(2*pi))
% plot_transfer_function(sys3,wo_poles(3)/(2*pi))
% plot_transfer_function(sys4,wo_poles(4)/(2*pi))
plot_transfer_function(sys,[fo f1 f2 f3 f4])
% plot_transfer_function(a,[fo f1 f2 f3 f4])



% Signal
Fs=1000000;
T=1/Fs;
t=0:T:0.014;
L=length(t);
f=Fs*(0:(L/2))/L;


vin = cos(2*pi*(fo-(fo-f1)/2)*t) + 0.8*cos(2*pi*(fo+(fo+f1)/3)*t) + 0.8*cos(2*pi*0.4*f3*t) + 0.6*cos(2*pi*2.5*f4*t) + 0.5*cos(2*pi*3*f4*t);


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
subplot(2,1,1)
plot(t,vin)
xlabel('sec')
ylabel('vin')
xlim([0.010 0.014])

subplot(2,1,2)
plot(t,vout)
xlabel('sec')
ylabel('vout')
xlim([0.010 0.014])


figure()
subplot(2,1,1);
plot(f,Vin(1:end-1))
xlabel('f(Hz)')
ylabel('|Vin|')
xlim([0 50000])

subplot(2,1,2);
plot(f,Vout(1:end-1))
xlabel('f(Hz)')
ylabel('|Vout|')
xlim([0 50000])











