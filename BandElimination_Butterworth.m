clear all
clc
close all

aem = [8 3 9 1];

fo = 1750;
f1 = 1000+50*aem(4);
f2 = fo^2/f1;
d = (fo^2-f1^2)/f1/3.5;
f3 = (-d+sqrt(d^2+4*fo^2))/2;
f4 = fo^2/f3;

amin = 32+aem(3)*5/9;
amax = 0.4+aem(4)*0.25/9;

wo = 2*pi*fo;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

bw = w2-w1;

Ws = (w2-w1)/(w4-w3);
Wp = 1;

n = log((10^(amin/10)-1)/(10^(amax/10)-1))/(2*log(Ws/Wp));
n = ceil(n);

% poloi Low Pass Buttherworth @n=5 apo pinaka 9.1
sk1_LP = 1;
sk2_LP = 0.8090;
wk2_LP = 0.5877;
sk3_LP = 0.3090;
wk3_LP = 0.9510;

% Hmiseias Isxuos @3dB
% Low Pass 
% Wo = Wp/((10^(amax/10)-1)^(1/(2*n)));
Wo = Ws/((10^(amin/10)-1)^(1/(2*n)));
% High Pass
Wo_HP = 1/Wo;

% poloi High Pass Butterworth
sk1_HP = Wo_HP*sk1_LP;
sk2_HP = Wo_HP*sk2_LP;
wk2_HP = Wo_HP*wk2_LP;
sk3_HP = Wo_HP*sk3_LP;
wk3_HP = Wo_HP*wk3_LP;

% M/S pragmatikou polou
qc = wo/bw;
S1 = sk1_HP;
Q1 = qc/S1;
y1 = acosd(1/(2*Q1));
wo1 = wo;

% M/S migadikou polou
% 1o zeygari
S2 = sk2_HP;
W2 = wk2_HP;
C = S2^2+W2^2;
D = 2*S2/qc;
E = 4+C/qc^2;
G = sqrt(E^2-4*D^2);
Q2 = sqrt((E+G)/2)/D;
kappa = S2*Q2/qc;
Wgeffe = kappa+sqrt(kappa^2-1);
wo2 = wo/Wgeffe;
wo3 = wo*Wgeffe;
% 2o zeygari
S3 = sk3_HP;
W3 = wk3_HP;
C = S3^2+W3^2;
D = 2*S3/qc;
E = 4+C/qc^2;
G = sqrt(E^2-4*D^2);
Q3 = sqrt((E+G)/2)/D;
kappa = S3*Q3/qc;
Wgeffe = kappa+sqrt(kappa^2-1);
wo4 = wo/Wgeffe;
wo5 = wo*Wgeffe;

% poloi
wo_poles = [wo1 wo2 wo3 wo4 wo5];
% zeros
wz = wo;
% Q
Q = [Q1 Q2 Q2 Q3 Q3];

% Transfer Function
sys1 = tf([1,0,wz^2],[1,wo1/Q(1),wo1^2]);
sys2 = tf([1,0,wz^2],[1,wo2/Q(2),wo2^2]);
sys3 = tf([1,0,wz^2],[1,wo3/Q(3),wo3^2]);
sys4 = tf([1,0,wz^2],[1,wo4/Q(4),wo4^2]);
sys5 = tf([1,0,wz^2],[1,wo5/Q(5),wo5^2]);

% MONADES
% monada 1 - Notch 7.21
wo_m1 = 1;
wz_m1 = wz/wo1;
k1_m1 = wo_m1^2/wz_m1^2-1;
k2_m1 = ((2+k1_m1)*Q(1)^2)/((2+k1_m1)*Q(1)^2+1);
k_m1 = k2_m1*(wo_m1^2/wz_m1^2);
R1_m1 = 1;
R2_m1 = Q(1)^2*(k1_m1+2)^2;
R3_m1 = 1;
R4_m1 = Q(1)^2*(k1_m1+2);
C_m1 = 1/(Q(1)*(k1_m1+2));
% klimakopoihsh
kf_m1 = wo1;
Creal_m1 = 10^(-7);
km_m1 = C_m1/(Creal_m1*kf_m1);
% pragmatikes times
R11 = R1_m1*km_m1;
R12 = R2_m1*km_m1;
R13 = R3_m1*km_m1;
R14 = R4_m1*km_m1;
C11 = Creal_m1;
k_m1;

% monada 2 - LPN 7.23
wo_m2 = 1;
wz_m2 = wz/wo2;
C_m2 = 1/(2*Q(2));
R1_m2 = 1;
R2_m2 = 4*Q(2)^2;
R3_m2 = wz_m2^2/(2*Q(2)^2);
R4_m2 = 1;
R5_m2 = (4*Q(2)^2)/(wz_m2^2-1);
k_m2 = 1/(R3_m2+1);
% klimakopoihsh
kf_m2 = wo2;
Creal_m2 = 10^(-7);
km_m2 = C_m2/(Creal_m2*kf_m2);
% pragmatikes times
R21 = R1_m2*km_m2;
R22 = R2_m2*km_m2;
R23 = R3_m2*km_m2;
R24 = R4_m2*km_m2;
R25 = R5_m2*km_m2;
C21 = Creal_m2;
C22 = Creal_m2;
k_m2;

% monada 3 - HPM 7.21
wo_m3 = 1;
wz_m3 = wz/wo3;
k1_m3 = wo_m3^2/wz_m3^2-1;
k2_m3 = ((2+k1_m3)*Q(3)^2)/((2+k1_m3)*Q(3)^2+1);
k_m3 = k2_m3*(wo_m3^2/wz_m3^2);
R1_m3 = 1;
R2_m3 = Q(3)^2*(k1_m3+2)^2;
R3_m3 = 1;
R4_m3 = Q(3)^2*(k1_m3+2);
C_m3 = 1/(Q(3)*(k1_m3+2));
C1_m3 = k1_m3*C_m3;
% klimakopoihsh
kf_m3 = wo3;
Creal_m3 = 10^(-7);
km_m3 = C_m3/(Creal_m3*kf_m3);
% pragmatikes times
R31 = R1_m3*km_m3;
R32 = R2_m3*km_m3;
R33 = R3_m3*km_m3;
R34 = R4_m3*km_m3;
C31 = Creal_m3;
C32 = Creal_m3;
C33 = C1_m3/(km_m3*kf_m3);
k_m3;

% monada 4 - LPN 7.23
wo_m4 = 1;
wz_m4 = wz/wo4;
C_m4 = 1/(2*Q(4));
R1_m4 = 1;
R2_m4 = 4*Q(4)^2;
R3_m4 = wz_m4^2/(2*Q(4)^2);
R4_m4 = 1;
R5_m4 = (4*Q(4)^2)/(wz_m4^2-1);
k_m4 = 1/(R3_m4+1);
% klimakopoihsh
kf_m4 = wo4;
Creal_m4 = 10^(-7);
km_m4 = C_m4/(Creal_m4*kf_m4);
% pragmatikes times
R41 = R1_m4*km_m4;
R42 = R2_m4*km_m4;
R43 = R3_m4*km_m4;
R44 = R4_m4*km_m4;
R45 = R5_m4*km_m4;
C41 = Creal_m4;
C42 = Creal_m4;
k_m4;

% monada 5 - HPN 7.21
wo_m5 = 1;
wz_m5 = wz/wo5;
k1_m5 = wo_m5^2/wz_m5^2-1;
k2_m5 = ((2+k1_m5)*Q(5)^2)/((2+k1_m5)*Q(5)^2+1);
k_m5 = k2_m5*(wo_m5^2/wz_m5^2);
R1_m5 = 1;
R2_m5 = Q(5)^2*(k1_m5+2)^2;
R3_m5 = 1;
R4_m5 = Q(5)^2*(k1_m5+2);
C_m5 = 1/(Q(5)*(k1_m5+2));
C1_m5 = k1_m5*C_m5;
% klimakopoihsh
kf_m5 = wo5;
Creal_m5 = 10^(-7);
km_m5 = C_m5/(Creal_m5*kf_m5);
% pragmatikes times
R51 = R1_m5*km_m5;
R52 = R2_m5*km_m5;
R53 = R3_m5*km_m5;
R54 = R4_m5*km_m5;
C51 = Creal_m5;
C53 = C1_m5/(km_m5*kf_m5);
k_m5;

% Gain
kd = 1;
A = k_m1*k_m2*k_m3*k_m4*k_m5*(wz/wo1)^2*(wz/wo2)^2*(wz/wo3)^2*(wz/wo4)^2*(wz/wo5)^2;
K = kd/A;

sys1 = sys1*k_m1;
sys2 = sys2*k_m2;
sys3 = sys3*k_m3;
sys4 = sys4*k_m4;
sys5 = sys5*k_m5;

% Total Transfer Function
sys = series(sys1,sys2);
sys = series(sys,sys3);
sys = series(sys,sys4);
sys = series(sys,sys5);

sys = sys*K;
a = inv(sys);

% Plots
% plot_transfer_function(sys1,wo1/(2*pi))
% plot_transfer_function(sys2,wo2/(2*pi))
% plot_transfer_function(sys3,wo3/(2*pi))
% plot_transfer_function(sys4,wo4/(2*pi))
% plot_transfer_function(sys5,wo5/(2*pi))
plot_transfer_function(sys,[ 1 fo f1 f2 f3 f4 50000])
% plot_transfer_function(a,[fo f1 f2 f3 f4])


% Signal
Fs=1000000;
T=1/Fs;
t=0:T:0.014;
L=length(t);
f=Fs*(0:(L/2))/L;


vin = 0.5*cos(2*pi*(fo-(fo-f3)/2)*t) + 0.8*cos(2*pi*(fo+(fo+f3)/3)*t) + 0.8*cos(2*pi*0.4*f1*t) + 0.6*cos(2*pi*2.5*f2*t) + 1.2*cos(2*pi*3*f2*t);


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



