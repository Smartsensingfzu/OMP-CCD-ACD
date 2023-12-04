clc;
clear;
N=1024;
M=256;
fs=1024;%采样频率
ts=1/fs;%采样间隔
Ts=1:N; %采样序列
f1=4;f2=8;f3=8;f4=12;
f=0.3*sin(2*pi*f1*Ts*ts)+0.6*sin(2*pi*f2*Ts*ts)+0.9*sin(2*pi*f3*Ts*ts)+0.1*sin(2*pi*f4*Ts*ts);
K=16;
[f,f4]=mapminmax(f,-1,1);
fft_f=fft(f)/sqrt(N);
load ('PN.mat');
load('guance.mat');
load('w1.mat');
Phi=guance*diag(PN);
y=Phi*f.';
[hat_x,hat_y,Aug_t,pos_array,Psi_t,aug_y,T,Psi,wn] = omp2(w1,y,K,M,N,Phi);
[pos_array_omp1,aug_y_omp1,hat_x_omp1,hat_y_omp1] = omp1(y,K,M,N,Phi);
figure(1)
plot(f,'k.-');
hold on;
% plot(hat_x,'r.-');
% hold on;
plot(hat_x_omp1,'g.--');
hold on;
axis([0 100 -1 1]);
xlabel('时间'),ylabel('Signal Amplitude ') 
% legend('原始信号','近似OMP','原始OMP')  
grid on;
