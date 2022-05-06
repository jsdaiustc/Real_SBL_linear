clear;
close all;
% randn('seed',0);
% rand('seed',0);

Snap=200;                  % Number of snapshots
SNR = 0;                   % SNR 
M1=2;                       
M2=2;
M=M1+M2;                   % Number of element nested array
position=[0:M1 [2:M2]*(M1+1)-1];
resolution=3;              % grid interval
etc=M2*(M1);               % Maximum number of active grid points

%% generate the signal
True_DOAs=10*rand(1,2) +   [-30,10];
N_alpha=length(True_DOAs);
[X]=signal(M,position,True_DOAs,SNR, Snap);

%% the proposed method
[Pm_our,search_area_our]=real_VBI(X,Snap,resolution,position,etc);   
figure; stem(search_area_our,Pm_our)
