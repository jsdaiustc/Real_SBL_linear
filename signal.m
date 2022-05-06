function [X]=signal(M,position,True_DOAs,SNR,Snap)

N_theta=length(True_DOAs);
A=exp(-1i*pi*position'*sin(True_DOAs*pi/180));
Pow=sqrt((   (10).^(SNR/10)   )/2);
S=Pow*(randn(N_theta,Snap)+1i*randn(N_theta,Snap));
noise=sqrt(1/2)*(randn(M,Snap)+1i*randn(M,Snap));
X=A*S+noise;
