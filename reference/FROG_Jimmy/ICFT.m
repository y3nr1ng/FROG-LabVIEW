% Continuous Inverse Fourier Transform f(t)=ICFT{F(w)} by f[n]=IFFT{F[k]}
% input:  w = equally spaced (dw) angular frequencies; F = spectrum, length=N
% output: t = equally spaced (dt=2pi/(N*dw)) sampling times; f = time domain function
% units of t & w are related, e.g. w(THz)->t(ps)
% t=0 always exists, no matter whether w=0 is sampled
% f(p*dt)=fftshift{IFFT(F[k])}*(N*dw/2pi)*exp[+j*w0*p*dt]: scaling & freq-shift (see memo.)
% v1. 2002/9/24; v2. 2006/04/29

function [t,f]=ICFT(w,F)

N=length(w);			% # of sampled points
dw=w(2)-w(1);			% spectral resolution; spectral window = N*dw, not (N-1)*dw
dt=2*pi/N/dw;           % temporal resolution

f=ifft(F);              % f[n]=IDFT{F[k]}; F[k]=F(w=w0+k*dw), n,k=0,1,...,(N-1)
                        % fftshift([1 2, 3])=[3, 1 2]; t=0 at center
                        % fftshift([1 2, 3 4])=[3 4, 1 2]; t=0 is biased right

if mod(N,2)==0			% if N is even
   M=N/2;               % g[M]=f[0]~f(t=0), if g[q]=fftshift{f[n]}, q=0,1,...N-1
else					% if N is odd
   M=(N-1)/2;
end
t=[-M:(N-M-1)]*dt;      % sampled angular frequencies: p*dt


f=fftshift(f)/dt.*exp(j*w(1)*t);