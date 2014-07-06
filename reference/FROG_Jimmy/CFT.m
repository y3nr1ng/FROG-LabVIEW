% Continuous Fourier Transform F(w)=CFT{f(t)} by F[k]=FFT{f[n]}
% input:  t = equally spaced (dt) sampling times; f = time domain function, length=N
% output: w = equally spaced (dw=2pi/(N*dt)) angular frequencies; F = spectrum
% units of t & w are related, e.g. t(ps)->w(THz)
% w=0 always exists, no matter whether t=0 is sampled
% F(p*dw)=fftshift{FFT(f[n])}*dt*exp[-j*p*dw*t0]: scaling & time-shift (see memo.)

% CAUTION: unlike actual CFT, FFT only takes care of a time window of f(t).
% FFT{f[n]} is equivalent to CFT{f(t)*rect(t)}=F(w) conv. sinc(w) in the best case.
% If the time window happens to be some specific caseS, FFT{f[n]} appears abnormal "nodes".
% For example, f(t)=exp(-(t-1).^2), t=-0.5~2.5, => FFT{f[n]} has nodes near the central peak.
% But t=0~2, or t=-1.5~3.5 won't show this bad phenomenon.
% This problem arises from finite time window, even one chooses trapz() to fulfill the CFT.

% v1. 2002/09/08; v2. 2006/04/29

function [w,F]=CFT(t,f)

N=length(t);			% # of sampled points
dt=t(2)-t(1);			% temporal resolution; temporal window = N*dt, not (N-1)*dt
dw=2*pi/N/dt;           % spectral resolution

F=fft(f);               % F[k]=DFT{f[n]}; f[n]=f(t=t0+n*dt), k,n=0,1,...,(N-1)
                        % fftshift([1 2, 3])=[3, 1 2]; w=0 at center
                        % fftshift([1 2, 3 4])=[3 4, 1 2]; w=0 is biased right

if mod(N,2)==0			% if N is even
   M=N/2;               % G[M]=F[0]~F(w=0), if G[q]=fftshift{F[k]}, q=0,1,...N-1
else					% if N is odd
   M=(N-1)/2;
end
w=[-M:(N-M-1)]*dw;      % sampled angular frequencies: p*dw
                        

F=fftshift(F)*dt.*exp(-j*w*t(1));