clear all;
close all;
clc;

load PGFROGdata;

dtau = tau(2) - tau(1);
w = 2.*pi.*f;

figure(1); % Plot the experimental data given by SDY
contourf(tau,f,I_PGFROG,'linestyle','none');
axis([-4 4 -2 2]);
title('Experimental Data of PG FROG Trace','fontsize',20);
xlabel('Delay \tau ( ps )','fontsize',20);
ylabel('Angular Frequency ( THz )','fontsize',20);

iteration_time = 200; % set the iteration time

A_noise = rand(1,length(w));
phase_noise = rand(1,length(w));
A_w = A_noise.*exp(-j.*phase_noise.*2.*pi); % give a random noise amplitude and random phase

[t,a_t]=ICFT(w, A_w); % initial guess a(t)

A_FROG = zeros(256,129);
I_FROG = zeros(256,129);
A_sig_new = zeros(256,129);
a_sig_new = zeros(256,129);
delay = (length(tau)-1)/2; % split the delay into half-left-half-right

for a = 1:iteration_time
    
    a
    
    for b = -delay : 1 : delay % Move the Delay Stage
        
        gate = circshift(abs(a_t).^2,b); % "Gate Constraint", gate(t,tau) = I(t-tau)
        a_sig = a_t.*gate;
        A_sig = fftshift(fft(a_sig))*(t(2)-t(1)); % Shift the pulse to the center, then do the Fourier Transform

        A_FROG(:,b+delay+1) = A_sig; % We'll get A_sig(w,tau)
        I_FROG(:,b+delay+1) = abs(A_sig).^2;
        
    end
    
    A_FROG = A_FROG./max(max(abs(A_FROG))); % Normalize max of A_FROG to 1
    I_FROG = I_FROG./max(max(abs(I_FROG))); % Normalize max of A_FROG to 1
    
    A_sig_new = A_FROG.*I_PGFROG.^0.5./abs(A_FROG); % Replace the magnitude with square root of I_FROG ( Data Constraint )
                                                    % We'll get A'_sig(w,tau)
    
    for c = 1:length(tau)
        
        a_sig_new(:,c) = ifft(ifftshift(A_sig_new(:,c)))/(t(2)-t(1)); % Shift the spectrum to the center then do the inverse Fourier Transform
        
    end
    
    for d = 1:length(t) % integrate a'_sig(t,tau) with respect to tau
        
        a_t(d) = sum(a_sig_new(d,:)).*dtau; % We'll get the next a_t
        
    end
            
    FROG_error(a) = (sum(sum((I_FROG - I_PGFROG))).^2./(length(f).*length(tau))).^0.5; % FROG error
    number(a) = a; % Iteration number

end

figure(2);
h1 = plot(number,10.*log10(FROG_error),'.','MarkerSize',10);
set(get(h1(1,:),'parent'),'linewidth',3,'fontsize',20);
text(50,0,'The minimum FROG error is ','fontsize',20);
text(100,0,[num2str(10.*log10(FROG_error(end)))],'fontsize',20,'color',[1,0,0]);
text(150,0,' dB','fontsize',20);
title('FROG Error G^{(k)}');
xlabel('Iteration Number k');
ylabel('G^{(k)} ( dB )');

t_new = linspace(-16,15.875,1e4);
a_t_new = interp1(t,a_t,t_new,'linear');

%% Find FWHM

I = abs(a_t_new).^2;
Ppeak = max(I);
half = 0.5*Ppeak;

move = length(t_new)./2 - find(I == Ppeak);
a_t_new = circshift(a_t_new,[0,move]);
I = abs(a_t_new).^2;

for m = 1:length(t_new)
    
    if I(m) > half
        index1 = m;
        break;
    end
    
end

for n = index1:length(t_new)
    
    if I(n) < half
        index2 = n;
        break;
    end
    
end

FWHM = abs(t_new(index2) - t_new(index1));

figure;
h2 = plot(t_new,I,'linewidth',3);
set(get(h2(1,:),'parent'),'linewidth',3,'fontsize',20);
title('Temperal Pulse Shape');
xlabel('Time ( ps )');
ylabel('Intensity ( a.u. )');
legend(['FWHM = ',num2str(FWHM),' ps']);
xlim([-4 4]);

figure;
contourf(tau,f,I_FROG,'linestyle','none');
axis([-4 4 -2 2]);
title('Retrieved PG FROG Trace I_{ FROG}(f,\tau)','fontsize',20);
xlabel('Delay \tau ( ps )','fontsize',20);
ylabel('Angular Frequency ( THz )','fontsize',20);

[w_new, A_new] = CFT(t_new, a_t_new);
phi = phase(A_new);% - max(phase(A_new));
f_new = w_new./(2.*pi); % 4937 ~ 5065

Apower = abs(A_new).^2;
Apeak = max(Apower);
half = 0.5*Apeak;

for m = 1:length(t_new)
    
    if Apower(m) > half
        index1 = m;
        break;
    end
    
end

for n = index1:length(t_new)
    
    if Apower(n) < half
        index2 = n;
        break;
    end
    
end

wFWHM = abs(w_new(index2) - w_new(index1))./(2.*pi);

figure;
[h3 spec wphase] = plotyy(f_new,Apower,f_new(4937:5065),phi(4937:5065));
title('Spectrum and Spectral Phase','fontsize',20);
set(get(h3(1),'ylabel'),'string','Power Spectrum |A( f )|^{ 2} ( a.u. )','fontsize',20);
set(get(h3(2),'ylabel'),'string','Phase \Psi ( f ) ( rad )','fontsize',20);
set(spec,'linewidth',3);
set(wphase,'linestyle','--','linewidth',3);
set(h3(1),'xlim',[-2,2],'fontsize',20,'LineWidth',3);
set(h3(1),'ylim',[0,4],'fontsize',20,'LineWidth',3);
set(h3(2),'xlim',[-2,2],'fontsize',20,'LineWidth',3);
%set(h3(2),'ylim',[0,-8000],'fontsize',20,'LineWidth',3);
xlabel('Frequency ( THz )','fontsize',20);
text(-1,3,['FWHM \Delta \nu = ',num2str(wFWHM),' Thz'],'fontsize',20);