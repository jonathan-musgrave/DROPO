clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----general definition-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c eps_0;  
c =  2.99792458E8;      % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;      % dielectric constant in vacuum, F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------GVM and GDD of PPLN--------%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------GV-------%%%%%%%%%
beta1_FF = 7421E-12;
beta1_SHG = 7421E-12 + 0E-12;
beta_offset4 = beta1_FF - beta1_FF;
beta_offset5 = beta1_SHG - beta1_FF;
%%%%%%%---------GVD-------%%%%%%%%%
beta2_FF = -325E-27;
beta2_SHG = 163E-27;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamFF = 2524E-9;               % wavelength of fudamental laser for SHG, m
lamSHG = lamFF/2;              % wavelength of SHG laser for SHG, m
global n_FF n_SHG
n_FF = beta1_FF * c;
n_SHG = beta1_SHG * c;

Anl = pi*(3E-6)^2;               % m^2, mode area in PPLN
deff = 14E-12;                   % nonlinear coefficient for PPLN,pm/V 
chi2_eff = deff/sqrt(Anl);       % nonlinear coefficient for PPLN,V 
w_FF = 2*pi*c/lamFF;                 
kappa = sqrt(2)*w_FF*chi2_eff/sqrt(c^3*n_FF^2*n_SHG*eps_0);
n2 = 2E-19;
gamma = w_FF*n2/c/Anl;     % for FF
Lnl = 1E-3;                % length of NL crystal, m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------Sim parameters---------%%%%%%%%%%%%%%%%%%%
Nw = 2001;             % slicing number for time
Nzcr= 21;              % slicing number for PPLN

tww = beta1_FF*Lnl*2;        % time window, related to the roundtrip time
t = linspace(-tww,tww,Nw);    
dt = mean(diff(t));                                                    
w = 2*pi*linspace(-1/2/dt,1/2/dt,Nw); % frequency 
dw = mean(diff(w));                                                                                                                        
zKTP = linspace(0,Lnl,Nzcr);   % PPLN2                            
dzKTP = mean(diff(zKTP));       % PPLN2                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alphact = pi/1600;% total loss
alphac = alphact/Lnl;
D_FF = exp(-alphac/2*(dzKTP/2)).*exp(1j*(beta_offset4.*w + beta2_FF/2.*w.^2)*dzKTP/2); 
D_SHG = exp(-alphac/2*(dzKTP/2)).*exp(1j*(beta_offset5.*w + beta2_SHG/2.*w.^2)*dzKTP/2);


RR1 = 1 - pi/1600;
alpha = ((1 - RR1) + alphact)/2;
detune = alpha*6.0;


Nrt = 3000;
PSHG0 = 5E-3;
tseed = 2E-12;
ASHG0 = sqrt(PSHG0)*ones(1,Nw);
AFF0 = sqrt(PSHG0)*exp(-2*sqrt(2)*(t/tseed).^2);
AFF = zeros(1,Nw);
ASHG = zeros(1,Nw);
AASHG = zeros(Nrt,Nw);
AAFF = zeros(Nrt,Nw);

dk = -2*atan(2*detune/alpha)/Lnl;

for indrt = 1:1:Nrt
    indrt
    
    if indrt<500
        AFF = sqrt(1-RR1)*AFF0 + sqrt(RR1)*AFF*exp(-1j*detune);
    else
        AFF =  sqrt(RR1)*AFF*exp(-1j*detune);
    end
    ASHG = sqrt(1-RR1)*ASHG0 + sqrt(RR1)*ASHG*exp(-1j*2*detune);
    
    for Z2 = zKTP  
     % Propagation (1st half), split-step Fourier method 
         sAFF = fftshift(ifft(ifftshift(AFF)));
         sASHG = fftshift(ifft(ifftshift(ASHG)));
         AFF = fftshift(fft(ifftshift(D_FF.*sAFF)));
         ASHG = fftshift(fft(ifftshift(D_SHG.*sASHG)));
        
     % nonlinear step using Runga-Kutta 4th order  
         [AFF, ASHG] = SHComb(AFF,ASHG,kappa,gamma,dk,Z2,dzKTP); 
         
     % Propagation (1st half), split-step Fourier method 
         sAFF = fftshift(ifft(ifftshift(AFF)));
         sASHG = fftshift(ifft(ifftshift(ASHG)));
         AFF = fftshift(fft(ifftshift(D_FF.*sAFF)));
         ASHG = fftshift(fft(ifftshift(D_SHG.*sASHG)));
    end
    
    AASHG(indrt,:) = ASHG;
    AAFF(indrt,:) = AFF;  
end
      
PFF = abs(AAFF).^2;
PSHG = abs(AASHG).^2;

AA1 = zeros(Nrt,1);
AA2 = zeros(Nrt,1);
for ind = 1:Nrt
    AA1(ind)= max(PFF(ind,:));
    AA2(ind)= max(PSHG(ind,:));
end
roundtrip=1:1:Nrt;
LW = 2;
FS = 20;
figure(1)
plot(roundtrip,AA1,'r','linewidth',LW)
hold on
plot(roundtrip,AA2,'g','linewidth',LW)
hold off
xlabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)


trange = 2;
ind = Nrt;
inst_w = -gradient(unwrap(angle(AAFF(ind,:))),dt)/2/pi/1E12;
figure(2)
yyaxis left
plot(t*1E12,PFF(ind,:),'r-','linewidth',LW);
hold on
plot(t*1E12,PSHG(ind,:),'g-','linewidth',LW);
hold off
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'xlim',[-trange trange])
%set(gca,'ylim',[0 0.6])
set(gca,'Ycolor','k')
%set(gca, 'YScale', 'log')

yyaxis right
plot(t*1E12,inst_w,'b--','linewidth',LW);
ylabel('chirp (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'Ycolor','b')
%set(gca,'ylim',[-0.8 0.8])
h = legend('signal','pump','chirp','location','northeast');legend boxoff 
set(h,'Fontsize',FS);


frange = 6;
figure(3)
SP = abs(fftshift(ifft(ifftshift(AAFF(ind,:))))).^2/max(abs(fftshift(ifft(ifftshift(AAFF(ind,:))))).^2);
plot(w/2/pi/1e12,10*log10(SP),'r','linewidth',LW)
hold on
SP1 = abs(fftshift(ifft(ifftshift(AASHG(ind,:))))).^2/max(abs(fftshift(ifft(ifftshift(AASHG(ind,:))))).^2);
plot(w/2/pi/1e12,10*log10(SP1),'g','linewidth',LW)
hold off
set(gca,'xlim',[-frange frange])
ylabel('power (dBm)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('frequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
h = legend('signal','pump','location','south');legend boxoff
set(h,'Fontsize',FS);clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----general definition-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c eps_0;  
c =  2.99792458E8;      % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;      % dielectric constant in vacuum, F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------GVM and GDD of PPLN--------%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------GV-------%%%%%%%%%
beta1_FF = 7421E-12;
beta1_SHG = 7421E-12 + 0E-12;
beta_offset4 = beta1_FF - beta1_FF;
beta_offset5 = beta1_SHG - beta1_FF;
%%%%%%%---------GVD-------%%%%%%%%%
beta2_FF = -325E-27;
beta2_SHG = 163E-27;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamFF = 2524E-9;               % wavelength of fudamental laser for SHG, m
lamSHG = lamFF/2;              % wavelength of SHG laser for SHG, m
global n_FF n_SHG
n_FF = beta1_FF * c;
n_SHG = beta1_SHG * c;

Anl = pi*(3E-6)^2;               % m^2, mode area in PPLN
deff = 14E-12;                   % nonlinear coefficient for PPLN,pm/V 
chi2_eff = deff/sqrt(Anl);       % nonlinear coefficient for PPLN,V 
w_FF = 2*pi*c/lamFF;                 
kappa = sqrt(2)*w_FF*chi2_eff/sqrt(c^3*n_FF^2*n_SHG*eps_0);
n2 = 2E-19;
gamma = w_FF*n2/c/Anl;     % for FF
Lnl = 1E-3;                % length of NL crystal, m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------Sim parameters---------%%%%%%%%%%%%%%%%%%%
Nw = 2001;             % slicing number for time
Nzcr= 21;              % slicing number for PPLN

tww= beta1_FF*Lnl*2;        % time window, related to the roundtrip time
t = linspace(-tww,tww,Nw);    
dt = mean(diff(t));                                                    
w = 2*pi*linspace(-1/2/dt,1/2/dt,Nw); % frequency 
dw = mean(diff(w));                                                                                                                        
zKTP = linspace(0,Lnl,Nzcr);   % PPLN2                            
dzKTP = mean(diff(zKTP));       % PPLN2                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alphact = pi/1600;% total loss
alphac = alphact/Lnl;
D_FF = exp(-alphac/2*(dzKTP/2)).*exp(1j*(beta_offset4.*w + beta2_FF/2.*w.^2)*dzKTP/2); 
D_SHG = exp(-alphac/2*(dzKTP/2)).*exp(1j*(beta_offset5.*w + beta2_SHG/2.*w.^2)*dzKTP/2);


RR1 = 1 - pi/1600;
alpha = ((1 - RR1) + alphact)/2;
detune = alpha*6.0;


Nrt = 5000;
PSHG0 = 5E-3;
tseed = 2E-12;
ASHG0 = sqrt(PSHG0)*ones(1,Nw);
AFF0 = sqrt(PSHG0)*exp(-2*sqrt(2)*(t/tseed).^2);
AFF = zeros(1,Nw);
ASHG = zeros(1,Nw);
AASHG = zeros(Nrt,Nw);
AAFF = zeros(Nrt,Nw);

dk = -2*atan(2*detune/alpha)/Lnl;

for indrt = 1:1:Nrt
    indrt
    
    if indrt<500
        AFF = sqrt(1-RR1)*AFF0 + sqrt(RR1)*AFF*exp(-1j*detune);
    else
        AFF =  sqrt(RR1)*AFF*exp(-1j*detune);
    end
    ASHG = sqrt(1-RR1)*ASHG0 + sqrt(RR1)*ASHG*exp(-1j*2*detune);
    
    for Z2 = zKTP  
     % Propagation (1st half), split-step Fourier method 
         sAFF = fftshift(ifft(ifftshift(AFF)));
         sASHG = fftshift(ifft(ifftshift(ASHG)));
         AFF = fftshift(fft(ifftshift(D_FF.*sAFF)));
         ASHG = fftshift(fft(ifftshift(D_SHG.*sASHG)));
        
     % nonlinear step using Runga-Kutta 4th order  
         [AFF, ASHG] = SHComb(AFF,ASHG,kappa,gamma,dk,Z2,dzKTP); 
         
     % Propagation (1st half), split-step Fourier method 
         sAFF = fftshift(ifft(ifftshift(AFF)));
         sASHG = fftshift(ifft(ifftshift(ASHG)));
         AFF = fftshift(fft(ifftshift(D_FF.*sAFF)));
         ASHG = fftshift(fft(ifftshift(D_SHG.*sASHG)));
    end
    
    AASHG(indrt,:) = ASHG;
    AAFF(indrt,:) = AFF;  
end
      
PFF = abs(AAFF).^2;
PSHG = abs(AASHG).^2;

AA1 = zeros(Nrt,1);
AA2 = zeros(Nrt,1);
for ind = 1:Nrt
    AA1(ind)= max(PFF(ind,:));
    AA2(ind)= max(PSHG(ind,:));
end
roundtrip=1:1:Nrt;
LW = 2;
FS = 20;
figure(1)
plot(roundtrip,AA1,'r','linewidth',LW)
hold on
plot(roundtrip,AA2,'g','linewidth',LW)
hold off
xlabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)


trange = 2;
ind = Nrt;
inst_w = -gradient(unwrap(angle(AAFF(ind,:))),dt)/2/pi/1E12;
figure(2)
yyaxis left
plot(t*1E12,PFF(ind,:),'r-','linewidth',LW);
hold on
plot(t*1E12,PSHG(ind,:),'g-','linewidth',LW);
hold off
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'xlim',[-trange trange])
%set(gca,'ylim',[0 0.6])
set(gca,'Ycolor','k')
%set(gca, 'YScale', 'log')

yyaxis right
plot(t*1E12,inst_w,'b--','linewidth',LW);
ylabel('chirp (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'Ycolor','b')
%set(gca,'ylim',[-0.8 0.8])
h = legend('signal','pump','chirp','location','northeast');legend boxoff 
set(h,'Fontsize',FS);


frange = 6;
figure(3)
SP = abs(fftshift(ifft(ifftshift(AAFF(ind,:))))).^2/max(abs(fftshift(ifft(ifftshift(AAFF(ind,:))))).^2);
plot(w/2/pi/1e12,10*log10(SP),'r','linewidth',LW)
hold on
SP1 = abs(fftshift(ifft(ifftshift(AASHG(ind,:))))).^2/max(abs(fftshift(ifft(ifftshift(AASHG(ind,:))))).^2);
plot(w/2/pi/1e12,10*log10(SP1),'g','linewidth',LW)
hold off
set(gca,'xlim',[-frange frange])
ylabel('power (dBm)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('frequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
h = legend('signal','pump','location','south');legend boxoff
set(h,'Fontsize',FS);