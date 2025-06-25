clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----general definition-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c eps_0;  
c =  2.99792458E8;      % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;      % dielectric constant in vacuum, F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------GVM and GDD of PPLN--------%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------GV-------%%%%%%%%%
beta1_Signal = 7421E-12;
beta1_Pump = 7421E-12 + 0E-12;
beta_offset4 = beta1_Signal - beta1_Signal;
beta_offset5 = beta1_Pump - beta1_Signal;
%%%%%%%---------GVD-------%%%%%%%%%
beta2_Signal = 104E-27*50;
beta2_Pump = 400E-27*50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamSignal = 1571.9E-9;               % wavelength of fudamental laser for SHG, m
lamPump = lamSignal/2;              % wavelength of SHG laser for SHG, m
global n_Signal n_Pump
n_Signal = beta1_Signal * c;
n_Pump = beta1_Pump * c;

Anl = pi*(3E-6)^2;               % m^2, mode area in PPLN
deff = 14E-12;                   % nonlinear coefficient for PPLN,pm/V 
chi2_eff = deff/sqrt(Anl);       % nonlinear coefficient for PPLN,V 
w_FF = 2*pi*c/lamSignal;                 
kappa = sqrt(2)*w_FF*chi2_eff/sqrt(c^3*n_Signal^2*n_Pump*eps_0);
n2 = 2E-19;
gamma = w_FF*n2/c/Anl;     % for FF
Lnl = 10E-3;                % length of NL crystal, m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------Sim parameters---------%%%%%%%%%%%%%%%%%%%
Nw = 500*2;             % slicing number for time
Nzcr= 21*2;              % slicing number for PPLN

tww= beta1_Signal*Lnl;        % time window, related to the roundtrip time
h = (-Nw/2:Nw/2-1);
dt = tww/Nw;
t = h.*dt;
dw = 2*pi/dt./Nw;
w = dw.*h;                                                                                                                
zKTP = linspace(0,Lnl,Nzcr);   % PPLN2                            
dzKTP = mean(diff(zKTP));      % PPLN2        
alphact = pi/1600;% total loss
alphac = alphact/Lnl;
D_Signal = exp(-alphac/2*(dzKTP/2)).*exp(1j*(beta_offset4.*w + beta2_Signal/2.*w.^2)*dzKTP/2); 
D_Pump = exp(-alphac/2*(dzKTP/2)).*exp(1j*(beta_offset5.*w + beta2_Pump/2.*w.^2)*dzKTP/2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR1 = 1 - pi/1600*100;
alpha = ((1 - RR1) + alphact)/2;

Finnese = 2*pi./(alpha);
FSR = c./(2*Lnl*n_Pump);
w_pump = c/(lamPump);
Q = Finnese*w_pump./FSR;