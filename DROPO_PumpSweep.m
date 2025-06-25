run('Simulation_Parameters.m')


%%%%%%%%%%%%
linewidth = alpha./tww./(4*pi); % loss/second == linewidth 
sweep_speed = linewidth./(100); % linewidths per Nrt
detune_win = alpha./tww./(4*pi)*([6,-7]);
Nrt = fix(abs(ceil(diff(detune_win))./(sweep_speed)))+1;
detune_ar = linspace(detune_win(1),detune_win(2),Nrt)./(alpha./tww./(4*pi)).*alpha;
noise = 1e-10; % amplitude of injected noise
%%%%%%%%%%%%

    
CCC_ss =  1.2420e-13;
    Q0 =1.8116e+07;
    Qfrac = Q0./Q;
PPump0 = 17.5E-3*Qfrac.^2;
PPump0 = 9.5E-3*Qfrac.^2;
    % PPump0 = 1.5E-3;
    % CCC = (PPump0).^2/(Q.^2*kappa.^2)
RNGs = [1:20];
AAPump_RNG = cell([length(RNGs),1]);
AASignal_RNG = cell([length(RNGs),1]);
for i = 1:length(RNGs)
    rng(RNGs(i))
    ASignal = zeros(1,Nw)+noise.*exp(-2*pi*1j*rand(1,Nw));
    tseed = 2E-12;
    APump0 = sqrt(PPump0)*ones(1,Nw);
    ASignal0 = sqrt(PPump0)*exp(-2*sqrt(2)*(t/tseed).^2);
    APump = zeros(1,Nw);
    AAPump = zeros(Nrt,Nw);
    AASignal = zeros(Nrt,Nw);
    
    for indrt = 1:1:Nrt
        detune = detune_ar(indrt);
        dk = -2*atan(2*detune/alpha)/Lnl*1;
    
        indrt
        
        if indrt<0
            ASignal = sqrt(1-RR1)*ASignal0 + sqrt(RR1)*ASignal*exp(-1j*detune);
        else
            ASignal =  sqrt(RR1)*ASignal*exp(-1j*detune);
        end
        APump = sqrt(1-RR1)*APump0 + sqrt(RR1)*APump*exp(-1j*2*detune);
        
        for Z2 = zKTP  
         % Propagation (1st half), split-step Fourier method 
             sASignal = fftshift(ifft(ifftshift(ASignal)));
             sAPump = fftshift(ifft(ifftshift(APump)));
             ASignal = fftshift(fft(ifftshift(D_Signal.*sASignal)));
             APump = fftshift(fft(ifftshift(D_Pump.*sAPump)));
            
         % nonlinear step using Runga-Kutta 4th order  
             [ASignal, APump] = SHComb(ASignal,APump,kappa,gamma,dk,Z2,dzKTP); 
             
         % Propagation (1st half), split-step Fourier method 
             sASignal = fftshift(ifft(ifftshift(ASignal)));
             sAPump = fftshift(ifft(ifftshift(APump)));
             ASignal = fftshift(fft(ifftshift(D_Signal.*sASignal)));
             APump = fftshift(fft(ifftshift(D_Pump.*sAPump)));
        end
        
        AAPump(indrt,:) = APump;
        AASignal(indrt,:) = ASignal;  
    end
    AAPump_RNG{i} = AAPump;
    AASignal_RNG{i} = AASignal;
    run('run_plots.m')
    f = figure(6);
    name = strcat('Signal Evolution',' RNG ',num2str(RNGs(i)));
    title(name)
    saveas(f,fullfile('RNG_Sweep',strcat(name,'.png')));

end

save("RNG_Sweep.mat",'AAPump_RNG','AASignal_RNG','t','w','detune_ar','alpha')