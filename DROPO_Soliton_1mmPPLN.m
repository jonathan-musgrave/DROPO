run('Simulation_Parameters.m')

Nrt = 500;
detune = -6*alpha;
CCC_ss =  1.2420e-13;
    Q0 =1.8116e+07;
    Qfrac = Q0./Q;
PPump0 = 17.0E-3*Qfrac.^2;
    % PPump0 = 1.5E-3;
    % CCC = (PPump0).^2/(Q.^2*kappa.^2)
tseed = 2E-12;
APump0 = sqrt(PPump0)*ones(1,Nw);
ASignal0 = sqrt(PPump0)*exp(-2*sqrt(2)*(t/tseed).^2);
ASignal = zeros(1,Nw);
APump = zeros(1,Nw);
AAPump = zeros(Nrt,Nw);
AASignal = zeros(Nrt,Nw);

dk = -2*atan(2*detune/alpha)/Lnl;

for indrt = 1:1:Nrt
    indrt
    
    if indrt<50
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
      
run('run_plots.m')