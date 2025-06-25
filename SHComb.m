function   [AFF, ASHG] = SHComb(AFF,ASHG,kappa,gamma,dk,Z2,dzKTP)                        
        
        kFF_1 = 1j*kappa*(conj(AFF).*ASHG)*exp(-1j*dk*Z2) ...
                + 1j*gamma*abs(AFF).^2.*AFF;    
        kSHG_1 = 1j*kappa*(AFF.^2)*exp(1j*dk*Z2)...
                + 1j*2*gamma*abs(ASHG).^2.*ASHG;           
        
        %A_half2 = A + k1*dz/2;
        AFF_half2 = AFF + kFF_1*dzKTP/2; 
        ASHG_half2 = ASHG + kSHG_1*dzKTP/2;
        kFF_2 = 1j*kappa*(conj(AFF_half2).*ASHG_half2)*exp(-1j*dk*(Z2+dzKTP/2))...
                + 1j*gamma*abs(AFF_half2).^2.*AFF_half2;           
        kSHG_2 = 1j*kappa*(AFF_half2.^2)*exp(1j*dk*(Z2+dzKTP/2))...
                 + 1j*2*gamma*abs(ASHG_half2).^2.*ASHG_half2;
         
        %A_half3 = A + k2*dz/2;
        AFF_half3 = AFF + kFF_2*dzKTP/2;
        ASHG_half3 = ASHG + kSHG_2*dzKTP/2;
        kFF_3 = 1j*kappa*(conj(AFF_half3).*ASHG_half3)*exp(-1j*dk*(Z2+dzKTP/2))...
                + 1j*gamma*abs(AFF_half3).^2.*AFF_half3; 
        kSHG_3 = 1j*kappa*(AFF_half3.^2)*exp(1j*dk*(Z2+dzKTP/2))...
                 + 1j*2*gamma*abs(ASHG_half3).^2.*ASHG_half3;
        
        %A_full = A + k3*dz;
        AFF_full = AFF + kFF_3*dzKTP;
        ASHG_full = ASHG + kSHG_3*dzKTP;
        kFF_4 = 1j*kappa*(conj(AFF_full).*ASHG_full)*exp(-1j*dk*(Z2+dzKTP))...
                + 1j*gamma*abs(AFF_full).^2.*AFF_full;
        kSHG_4 = 1j*kappa*(AFF_full.^2)*exp(1j*dk*(Z2+dzKTP))...
                 + 1j*2*gamma*abs(ASHG_full).^2.*ASHG_full;
        
        %A = A + dz*(k1 + 2*k2 + 2*k3 + k4)/6;
        AFF = AFF + dzKTP*(kFF_1 + 2*kFF_2 + 2*kFF_3 + kFF_4)/6;
        ASHG = ASHG + dzKTP*(kSHG_1 + 2*kSHG_2 + 2*kSHG_3 + kSHG_4)/6;
        
end