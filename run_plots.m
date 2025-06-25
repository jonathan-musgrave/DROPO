%% Plots      
PSignal = abs(AASignal).^2;
PPump = abs(AAPump).^2;


AA1 = zeros(Nrt,1);
AA2 = zeros(Nrt,1);
for ind = 1:Nrt
    AA1(ind)= max(PSignal(ind,:));
    AA2(ind)= max(PPump(ind,:));
end
roundtrip=1:1:Nrt;
LW = 2;
FS = 20;
figure(1);clf;
plot(roundtrip,AA1,'r','linewidth',LW)
hold on
plot(roundtrip,AA2,'g','linewidth',LW)
hold off
xlabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('Peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)

figure(4);clf;
plot(roundtrip,sum(PSignal,2),'r','linewidth',LW)
hold on
plot(roundtrip,sum(PPump,2),'g','linewidth',LW)
hold off
xlabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('Average power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)

%
trange = 5;
ind = Nrt;

%%%%%%%%%% Specific plot ind %%%%%%%
% ind = 11175;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inst_w = -gradient(unwrap(angle(AASignal(ind,:))),dt)/2/pi/1E12;
figure(2);clf;
%yyaxis left
plot(t*1E12,PSignal(ind,:),'r-','linewidth',LW);
hold on
plot(t*1E12,PPump(ind,:),'g-','linewidth',LW);
hold off
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
ylabel('peak power (W)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%set(gca,'xlim',[-5 5])
%set(gca,'ylim',[-0.01 0.13],'ytick',0:0.1:0.2)
%set(gcf, 'position', [0 0 530 300]);
set(gca,'Ycolor','k')
%{
yyaxis right
plot(t*1E12,inst_w,'b-','linewidth',LW);
ylabel('FF chirp (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'Ycolor','b')
%}
legend('Signal','Pump','chirp');legend boxoff 

frange = 1;
figure(3);clf;
BV = -130;
SP = abs(fftshift(ifft(ifftshift(AASignal(ind,:))))).^2/max(abs(fftshift(ifft(ifftshift(AASignal(ind,:))))).^2);
stem(w/2/pi/1e12,10*log10(SP)+20,'r','linewidth',LW,'marker','none','basevalue',BV)
hold on
SP1 = abs(fftshift(ifft(ifftshift(AAPump(ind,:))))).^2/max(abs(fftshift(ifft(ifftshift(AAPump(ind,:))))).^2);
stem(w/2/pi/1e12,10*log10(SP1)-00,'g','linewidth',LW,'marker','none','basevalue',BV)
hold off
%set(gcf, 'position', [0 0 530 300]);
set(gca,'xlim',[-frange frange])
ylabel('power (dBm)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('frequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
legend('Signal','Pump');legend boxoff 
ylim([-50,24])


% figure(33)
% miu = c/lamSignal + w/2/pi;
% lambda = c./miu;
% plot(lambda*1E9,10*log10(SP),'r','linewidth',LW)
% set(gca,'ylim',[-80 0])
% 
% figure(44)
% miu_SH = c/lamPump + w/2/pi;
% lambda_SH = c./miu_SH;
% plot(lambda_SH*1E9,10*log10(SP1),'g','linewidth',LW)
% set(gca,'ylim',[-80 0])


%%%%%%%%%%%%%%%%%%%%%%%%%-----pulse evolutionfor each roundtrip-----%%%%%%%%%%%%%%%%%%%%%%%%%
roundtrip = 1:1:Nrt;
tt = t*1E12;
[TM,RTM] = meshgrid(tt,roundtrip);
tickrange_Nrt = 1000:1000:Nrt;
figure(6)
hh=pcolor(TM,RTM,PSignal);
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
%load seismic;
%colormap(seismic);
colorbar
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%set(gca,'ytick',1000:1000:5000);
set(gca,'ytick',tickrange_Nrt);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
title('Signal Temporal Evolution')
%set(get(hco,'Title'),'string','signal W');


%%%%%%%%%%%%%%%%%%%%%%%%%-----pulse evolutionfor each roundtrip-----%%%%%%%%%%%%%%%%%%%%%%%%%
roundtrip = 1:1:Nrt;
tt = t*1E12;
[TM,RTM] = meshgrid(tt,roundtrip);
tickrange_Nrt = 1000:1000:Nrt;
figure(7);clf;
hh=pcolor(TM,RTM,PPump);
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
%load seismic;
colormap(jet);
colorbar
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('time (ps)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
%set(gca,'ytick',1000:1000:5000);
set(gca,'ytick',tickrange_Nrt);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
%set(get(hco,'Title'),'string','signal W');
title('Pump Temporal Evolution')


% Spectral evolution Signal
roundtrip = 1:1:Nrt;
frange = 3E12;
wtickrange = floor((Nw+1)/2)-floor(frange/(dw/2/pi)):1:floor((Nw+1)/2)+floor(frange/(dw/2/pi));
ww = w(wtickrange)/2/pi/1E12;
[WM,RTM] = meshgrid(ww,roundtrip);
tickrange_Nrt = 1000:1000:Nrt;
SPEC = zeros(Nrt,Nw);
for ind = 1:Nrt
    
    SPEC(ind,:) = abs(fftshift(ifft(ifftshift(AASignal(ind,:)))))./max(abs(fftshift(ifft(ifftshift(AASignal(ind,:))))).^2);
end
% SPEC = SPEC./max(max(SPEC));

figure(8)
hh=pcolor(WM,RTM,10*log10(SPEC(roundtrip,wtickrange)));
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('rrequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'ytick',tickrange_Nrt);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
title('Signal Spectrum Evolution')
set(gca,'clim',[-80,0])
% Spectral evolution Signal
roundtrip = 1:1:Nrt;
frange = 3E12;
wtickrange = floor((Nw+1)/2)-floor(frange/(dw/2/pi)):1:floor((Nw+1)/2)+floor(frange/(dw/2/pi));
ww = w(wtickrange)/2/pi/1E12;
[WM,RTM] = meshgrid(ww,roundtrip);
tickrange_Nrt = 1000:1000:Nrt;
SPEC = zeros(Nrt,Nw);
for ind = 1:Nrt
    SPEC(ind,:) = abs(fftshift(ifft(ifftshift(AAPump(ind,:))))).^2./max(abs(fftshift(ifft(ifftshift(AASignal(ind,:))))).^2);
end
% SPEC = SPEC./max(max(SPEC));
figure(9);clf;
hh=pcolor(WM,RTM,10*log10(SPEC(roundtrip,wtickrange)));
set(hh,'edgecolor','k','Marker','*')
shading interp
set(gca,'LineWidth',LW,'FontSize',FS)
ylabel('roundtrip','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
xlabel('rrequency (THz)','FontName','Times New Roman','FontSize',FS,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold','linewidth',LW)
set(gca,'ytick',tickrange_Nrt);
hco=colorbar;
set(hco,'FontName','Times New Roman','FontSize',FS,'FontWeight','bold');
title('Signal Spectrum Evolution')
set(gca,'clim',[-60,-20])
nd = 2.^9;
D1 = beta_offset5*Lnl./alpha;
D2 = beta2_Signal*Lnl./alpha;
delta_2 = linspace(-12,12,nd);
[Omega,Delta_2] = meshgrid(w,delta_2);
X = (alpha*(1+(Delta_2-D1.*Omega-D2*Omega.^2).^2)).^(-1);
Y = (Delta_2-D1.*Omega-D2*Omega.^2).*X;
for pp = 1:length(delta_2)
    X(pp,:) = X(pp,:)./max(X(pp,:));
    Y(pp,:) = normalize(Y(pp,:),'range',[-1,1]);
    % Y(i,:) = Y(i,:)./max(Y(i,:));

end
figure(10);clf;
subplot(1,2,1)
imagesc(w./(2*pi*1e12),delta_2./2,X)
title('X(\nu)')
xlim([-frange,frange]./1e12);colorbar
xlabel('\nu')
ylabel('\delta_1')
subplot(1,2,2)
imagesc(w./(2*pi*1e12),delta_2./2,Y)
title('Y(\nu)')
xlim([-frange,frange]./1e12)
xlabel('\nu')
ylabel('\delta_1')
colorbar

