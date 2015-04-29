clear
clc
close all
pause on

%constants
%E0 = 8.85419E-12; %F/m, permittivity of free space
%q = 1.6022e-19; %C, charge of a particle
k = 1.6022e-19; %eV/(degrees C)
Z = 1;  % charge state
lnA = 10; % coulomb logarithm
mp = 1.6726e-27; %kg, proton mass
me = 9.1094e-31; %kg, electron mass
nxsec = 3e-16;  %cm^2, neutral collision cross-section from Bellan

%Input values
mew = 1;  %ave ratio of ion mass to proton mass
Efield = 30000; % V/cm, electric field
Nset = [1e13, 1e14, 1e15, 1e16]; %cm^-3, electron and ion densities
nNset = 2.44627E+19;  % cm^-3 neutral density
Epart = .1*2; %the energy of particles immediately after ionization
setTi = .033; %eV, ion Temperature

%font sizes
Legendsize = 12;
AxisSize = 12;
xLabelSize = 13;
yLabelSize = 13;
TitleSize = 13;


%Electron temperature range in eV
Tstart = .01; %first Te value 
Tstepsize = .01;  %the step size of each Te step
Tend = 10;  %last Te value 

Tintstart = 0;
Tintend = 50;

% calculated valuesmew
Teset = flipud(rot90(Tstart:Tstepsize:Tend)); %sets up Temperature array
Teintset = flipud(rot90(Tintstart:Tstepsize:Tintend));
debyeL = 7.43*10^2*Teset.^.5*Nset.^(-.5);  %cm, debye length
debyeLint = 7.43*10^2*Teintset.^.5*Nset.^(-.5);
mion = mew*mp;  %kg, ion mass 
redmei = mion*me/(mion+me); %kg, reduced mass between an electron and ion
redmii = (mion)^2/(2*mion);  %kg, reduced mass between two ions
redmee = me^2/(2*me);  %kg, reduced mass between two electrons

[~, ArraysizeX] = size(Nset);
[ArraysizeY,~] = size(Teset);
[Arraysizeint,~] = size(Teintset);

%setting up variables
nN = zeros(ArraysizeY, ArraysizeX);
N = zeros(ArraysizeY, ArraysizeX);
nxsec3 = zeros(ArraysizeY, ArraysizeX);
Ti = ones(ArraysizeY, ArraysizeX)*setTi;
Tiint = ones(Arraysizeint, ArraysizeX)*setTi;

for K = 1:ArraysizeX
    nN2(:,K) = nNset;
    nxsec2(:,K) = nxsec;
    Te(:,K) = Teset;
    Teint(:,K) = Teintset;
end

for K = 1:ArraysizeY
    nxsec4(K,:) = nxsec2;
    N(K,:) = Nset;
    nN(K,:) = nN2;
end

for K = 1:Arraysizeint
    N2(K,:) = Nset;
    nxsec3(K,:) = nxsec2;
end

%calculating the cross-sections
[iixsec] = xseccalc(N,Ti,debyeL,redmii,nxsec);
[iexsec] = xseccalc(N,Te,debyeL,redmei,nxsec);
[eixsec] = xseccalc(N,Te,debyeL,redmei,nxsec);
[eexsec] = xseccalc(N,Te,debyeL,redmee,nxsec);

[iixsecint] = xseccalc(N2,Tiint,debyeLint,redmii,nxsec);
[iexsecint] = xseccalc(N2,Teint,debyeLint,redmei,nxsec);
[eixsecint] = xseccalc(N2,Teint,debyeLint,redmei,nxsec);
[eexsecint] = xseccalc(N2,Teint,debyeLint,redmee,nxsec);

[iicollfreq, iiEpkup] = calcEpickup(iixsecint, N, Ti,Tiint,Tstepsize,redmii,Efield,Epart,mion);
[iecollfreq, ieEpkup] = calcEpickup(iexsecint, N, Te,Teint,Tstepsize,redmei,Efield,Epart,mion);
[eicollfreq, eiEpkup] = calcEpickup(eixsecint, N, Te,Teint,Tstepsize,redmei,Efield,Epart,me);
[eecollfreq, eeEpkup] = calcEpickup(eexsecint, N, Te,Teint,Tstepsize,redmee,Efield,Epart,me);
[encollfreq, enEpkup] = calcEpickup(nxsec3, nN, Te,Teint,Tstepsize,redmei,Efield,Epart,me);
[incollfreq, inEpkup] = calcEpickup(nxsec3, nN, Ti,Tiint,Tstepsize,redmii,Efield,Epart,mion);

% [iicollfreq2, iiEpkup2] = calcEpickup2(iixsecint, N, Ti,Tiint,Tstepsize,redmii,Efield,Epart,mion);
% [iecollfreq2, ieEpkup2] = calcEpickup2(iexsecint, N, Te,Teint,Tstepsize,redmei,Efield,Epart,mion);
% [eicollfreq2, eiEpkup2] = calcEpickup2(eixsecint, N, Te,Teint,Tstepsize,redmei,Efield,Epart,me);
% [eecollfreq2, eeEpkup2] = calcEpickup2(eexsecint, N, Te,Teint,Tstepsize,redmee,Efield,Epart,me);
% [encollfreq2, enEpkup2] = calcEpickup2(nxsec3, nN, Te,Teint,Tstepsize,redmei,Efield,Epart,me);
% [incollfreq2, inEpkup2] = calcEpickup2(nxsec3, nN, Ti,Tiint,Tstepsize,redmii,Efield,Epart,mion);

iTotcollfreq = iicollfreq + iecollfreq + incollfreq;
eTotcollfreq = eicollfreq + eecollfreq + encollfreq;

iTotEpkup = 1./(1./iiEpkup + 1./ieEpkup + 1./inEpkup);
eTotEpkup = 1./(1./eiEpkup + 1./eeEpkup + 1./enEpkup);

% iTotcollfreq2 = iicollfreq2 + iecollfreq2 + incollfreq2;
% eTotcollfreq2 = eicollfreq2 + eecollfreq2 + encollfreq2;
% 
% iTotEpkup2 = 1./(1./iiEpkup2 + 1./ieEpkup2 + 1./inEpkup2);
% eTotEpkup2 = 1./(1./eiEpkup2 + 1./eeEpkup2 + 1./enEpkup2);

%Plotting the graphs
figure('name', 'collision cross section, cm^2')
set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 1;0 0 0],'DefaultAxesLineStyleOrder','-|--|:');

hold('all');

for P = 1:ArraysizeX
    subplot(2,ArraysizeX,P)
    semilogy(Te(:,1), iixsec(:,P),'-', Te(:,1), iexsec(:,P),'--', Te(:,1), nxsec4(:,P),':','LineWidth',2)
%     semilogy(Te(:,1), iixsec(:,P), Te(:,1), iexsec(:,P), Te(:,1), nxsec3(:,P), 'g')
    title(['N_e=10^{',(num2str(log10(N(1,P)))),'} cm^{-3}'],'fontSize', TitleSize)
    axis([Tstart, Tend, 10^-16, 10^-8])
    if P == 1
        ilegend = {'ion-ion','ion-e','ion-neut'};
        legend(ilegend,'Location','northwestoutside','fontSize', Legendsize)
        ylabel('Ion Collisional X-section, cm^2','fontSize', yLabelSize)
    end
    subplot(2,ArraysizeX,P+ArraysizeX)
    semilogy(Te(:,1), eexsec(:,P),'-', Te(:,1), eixsec(:,P),'--', Te(:,1), nxsec4(:,P),':','LineWidth',2)
%     loglog(Te, eexsec(:,P), Te, eixsec(:,P), Te, nxsec3(:,P))
    axis([Tstart, Tend, 10^-16, 10^-8])
    xlabel('T_e (eV)','fontSize', xLabelSize)
    if P == 1
        elegend = {'e-e','e-ion','e-neut'};
        legend(elegend,'Location','southwestoutside','fontSize', Legendsize)
        ylabel('e- Collisional X-section, cm^2','fontSize', yLabelSize)
    end
end

% axis tight
% set(gca,'looseinset',[0 0 0 0])
% set(gcf,'units','normalized','outerposition',[0 0 1 1])

cd C:\Users\Rufus\Dropbox\Dissertation\DissPlots

Filename = 'CollXsec'
set(gcf,'PaperPositionMode','auto')
Print(gcf, '-r300', '-dpng',Filename)

figure('name','collision frequency')

multiplier = 1;

for P = 1:ArraysizeX
    subplot(2,ArraysizeX,P)
    semilogy(Te(:,1), iicollfreq(:,P),'-', Te(:,1), iecollfreq(:,P),'--', Te(:,1), incollfreq(:,P),':', Te(:,1), iTotcollfreq(:,P),'-.','LineWidth',2)
    title(['N_e=10^{',num2str(log10(N(1,P))),'} cm^{-3}'],'fontSize', TitleSize)
    axis([Tstart, Tend, 10^8, 10^13])
    set(gca,'FontSize',AxisSize)
    if P == 1
        ilegend = {'ion-ion','ion-e','ion-neut','tot ion coll freq'};
        legend(ilegend,'Location','northwestoutside','fontSize', Legendsize)
        ylabel('Ion Collision Freq, Hz','fontSize', yLabelSize)
    end
    subplot(2,ArraysizeX,P+ArraysizeX)
    semilogy(Te(:,1), eecollfreq(:,P),'-', Te(:,1), eicollfreq(:,P),'--', Te(:,1), encollfreq(:,P),':', Te(:,1), eTotcollfreq(:,P),'-.','LineWidth',2)
    set(gca,'FontSize',AxisSize)
    if P == 1
        elegend = {'e-e','e-ion','e-neut','tot e coll freq'};
        legend(elegend,'Location','southwestoutside','fontSize', Legendsize)
        ylabel('e- Collision Freq, Hz','fontSize', yLabelSize)
    end
%    axis([Tstart, Tend, 10^8, 10^13])
    xlabel('T_e (eV)','fontSize', xLabelSize)
end
% axis tight
% set(gca,'looseinset',[0 0 0 0])
% set(gcf,'units','normalized','outerposition',[0 0 1 1])

Filename = 'CollFreq'
set(gcf,'PaperPositionMode','auto')
Print(gcf, '-r300', '-dpng',Filename)

% for P = 1:ArraysizeX
%     subplot(2,ArraysizeX,P)
%     hold on
%     semilogy(Te(:,1), iicollfreq2(:,P)*multiplier, Te(:,1), iecollfreq2(:,P)*multiplier, Te(:,1), incollfreq2(:,P)*multiplier, Te(:,1), iTotcollfreq2(:,P)*multiplier,'linewidth',1.5)
%     hold off
%     subplot(2,ArraysizeX,P+ArraysizeX)
%     hold on
%     semilogy(Te(:,1), eecollfreq2(:,P)*multiplier, Te(:,1), eicollfreq2(:,P)*multiplier, Te(:,1), encollfreq2(:,P)*multiplier, Te(:,1), eTotcollfreq2(:,P)*multiplier,'linewidth',1.5)
%     hold off
% end

figure('name','Average particle energy pickup before collision, eV')
for P = 1:ArraysizeX
    subplot(2,ArraysizeX,P)   
%     loglog(Te, iiEpkup(:,P), Te, ieEpkup(:,P), Te, inEpkup(:,P),Te, iTotEpkup(:,P))
    semilogy(Te(:,1), iiEpkup(:,P),'-', Te(:,1), ieEpkup(:,P),'--', Te(:,1), inEpkup(:,P),':',Te(:,1), iTotEpkup(:,P),'-.','LineWidth',2)
    title(['N_e=10^{',(num2str(log10(N(1,P)))),'} cm^{-3}'],'fontSize', TitleSize)
    axis([Tstart, Tend, 10^-3, 10^4])
    set(gca,'FontSize',AxisSize)
    if P == 1
        ilegend = {'ion-ion','ion-e','ion-neut','ave ion E pickup'};
        ylabel('Ion Ave E Pickup, eV','fontSize', yLabelSize)
        legend(ilegend,'Location','northwestoutside','fontSize', Legendsize)
    end
    
    subplot(2,ArraysizeX,P+ArraysizeX)  
%     loglog(Te, eeEpkup(:,P), Te, eiEpkup(:,P), Te, enEpkup(:,P),Te, eTotEpkup(:,P))
    semilogy(Te(:,1), eeEpkup(:,P),'-', Te(:,1), eiEpkup(:,P),'--', Te(:,1), enEpkup(:,P),':',Te(:,1), eTotEpkup(:,P),'-.','LineWidth',2)
    axis([Tstart, Tend, 10^-3, 10^4])
    xlabel('T_e (eV)','fontSize', xLabelSize)
    set(gca,'FontSize',AxisSize)
    if P == 1
        elegend = {'e-e','e-ion','e-neut','ave e E pickup'};
        legend(elegend,'Location','southwestoutside','fontSize', Legendsize)
        ylabel('e- Ave E Pickup, eV','fontSize', yLabelSize)
    end
end
%axis tight
%set(gca,'looseinset',[0 0 0 0])
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

Filename = 'EnergyPickUp'
set(gcf,'PaperPositionMode','auto')
Print(gcf, '-r300', '-dpng',Filename)

% for P = 1:ArraysizeX
%     subplot(2,ArraysizeX,P)
%     hold on
% %     loglog(Te, iiEpkup(:,P), Te, ieEpkup(:,P), Te, inEpkup(:,P),Te, iTotEpkup(:,P))
%     semilogy(Te(:,1), iiEpkup2(:,P), Te(:,1), ieEpkup2(:,P), Te(:,1), inEpkup2(:,P),Te(:,1), iTotEpkup2(:,P),'linewidth',1.5)
%     hold off
% 
%     subplot(2,ArraysizeX,P+ArraysizeX)
%     hold on
% %     loglog(Te, eeEpkup(:,P), Te, eiEpkup(:,P), Te, enEpkup(:,P),Te, eTotEpkup(:,P))
%     semilogy(Te(:,1), eeEpkup2(:,P), Te(:,1), eiEpkup2(:,P), Te(:,1), enEpkup2(:,P),Te(:,1), eTotEpkup2(:,P),'linewidth',1.5)
%     hold off
% end

pause off