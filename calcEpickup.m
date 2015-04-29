function [collfreq, aveEpkup] = calcEpickup(xsec,N, Trel, Tint, Tstepsize, redmass, Efield, Epart,mass)
% this function calculates the collision frequency between particles and
% using that determines the average energy that a particle gains in an
% electric field given the collision frequency.

k = 1.6022e-19; %J/eV
q = 1.6022e-19; %Coulombs per particle
standardize =0;
xsectot = 0;

[TrelsizeY,TrelsizeX] = size(Trel);
[IntsizeY,~] = size(Tint);

vrel = (2*Tint*k/redmass).^.5*100;  %relative velocity of particles colliding

for X = 1:TrelsizeX
   for Y = 1:TrelsizeY
       for Z = 1:IntsizeY
           distfunction = exp(-(Tint(Z,X))^2/(Trel(Y,X))^2)*Tstepsize; %Druyvesteyn distribution
           standardize = standardize + distfunction;
           xsectot = xsectot + vrel(Z,X)*xsec(Z,X)*distfunction;
       end
       xsecfin(Y,X) = xsectot/standardize;
       standardize = 0;
       xsectot = 0;
   end
end

collfreq = xsecfin.*N; %particle collision frequency

vpart = (2*Epart*k/mass).^.5*100; %cm/s, the particle velocity after ionization

aveEpkup = q^2*Efield^2.*collfreq.^-2/(2*mass*k)*100^2+vpart*Efield./collfreq;