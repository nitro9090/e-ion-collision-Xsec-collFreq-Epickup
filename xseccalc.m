function [xsection]  = xseccalc(N,Trel,debyeL,redmass,nxsec)
% Trel is the relative temperature between the two colliding objects, since
% electrons travel much faster than ions, the electrons dominate this value

E0 = 8.85419E-12; %F/m, permittivity of free space
q = 1.6022e-19; %C, charge of a particle
k = 1.6022e-19; %J/eV

v = (2*Trel*k/redmass).^.5;  %m/s relative velocity between particles
bpi2 = q^2./(4*pi*E0*redmass*v.^2)*100;  %cm
xseclargeangle = pi*bpi2.^2;  %the large angle cross-section
xsection = (8*log(debyeL./bpi2)+1).*xseclargeangle;  % including the small angle coulomb collisions with the large angle collisions

[ArraysizeX, ArraysizeY] = size(N);

for X = 1:ArraysizeX
    for Y = 1:ArraysizeY
        if xsection(X,Y) <= nxsec
            xsection(X,Y) = nxsec;
        end
    end
end

