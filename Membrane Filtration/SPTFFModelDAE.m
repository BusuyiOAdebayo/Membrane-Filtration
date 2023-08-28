% SPTFFModel
function dcTFFdt = SPTFFModelDAE(~,cTFF,noOfComps,cTFFfeed,vdotfeed,Patm,L,Nz,d,H,k,rho,mu,De)
% cTFF and thus dcTFFdt here are column vectors of length noOfComps*Nz
% function dcTFFdt = SPTFFModel(~,cTFF,theta)
% noOfComps = theta(1)
% cTFFfeed = theta(2:noOfComps+1);
% vdotfeed = theta(noOfComps+2);
% Patm = theta(noOfComps+3);
% L = theta(noOfComps+4);
% Nz = theta(noOfComps+5);
% d = theta(noOfComps+6);
% H = theta(noOfComps+7);
% k = theta(noOfComps+8);
% rho = theta(noOfComps+9);
% mu = theta(noOfComps+10);
% De = theta(noOfComps+11:2*noOfComps+10);
% theta = [cTFFfeed,vdotfeed,Patm,L,Nz,d,H,k,rho,mu,De];

% Feed properties not same as inlet (z = 0) properties!
% Parameters to optimize: k, rho, mu, and De
% Mm_IgG = ?
Ax = pi*d^2/4;
Rmem = H/k; % Intrinsic resistance of the membrane, It is determine using DI water flux, H is membrane thickness, k is membrane permeability on DI water.
% H = membrane thickness [m]
% k = membrane permeability [m2]
Rrev = 0;%(m/A_cross)*alpha; % Reversible resistance due to concentration polarization (cake formation)
% alpha = alpha0*DeltaP_rad^k_compr
% alpha0 = 6.4e10; % [m/kg Pa^-0.268] A constant that depends on the cake particle size and shape
% k_compr = 0.268; % Cake compresibility index
Rirr = 0; %E1*(Re/CTFFA_wall*(PTMP-POsm))^E2;% Irreversible resistance due to fouling, E1 and E2 are experimental parameters, which must be determined for each filter type: capsule, cassette, spiral, etc
% E1 = 8694 g2s2/dm^5;
% E2 = -0.84;
% R2 = 0.97
% Alternatively, Rirr = 0;%E1*PTMP^E2;
% E1 = 5.31;
% E2 = 0.65;
% R2 = 0.99, implies this is more accurate.
g = 9.8; % m/s^2
% vdotinlet = Ax*uinlet; % Volumtric flow rate
% uinlet = vdotinlet/Ax; % Volumtric flow rate
Pinlet = 5*Patm; % Inlet/feed pressure
Pper = Patm;
% Pret = Patm;
% PTMB = (Pfeed+Pret)/2-Pper;
% Rg = 8.314; % [J/mol/K]
% T = 298.15; % [K]
% POsm = cTFFA*Rg*T; %Simplest is the van't Hoff approximation. Not good for any unusually high concentration operation, or where accuracy is important.
% A1 = 3.7e-1; A2 = -2.98e-3; A3 = 1e-5;
% POsm = 0.25*Patm; % A1*CTFFA+A2*CTFFA^2+A3*CTFFA^3; For Bovine Serum Albumin (BSA), Experimental fitting to estimate A1, A2 and A3
% B1 = 4.11e-2; B2 = -3.10e-4; B3 = 1.10e-6;
% viscosity_sol = B1*CTFFA+B2*CTFFA^2+B3*CTFFA^3; % For Bovine Serum Albumin (BSA)
% A1 = 7.0e-3; A2 = 2.60e-4; A3 = 4.0e-7;
POsm = 0.25*Patm; % Rg*T*(A1*CTFFA+A2*CTFFA^2+A3*CTFFA^3); For immunoglobulin G (IgG), Experimental fitting to estimate A1, A2 and A3
% CTFFA_max = 800g/L; B = 1.19; viscosity_sol_0 = viscosity_per = 1.75e-5; % This is equivalent to permeate viscosity, i.e, viscosity of the pure solvent, here we could use DI water
% viscosity_sol = viscosity_sol_0*exp(B*CTFFA/(1-CTFFA/CTFFA_max)); % For immunoglobulin G (IgG)
% Density_sol and Density_per
% DeltaP_rad = PTMB-POsm;

% Polarization modulus
% Polar_mod = CTFFA_wall/CTFFA_bulk = exp(Jr_A/kf); % note that CTFFA_bulk = CTFFA;

%     % Determination of mass transfer coefficient, kf
%      Re_crit = 2000;
%      ubulk = vdot/Ax;
%      Re = ubulk*density_sol*L/viscosity_sol;
%      Alternatively
%      Re = ueff*density_sol*d_h/viscosity_sol;
%      F1 = 1.620-1.664; F2 = 0.33; f1 = 0.33; f2 = 0.33;
%      De_IgG = 8.314e-8*T/viscosity_sol/Mm_A^1/3;
%      Sc = viscosity_sol/density_sol/De_IgG;
%      Sh = F1+F2*Sc^f1*Re^f2;
%      Alternatively, F1 = 0.082; f1 = 0.33; f2 = 0.69;
%      Sh = F1*Sc^f1*Re^f2;
%       if Re < Re_crit
%           kf = kf0*(wall_shear_rate*d^2/L)^1/3
%                 if conf == "Circular"
%                     wall_shear_rate = 8*ubulk/d
%                 elseif conf == "Rectangular"
%                     wall_shear_rate = 3*ubulk/d_h
%                 else
%                     wall_shear = 90;
%                 end
%       else
%           kf = Sh*D_A/L;
%       end

% Step for position
z = linspace(0,L,Nz);
dz = z(2)-z(1);

dcTFFdt = zeros(noOfComps*Nz,1);
P_axi = zeros(Nz,1);
vdot = zeros(Nz,1);
c_factor = 7; % cTFF(Nz)/cTFFfeed(1); % Volumetric concentration factor (VCF)

% Initial boundary equations (Conditions) for cTFF, P_axi, and vdot
for i = 1:noOfComps
    cTFF((i-1)*Nz+1) = cTFFfeed(i);
    % cTFF(i*Nz) = (4*cTFF(i*Nz-1)-cTFF(i*Nz-2))/3;
    if i == 1
        cTFF(i*Nz) = c_factor*cTFFfeed(i);
    else
        cTFF(i*Nz) = (4*cTFF(i*Nz-1)-cTFF(i*Nz-2))/3;
    end
end
P_axi(1) = Pinlet;
% P_axi(Nz) = (4*P_axi(Nz-1)-P_axi(Nz-2))/3;
vdot(1) = vdotfeed;
% vdot(Nz) = (4*vdot(Nz-1)-vdot(Nz-2))/3;

% Initalize dcTFFdz, d2cTFFdz2, dP_axidz, and dvdotdz
dcTFFdz = zeros(noOfComps*Nz,1);
d2cTFFdz2 = zeros(noOfComps*Nz,1);
% for i = 1:noOfComps*Nz
%     dcTFFdz(i) = zeros(noOfComps*Nz,1);
%     d2cTFFdz2(i) = zeros(noOfComps*Nz,1);
% end
dP_axidz = zeros(Nz,1);
dvdotdz = zeros(Nz,1);

% Equations for Jr: % Volumetric flux, [m3-volume-flow/s-time/m2-cross-sectional area
Jr = zeros(Nz,1);
for i = 1:Nz
    Jr(i) = (P_axi(i)-Pper-POsm)/(mu*(Rmem+Rrev+Rirr));
end

% Semi-interior equations for P_axi and vdot
% Rewrite in DAE format
for i = 2:Nz
    y(4) = dP_axidz(i)+(128*mu/pi/d^4)*vdot(i)-rho*g; % For cylindrical channel
    % dP_axidz(j) = -(3*mu/2/d^3/Ww)*vdot(j)+rho*g; % For rectangular channel
    % Alternative:
    % dP_axidz(j) = cd*density_sol*ueff^2/2/d_h;
    % cd = D1*Re^D2;
    % D1 and D2 are experimental parameters, which must be determined for each filter type: capsule, cassette, spiral, etc
    % ueff = vdot/Ax
    % Alternative
    % ueff = vdotinlet/Ax-Jr; % This assumes that the mass transport through the membrane leads to a reduced effective velocity while the volume of the membrane channel is constant
    % Alternative
    % ueff = vdot/Ax-Jr; % In my case here, this is the one to use bcos channel volume is not constant, i.e., volumetric flow rate is not constant and that Jr affects my volumetric flow rate
    y(5) = dvdotdz(i)+pi*d*Jr(i); % For cylindrical channel
    % dvdotdz(i) = -Nw*Ww*Jr(i); % For rectangular channel
end

% Interior equations for cTFF
for i = 1:noOfComps
    for j = (i-1)*Nz+2:i*Nz-1
        dcTFFdz(j) = (cTFF(j+1)-cTFF(j-1))/2/dz;
        d2cTFFdz2(j) = (cTFF(j+1)-2*cTFF(j)+cTFF(j-1))/dz/dz;
        dcTFFdt(j) = -(vdot(j-(i-1)*Nz)/Ax)*dcTFFdz(j)-(cTFF(j)/Ax)*dvdotdz(j-(i-1)*Nz)+De(i)*d2cTFFdz2(j)-4*Jr(j-(i-1)*Nz)*cTFF(j)/d;
        % % dmdt(i) = cTFFAdt(i)*vdot_per(i)-k_attr*vdot_feed(i)*m(i);
        % % k_attr0 = 300;
        % % k_attr1 = 10e9;
        % % k_attr = k_attr0; % k_attr0+k_attr1*vdot_per; Cake attrition factor
    end
end
end