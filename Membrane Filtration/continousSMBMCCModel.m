% PC Model
function dcqPCdt = continousSMBMCCModel(~,cqPC,noOfComps,cPCfeed,vdotfeed,Patm,L,Nz,d,evoid,qs,Kc,KL,mu,rhob,De)
% cqPC and thus dcqPCdt here are column vectors of length noOfComps*Nz
% function dcqPCdt = continousPCModel(~,cqPC,theta)
% noOfComps = theta(1)
% cPCfeed = theta(noOfComps);
% vdotfeed = theta(noOfComps+2);
% Patm = theta(noOfComps+3);
% L = theta(noOfComps+4);
% Nz = theta(noOfComps+5);
% d = theta(noOfComps+6);
% evoid = theta(noOfComps+6);
% qs = theta(noOfComps+6);
% Kc = theta(noOfComps+6);
% KL = theta(noOfComps+6);
% mu = theta(noOfComps+9);
% rhob = theta(noOfComps+9);
% De = theta(noOfComps+11:2*noOfComps+10);
% theta = [cPCfeed,vdotfeed,Patm,L,Nz,r,evoid,qs,Kc,KL,rhob,De];

% Feed properties not same as inlet (z = 0) properties!
% Parameters to optimize: Kc, KL, mu, rhob and De
% Mm_IgG = ?
Ax = pi*d^2/4;
g = 9.8; % m/s^2
% vdotinlet = Ax*uinlet; % Volumtric flow rate
% uinlet = vdotinlet/Ax; % Volumtric flow rate
Pinlet = 5*Patm; % Inlet/feed pressure
% Rg = 8.314; % [J/mol/K]
% T = 298.15; % [K]
% viscosity_sol = B1*cPC_IgG+B2*cPC_IgG^2+B3*cPC_IgG^3; % For Bovine Serum Albumin (BSA)
% B1 = 4.11e-2; B2 = -3.10e-4; B3 = 1.10e-6;
% viscosity_sol = viscosity_sol_0*exp(B*cPC_IgG/(1-cPC_IgG/cPC_IgG_max)); % For immunoglobulin G (IgG)
% cPC_IgG_max = 800g/L; B = 1.19; viscosity_sol_0 = viscosity_per = 1.75e-5; % This is equivalent to permeate viscosity, i.e, viscosity of the pure solvent, here we could use DI water
% Density_sol = density of solution

% % Determination of mass transfer coefficient, kf
% Re_crit = 2000;
% ubulk = vdot/Ax;
% Re = ubulk*density_sol*L/viscosity_sol;
% Alternative:
% Re = ueff*density_sol*d_h/viscosity_sol;
% F1 = 1.620-1.664; F2 = 0.33; f1 = 0.33; f2 = 0.33;
% De_IgG = 8.314e-8*T/viscosity_sol/Mm_A^1/3;
% Sc = viscosity_sol/density_sol/De_IgG;
% Sh = F1+F2*Sc^f1*Re^f2;
% Alternative:
% F1 = 0.082; f1 = 0.33; f2 = 0.69;
% Sh = F1*Sc^f1*Re^f2;
% if Re < Re_crit
%     kf = kf0*(wall_shear_rate*d^2/L)^1/3
%     if conf == "Circular"
%         wall_shear_rate = 8*ubulk/d
%     elseif conf == "Rectangular"
%         wall_shear_rate = 3*ubulk/d_h
%     else
%         wall_shear = 90;
%     end
% else
%     kf = Sh*D_A/L;
% end

% Step for position
z = linspace(0,L,Nz);
dz = z(2)-z(1);

% Initialize dcPCdt, dqPCdt, dvdotdt, P_axi and vdot
dcPCdt = zeros(noOfComps*Nz,1);
dqPCdt = zeros(noOfComps*Nz,1);
dvdotdt = zeros(Nz,1);
P_axi = zeros(Nz,1);
vdot = zeros(Nz,1);

% Define values
cPC = cqPC(1:noOfComps*Nz);
qPC = cqPC(noOfComps*Nz+1:2*noOfComps*Nz);

% Initial boundary equations (conditions) for cPC, P_axi, and vdot
for i = 1:noOfComps
    cPC((i-1)*Nz+1) = cPCfeed(i);
    cPC(i*Nz) = (4*cPC(i*Nz-1)-cPC(i*Nz-2))/3;
end
% P_axi(1) = Pinlet;
% P_axi(Nz) = (4*P_axi(Nz-1)-P_axi(Nz-2))/3;
vdot(1) = vdotfeed;
vdot(Nz) = (4*vdot(Nz-1)-vdot(Nz-2))/3;

% Initalize dcqPCdz, d2cqPCdz2, dvdotdz, d2vdotdz2 and dP_axidz
dcPCdz = zeros(noOfComps*Nz,1);
d2cPCdz2 = zeros(noOfComps*Nz,1);
dvdotdz = zeros(Nz,1);
d2vdotdz2 = zeros(Nz,1);
% for i = 1:noOfComps*Nz
%     dcqPCdz(i) = zeros(noOfComps*Nz,1);
%     d2cqPCdz2(i) = zeros(noOfComps*Nz,1);
% end
dP_axidz = zeros(Nz,1);

% Equations for P_axi
for i = i:Nz
    if i == 1
        P_axi(i) = Pinlet;
    else
        dP_axidz(i) = -150*(1-evoid)^2*mu*vdot/evoid/Ax/evoid^3/dp^2-1.75*(1-evoid)*rhob*(vdot/evoid/Ax)^2/evoid^3/dp+rhob*g;
    end
end

% Interior equations for vdot
for i = 2:Nz-1
    %     dudt(i) = -dP_axidz(i)-u(i)*dudz(i)+mu*d2udz2(i)+rhob*g;
    dvdotdz(i) = (vdot(i+1)-vdot(i-1))/2/dz;
    d2vdotdz2(i) = (vdot(i+1)-2*vdot(i)+vdot(i-1))/dz/dz;
    dvdotdt(i) = -dP_axidz(i)/rhob*evoid*Ax-vdot(i)*dvdotdz(i)*evoid*Ax+mu/rhob*d2vdotdz2(i)*evoid*Ax+g*evoid*Ax;
end

% Interior equations for cPC
for i = 1:noOfComps
    for j = (i-1)*Nz+2:i*Nz-1
        dcPCdz(j) = (cPC(j+1)-cPC(j-1))/2/dz;
        d2cPCdz2(j) = (cPC(j+1)-2*cPC(j)+cPC(j-1))/dz/dz;
        dcPCdt(j) = -(vdot(j-(i-1)*Nz)/evoid/Ax)*dcPCdz(j)-(cPC(j)/Ax)*dvdotdz(j-(i-1)*Nz)+De(i)/evoid*d2cPCdz2(j)-((1-evoid)/(evoid))*rhob*dqPCdt(j);
    end
end

% Equations for qPC
for i = 1:noOfComps
    for j = (i-1)*Nz+1+noOfComps*Nz:i*Nz+noOfComps*Nz
        dqPCdt(j) = 6/d*Kc(i)*(((qs(i)*cPC(j-noOfComps*Nz))/(KL(i)+cPC(j-noOfComps*Nz)))-qPC(j));
    end
end
dcqPCdt = [dcPCdt;dqPCdt;dvdotdt];
end