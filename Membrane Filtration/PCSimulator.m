% PC Simulator
% function cqPC = PCSimulator
close all; clear; clc

% Input data
noOfComps = 3;
noOfGroupVariables = 2;
cPCfeed = [3.75, 1.5, 1.5]; % kg/m^3, Feed concentrations
vdotPCfeed = 1e-5; % m^3/s, Feed volumetric flowrate
Patm = 101325;
L = 0.1;
Nz = 50;
z = linspace(0,L,Nz);
db = 0.05;
dp = 1e-5;
evoid = 0.5;
qs = [2.5e-15, 2.5e-2, 2.5e-2];
kc = [2.5e-6, 2.5e-2, 2.5e-2];
KL = [2.5e-2, 2.5e-6, 2.5e-6];
rhob = 2.5e3;
mu = 2.5e-5;
De = [6.25e-4, 6.25e-6, 6.25e-6]; % Axial Dispersion coefficient

% Inlet (z = 0) values of spatial differential variables (at all times)
cqvdotPC = zeros(1,noOfGroupVariables*noOfComps*Nz+Nz);
for i = 1:noOfComps
    cqvdotPC((i-1)*Nz+1) = cPCfeed(i); % kg/m3
end
cqvdotPC(noOfGroupVariables*noOfComps*Nz+1) = vdotPCfeed; % kg/m3

% Step for time
% Take time 80 hour
t0 = 0;
tf = 10;
% tspan = t0:tf;
tspan = [t0 tf]';

% Initial (t = 0) conditions of time differential variables, i.e., state variables (at all locations)
cqvdotPC0 = zeros(1,noOfGroupVariables*noOfComps*Nz+Nz);
for i = 1:noOfComps
    cqvdotPC0((i-1)*Nz+1) = cqvdotPC((i-1)*Nz+1); % kg/m3
end
cqvdotPC0(noOfGroupVariables*noOfComps*Nz+1) = cqvdotPC(noOfGroupVariables*noOfComps*Nz+1);

% ODE Solver
odefun = @(t,cqvdotPC) continousPCModel(t,cqvdotPC,noOfComps,cPCfeed,vdotPCfeed,Patm,L,Nz,db,dp,evoid,qs,kc,KL,rhob,mu,De); % odefun is my function handle!
options = [];
[t,cqvdotPC] = ode45(odefun,tspan,cqvdotPC0,options);

varName = ["IgG", "Substrates", "Metabolites"];%, "GLC", "GLN", "ASP", "LAC", "AMM", "IgG", "GLU"];
varUnit = ["[g/L]", "[mM]", "[mM]"];%, "[g/L]", "[mM]", "[mM]", "[g/L]", "[mM]", "[g/L]", "[mM]"];

% Boundary Conditions
for i = 1:noOfComps
    cqvdotPC(:,(i-1)*Nz+1) = cPCfeed(i);
    cqvdotPC(:,i*Nz) = (4*cqvdotPC(:,i*Nz-1)-cqvdotPC(:,i*Nz-2))/3;
    %     if i == 1
    %         cTFF(:,i*Nz) = c_factor*cTFFfeed(i);
    %     else
    %         cTFF(:,i*Nz) = (4*cTFF(:,i*Nz-1)-cTFF(:,i*Nz-2))/3;
    %     end
end

% Performance variables
% c_factor = cTFF(:,Nz)/cTFFfeed(1); % Volumetric concentration factor (VCF)
% Ret = 1-cqPC(:,Nz)/cTFFfeed(1); % Retention
% Sep_fac = (cqPC(:,Nz)./(cqPC(:,2*Nz)+cqPC(:,3*Nz)))/(cTFFfeed(1)/(cqPC(2)+cqPC(3))); % Separation factor

% Visual 2D
for i = 1:noOfComps
    figure(i)
    imagesc(t,z,cqvdotPC(:,(i-1)*Nz+1:i*Nz))
    title(strcat(varName(i),"(Outlet), ",varUnit(i)," vs t, [hr] and z, [m]"))
%     xlabel('Time, t [hr]')
%     ylabel(strcat(varName(i),"(Outlet), ",varUnit(i)))
    xlabel('Time, t [hr]')
    ylabel('Space, z [m]')
    colormap jet
    colorbar
    grid on
end

% Visual 1D
for i = 1:noOfComps
    figure(noOfComps+i)
    plot(t,cqvdotPC(:,i*Nz))
    title(strcat(varName(i),"(Outlet), ",varUnit(i)," vs t, [hr]"))
    xlabel('Time, t [hr]')
    ylabel(strcat(varName(i),"(Outlet), ",varUnit(i)))
end

% figure(2*noOfComps+1)
% plot(t,P_axi(:,Nz))
% title('P_{axi}(Outlet), [N/m^2] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('P_{axi}(Outlet), [N/m^2]')
%
% figure(2*noOfComps+2)
% plot(t,vdot(:,Nz))
% title('vdot, [m^3/h] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('vdot, [m^3/h]')
% end
% end