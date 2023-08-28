% SPTFF Simulator
% function cTFF = SPTFFSimulator
close all; clc

% Input data
noOfComps = 3;
cTFFfeed = [10, 5, 2.5]; % kg/m^3, Feed concentrations
vdotfeed = 1e-5; % m^3/s, Feed volumetric flowrate
Patm = 101325;
L = 0.1;
Nz = 100;
z = linspace(0,L,Nz);
d = 0.05;
H = 2.5e-6;
k = 2.5e-2;
rho = 2.5e3;
mu = 2.5e-5;
De = [6.25e-4, 6.25e-6, 6.25e-6]; % Axial Dispersion coefficient

% Inlet (z = 0) concentrations (at all times)
cTFF = zeros(1,noOfComps*Nz);
for i = 1:noOfComps
    cTFF((i-1)*Nz+1) = cTFFfeed(i); % kg/m3
end

% Step for time
% Take time 80 hour
t0 = 0;
tf = 2;
% tspan = t0:tf;
tspan = [t0 tf]';

% Initial (t = 0) conditions (for all state variables)
cTFF0 = zeros(1,noOfComps*Nz);
for i = 1:noOfComps
    cTFF0((i-1)*Nz+1) = cTFF((i-1)*Nz+1); % kg/m3
end

% ODE Solver
odefun = @(t,cTFF) SPTFFModel(t,cTFF,noOfComps,cTFFfeed,vdotfeed,Patm,L,Nz,d,H,k,rho,mu,De); % odefun is my function handle!
options = [];
[t,cTFF] = ode45(odefun,tspan,cTFF0,options);

% Boundary Conditions
for i = 1:noOfComps
    cTFF(:,(i-1)*Nz+1) = cTFFfeed(i);
    cTFF(:,i*Nz) = (4*cTFF(:,i*Nz-1)-cTFF(:,i*Nz-2))/3;
    %     if i == 1
    %         cTFF(:,i*Nz) = c_factor*cTFFfeed(i);
    %     else
    %         cTFF(:,i*Nz) = (4*cTFF(:,i*Nz-1)-cTFF(:,i*Nz-2))/3;
    %     end
end

% Performance variables
% c_factor = cTFF(:,Nz)/cTFFfeed(1); % Volumetric concentration factor (VCF)
Ret = 1-cTFF(:,Nz)/cTFFfeed(1); % Retention
Sep_fac = (cTFF(:,Nz)./(cTFF(:,2*Nz)+cTFF(:,3*Nz)))/(cTFFfeed(1)/(cTFF(2)+cTFF(3))); % Separation factor

% Visual 2D
figure(1)
imagesc(z,t,cTFF(:,1:Nz))
title('cTFF_A in Filteration Stream, [Cells/m^3]')
xlabel('Position, x [m]')
ylabel('Time, [hr]')
colormap jet
colorbar
grid on

figure(2)
imagesc(z,t,cTFF(:,Nz+1:2*Nz))
title('cTFF_B in Filteration Stream, [kg/m^3]')
xlabel('Position, x [m]')
ylabel('Time, [hr]')
colormap jet
colorbar
grid on

figure(3)
imagesc(z,t,cTFF(:,2*Nz+1:3*Nz))
title('cTFF_C in Filteration Stream, [kg/m^3]')
xlabel('Position, x [m]')
ylabel('Time, [hr]')
colormap jet
colorbar
grid on

figure(4)
plot(t,cTFF(:,Nz))
title('cSPTFF_m_A_b(Retentate), [kg/m^3] vs t, [hr]')
xlabel('Time, t [hr]')
ylabel('cSPTFF_m_A_b(Retentate), [kg/m^3]')

figure(5)
plot(t,cTFF(:,2*Nz))
title('cSPTFF_S_u_b_s_t_r_a_t_e_s(Retentate), [kg/m^3] vs t, [hr]')
xlabel('Time, t [hr]')
ylabel('cSPTFF_S_u_b_t_r_a_t_e_s(Retentate), [kg/m^3]')

figure(6)
plot(t,cTFF(:,3*Nz))
title('cSPTFF_M_e_t_a_b_o_l_i_t_e_s(Retentate), [kg/m^3] vs t, [hr]')
xlabel('Time, t [hr]')
ylabel('cSPTFF_M_e_t_a_b_o_l_i_t_e_s(Retentate), [kg/m^3]')

figure(7)
plot(t,P_axi(:,Nz))
title('P_{axi}(Retentate), [N/m^2] vs t, [hr]')
xlabel('Time, t [hr]')
ylabel('P_{axi}(Retentate), [N/m^2]')

figure(8)
plot(t,vdot(:,Nz))
title('vdot, [m^3/h] vs t, [hr]')
xlabel('Time, t [hr]')
ylabel('vdot, [m^3/h]')
% end