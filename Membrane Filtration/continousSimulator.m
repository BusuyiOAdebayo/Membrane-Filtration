% Continous Simulator
close all; clear; clc
noOfComps = 10;
% Simulate Perfusion Production Bioreactor
% Initial Condition
c1initial = 1.75; % Viable cells density, [E6 Cells/mL]
% c1target = 100; % Viable cells density, [E6 Cells/mL]
c2initial = 1e-2; % Dead cells density, [E6 Cells/mL]
c3initial = c1initial+c2initial; % Total cells density, [E6 Cells/mL]
%c4initial = 4.94*1000/180.156; % Glucose conc, [mM]
c4initial = 4.94; % Glucose conc, [g/L]
%c4to7 = 0.567*c4initial;
c5initial = 0.00; %0.50*c4initial; % Glutamine conc, [mM]
c6initial = 0.25*c4initial; % Asparagine conc, [mM]
%c7initial = 0.51*1000/89.07; % Lactate conc, [mM]
c7initial = 0.51; % Lactate conc, [g/L]
c8initial = 0.96; % Ammonia conc, [mM]
c9initial = 0.00; % Titre, [g/L]
c10initial = 2.46; % Glutamate conc, [mM]
c11initial = 2500; % Culture volume, [mL]
c12initial = (c1initial/c3initial)*100; % Culture cell viability
% c13initial = 78; % Dissolved N2 conc
% c14initial = 21; % Dissolved O2 conc
% c15initial = 1; % Dissolved CO2 conc
% c16initial = 7.0; % pH
% c17initial = 250; % Osmolality, mOsmo/kg

% Variables
% c(1) = Cvc = Viable cells concentration
% c(2) = Cdc = Dead cells concentration
% c(3) = Ctc = Total cells concentration
% c(4) = CGlc = Glucose concentration
% c(5) = CGln = Glutamine concentration
% c(6) = CAsp = Asparagine concentration
% c(7) = CLac = Lactate concentration
% c(8) = CAmm = Ammonia concentration
% c(9) = CmAb = Monoclonal Antibody concentration
% c(10) = CGlu = Glutamate concentration
% c(11) = V = Reactor culture volume
% c(12) = Vi = Reactor culture cell viability
% c(13) = CDN2 = Dissolved N2 concentration
% c(14) = CDO2 = Dissolved O2 concentration
% c(15) = CDCO2 = Dissolved CO2 concentration
% c(16) = pH = Reactor culture pH
% c(17) = Osmo = Reactor culture osmolarity


tinit = 0;
tfin = 30;
% numberOfSteps = 30;
tSpan = tinit:1:tfin;%[tinit tfin]';
% tSpan = linspace(tinit,tfin);
% tSpan = linspace(tinit,tfin,numberOfSteps+1);

%%% For control purpose
% SimTime = zeros(0,1); % [];
% SimState = zeros(0,12); % [];
%
% tinit = 0;
% tfin = 30;
% numberOfSteps = tfin;
% timeV = linspace(tinit,tfin,numberOfSteps+1);
% c0 = [c1initial c2initial c3initial c4initial c5initial c6initial c7initial c8initial c9initial c10initial c11initial c12initial]; % c13initial c14initial c15initial c16initial];
%
% Kp = -0.2;
% Ti = 0.5;
% % Td = 0.5;
% % delta_t = 1;
%
% % previousError = 0; % Needed for derivative control part of PID
% previousIntegral = 0;
% previousFbleed = 0;
% previousc1 = c1initial;
%
% for j = 2:length(timeV)
%     % tSpan = [tinit,timeV(j)];
%     tSpan = linspace(tinit,timeV(j));
%     delta_t(j) = tSpan(j)-tSpan(j-1);
%     error(j) = c1target - previousc1; % c(1) error, assuming no plantmodel mismatch(i.e., simulated values = measured values)
%     proportional(j) = error(j);
%     integral(j) = previousIntegral + error(j)*delta_t(j);
%     %     derivative(j) = (error(j)-previousError)/delta_t(j);
%     output(j) = Kp*(proportional(j) + 1/Ti*previousIntegral(j)); % For PI
%     %     output(j) = Kp*(proportional(j) + 1/Ti*integral(j) + Td*derivative(j)); % For PID
%     Fbleed(j) = max(0, min(Fmedia(j), previousFbleed + output(j)));
%     [t,c] = ode45(@perfusionDesignCase1_1, tSpan, c0);
%     SimTime = [SimTime;t];
%     SimState = [SimState;c(:,1:12)];
%     tinit = timeV(j);
%     c0 = c(end,1:12);
%     %     previousError = error(j); % Needed for derivative control part of PID
%     previousIntegral = integral(j);
%     previousFbleed = Fbleed(j);
%     previousc1 = c(j,1);
% end
odefun = @(t,c) continousPerfusionBioreactorModel(t,c);
c0 = [c1initial c2initial c3initial c4initial c5initial c6initial c7initial c8initial c9initial c10initial c11initial c12initial]; % c13initial c14initial c15initial c16initial];
options = [];
[t,c] = ode45(odefun, tSpan, c0, options);
%[t,c] = ode15s(@perfusionDesignCase1_1, tSpan, c0)
%plot(t,c(:,1), '-o', t,c(:,2),'-o', t,c(:,3), '-o', t,c(:,4),'-o', t,c(:,5), '-o', t,c(:,6),'-o', t,c(:,7), '-o', t,c(:,8),'-o', t,c(:,9), '-o', t,c(:,10),'-o', t,c(:,11), '-o', t,c(:,12),'-o', t,c(:,13), '-o', t,c(:,14),'-o', t,c(:,15),'-o', t,c(:,16),'-o')

varName = ["Xv", "Xd", "Xt", "GLC", "GLN", "ASP", "LAC", "AMM", "IgG", "GLU", "Vol", "Via"];
varUnit = ["[E6 Cells/mL]", "[E6 Cells/mL]", "[E6 Cells/mL]", "[g/L]", "[mM]", "[mM]", "[g/L]", "[mM]", "[g/L]", "[mM]", "[mL]", "[%]"];

figure(1);
% % Plot Xv vs t
% subplot(4,4,1);
% plot(t,c(:,1),'ob')
% xlabel('t, [days]')
% ylabel('Xv, [E6 Cells/mL]')
% legend('Xv vs t');
% % Plot Xd vs t
% subplot(4,4,2);
% plot(t,c(:,2),'og')
% xlabel('t, [days]')
% ylabel('Xd, [E6 Cells/mL]')
% legend('Xd vs t');
% % Plot Xt vs t
% subplot(4,4,3);
% plot(t,c(:,3),'or')
% xlabel('t, [days]')
% ylabel('Xt, [E6 Cells/mL]')
% legend('Xt vs t');
% % Plot [GLC] vs t
% subplot(4,4,4);
% plot(t,c(:,4),'ob')
% xlabel('t, [days]')
% ylabel('[GLC], [g/L]')
% legend('[GLC] vs t');
% % Plot [GLN] vs t
% subplot(4,4,5);
% plot(t,c(:,5),'og')
% xlabel('t, [days]')
% ylabel('[GLN], [mM]')
% legend('[GLN] vs t');
% % Plot [ASP] vs t
% subplot(4,4,6);
% plot(t,c(:,6),'or')
% xlabel('t, [days]')
% ylabel('[ASP], [mM]')
% legend('[ASP] vs t');
% % Plot [LAC] vs t
% subplot(4,4,7);
% plot(t,c(:,7),'ob')
% xlabel('t, [days]')
% ylabel('[LAC], [g/L]')
% legend('[LAC] vs t');
% % Plot [AMM] vs t
% subplot(4,4,8);
% plot(t,c(:,8),'og')
% xlabel('t, [days]')
% ylabel('[AMM], [mM]')
% legend('[AMM] vs t');
% % Plot Titer vs t
% subplot(4,4,9);
% plot(t,c(:,9),'or')
% xlabel('t, [days]')
% ylabel('Titer, [g/L]')
% legend('Titer vs t');
% % Plot [GLU] vs t
% subplot(4,4,10);
% plot(t,c(:,10),'og')
% xlabel('t, [days]')
% ylabel('[GLU], [mM]')
% legend('[GLU] vs t');
% % Plot V vs t
% subplot(4,4,11);
% plot(t,c(:,11),'or')
% xlabel('t, [days]')
% ylabel('V, [mL]')
% legend('V vs t');
% % Plot Cell Viability vs t
% subplot(4,4,12);
% plot(t,c(:,12),'og')
% %plot(t,c(:,(c(1)/c(3))*100),'og') % Code ran but could not plot this. Error msg: Index in position 2 is invalid. Array indices must be positive integers or logical values.
% xlabel('t, [days]')
% ylabel('Cell Viability, [%]')
% legend('Cell Viability vs t');
for i = 1:length(varName)
    subplot(4,4,i);
    plot(t,c(:,i),'ob')
    legend(strcat(varName(i),"(Retentate), ",varUnit(i)," vs t, [hr]"))
    xlabel('Time, t, [days]')
    ylabel(strcat(varName(i),"(Retentate), ",varUnit(i)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate SPTFF
% Input data
% noOfComps = 10;
cTFFfeed = zeros(1,noOfComps);
%for ii = 2:length(t)
for i = 1:noOfComps
    if 1 <= i && i <= 3
        cTFFfeed(i) = 0; % kg/m^3, Feed concentrations
    else
        cTFFfeed(i) = c(end,i); % kg/m^3, Feed concentrations
    end
end
% cTFFfeed = c(end,1:10); % kg/m^3, Feed concentrations
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
De = [6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-4, 6.25e-6]; % Axial Dispersion coefficient

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
odefun = @(t,cTFF) continousSPTFFModel(t,cTFF,noOfComps,cTFFfeed,vdotfeed,Patm,L,Nz,d,H,k,rho,mu,De); % odefun is my function handle!
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
% c_factor = cTFF(:,9*Nz)/cTFFfeed(9); % Volumetric concentration factor (VCF)
% Ret_factor = 1-cTFF(:,9*Nz)/cTFFfeed(9); % Retention
% Sep_factor = (cTFF(:,9*Nz)./(cTFF(:,2*Nz)+cTFF(:,3*Nz)))/(cTFFfeed(9)/(cTFF(2)+cTFF(3))); % Separation factor

varName = ["Xv", "Xd", "Xt", "GLC", "GLN", "ASP", "LAC", "AMM", "IgG", "GLU"];
varUnit = ["[E6 Cells/mL]", "[E6 Cells/mL]", "[E6 Cells/mL]", "[g/L]", "[mM]", "[mM]", "[g/L]", "[mM]", "[g/L]", "[mM]"];

% Visual 2D
for i = 1:noOfComps
    figure(1+i)
    imagesc(t,z,cTFF(:,(i-1)*Nz+1:i*Nz))
%     title('cTFF_A in Filteration Stream, [Cells/m^3]')
%     xlabel('Position, x [m]')
%     ylabel('Time, [hr]')
    title(strcat(varName(i),"(Retentate), ",varUnit(i)," vs t, [hr] and z, [m]"))
    xlabel('Time, t [hr]')
    ylabel('Space, z [m]')
    colormap jet
    colorbar
    grid on
end
% figure(2)
% imagesc(z,t,cTFF(:,1:Nz))
% title('cTFF_A in Filteration Stream, [Cells/m^3]')
% xlabel('Position, x [m]')
% ylabel('Time, [hr]')
% colormap jet
% colorbar
% grid on
%
% figure(3)
% imagesc(z,t,cTFF(:,Nz+1:2*Nz))
% title('cTFF_B in Filteration Stream, [kg/m^3]')
% xlabel('Position, x [m]')
% ylabel('Time, [hr]')
% colormap jet
% colorbar
% grid on
%
% figure(4)
% imagesc(z,t,cTFF(:,2*Nz+1:3*Nz))
% title('cTFF_C in Filteration Stream, [kg/m^3]')
% xlabel('Position, x [m]')
% ylabel('Time, [hr]')
% colormap jet
% colorbar
% grid on

% Visual 1D
for i = 1:noOfComps
    figure(1+noOfComps+i)
    plot(t,cTFF(:,i*Nz))
%     title('cSPTFF_m_A_b(Retentate), [kg/m^3] vs t, [hr]')
%     xlabel('Time, t [hr]')
%     ylabel('cSPTFF_m_A_b(Retentate), [kg/m^3]')
    title(strcat(varName(i),"(Retentate), ",varUnit(i)," vs t, [hr]"))
    xlabel('Time, t [hr]')
    ylabel(strcat(varName(i),"(Retentate), ",varUnit(i)))
end
% figure(5)
% plot(t,cTFF(:,Nz))
% title('cSPTFF_m_A_b(Retentate), [kg/m^3] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('cSPTFF_m_A_b(Retentate), [kg/m^3]')
%
% figure(6)
% plot(t,cTFF(:,2*Nz))
% title('cSPTFF_S_u_b_s_t_r_a_t_e_s(Retentate), [kg/m^3] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('cSPTFF_S_u_b_t_r_a_t_e_s(Retentate), [kg/m^3]')
%
% figure(7)
% plot(t,cTFF(:,3*Nz))
% title('cSPTFF_M_e_t_a_b_o_l_i_t_e_s(Retentate), [kg/m^3] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('cSPTFF_M_e_t_a_b_o_l_i_t_e_s(Retentate), [kg/m^3]')

% figure(8)
% plot(t,P_axi(:,Nz))
% title('P_{axi}(Retentate), [N/m^2] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('P_{axi}(Retentate), [N/m^2]')
%
% figure(9)
% plot(t,vdot(:,Nz))
% title('vdot, [m^3/h] vs t, [hr]')
% xlabel('Time, t [hr]')
% ylabel('vdot, [m^3/h]')
% % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate PC
% Input data
%noOfComps = 3;
noOfGroupVariables = 2;
%cPCfeed = [3.75, 1.5, 1.5]; % kg/m^3, Feed concentrations
%cPCfeed = zeros(1,noOfComps);
%cPCfeed(i) = c(end,i); % kg/m^3, Feed concentrations
cPCfeed = cTFF(end,:);
vdotPCfeed = 1e-5; % m^3/s, Feed volumetric flowrate
Patm = 101325;
L = 0.1;
Nz = 50;
z = linspace(0,L,Nz);
db = 0.05;
dp = 1e-5;
evoid = 0.5;
qs = [2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 0, 2.5e-2];
kc = [2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6, 2.5e-6, 2.5e-2, 2.5e-6];
KL = [2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-2, 2.5e-6, 2.5e-2];
rhob = 2.5e3;
mu = 2.5e-5;
De = [6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-6, 6.25e-4, 6.25e-6]; % Axial Dispersion coefficient

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

% varName = ["IgG", "Substrates", "Metabolites"];%, "GLC", "GLN", "ASP", "LAC", "AMM", "IgG", "GLU"];
% varUnit = ["[g/L]", "[mM]", "[mM]"];%, "[g/L]", "[mM]", "[mM]", "[g/L]", "[mM]", "[g/L]", "[mM]"];

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
    figure(1+2*noOfComps+i)
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
    figure(1+3*noOfComps+i)
    plot(t,cqvdotPC(:,i*Nz))
    title(strcat(varName(i),"(Outlet), ",varUnit(i)," vs t, [hr]"))
    xlabel('Time, t [hr]')
    ylabel(strcat(varName(i),"(Outlet), ",varUnit(i)))
end
%end