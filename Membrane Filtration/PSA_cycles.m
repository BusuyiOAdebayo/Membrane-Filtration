% Clear previous files
% clear all
close all
clear
clc
format short
%
% parameters shared with the ODE routine
%
global Mg1 Mg2 De1 De2 taup adsorptionORpurge dYbardZbarDifferentialMethod rhog0 T0 ncall n zbar dzbar tbar epsilon epsilonbar tc tt tstep k1 k2 rhob n10 n20 dp L D A ta td tp tdep y1a y1d y1feed y10 Y0 qa Qa Qd P0 Pa Pd Ta Pmax Pmin Tmax Tmin Tref Q0 MassAxial tbar tbara tbarp tbardep tbard nout nouta noutp noutdep noutd noutcycle noutstep cycle miu dtbara dtbarp dtbardep dtbard
%
% Input data PSA_inputs
%
% ODE integration
%
ncall = 0;
reltol = 1.0e-05;
abstol = 1.0e-05;
options = odeset('RelTol',reltol,'AbsTol',abstol);
%
% Cycle #1 pressurization
%
% Initial conditions for y1bar, n1bar, n2bar for i = 1:n
%
for i = 1:n
    Y0(0*n+i) = y10*rhog0*L/(Qa*tt); % for Y1bar
    Y0(1*n+i) = (k1*rhob*L/Qa)*n10; % for n1bar
    Y0(2*n+i) = (k2*rhob*L/Qa)*n20; % for n2bar
    Y0(3*n+i) = P0/(Pmax-Pmin); % for Pbar
end
for i = 1
    Y0(3*n+1) = Pa/(Pmax-Pmin);
end
%
[~,Ys] = ode15s(@PSA_pressurization, tbarp, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_pressurization, tbarp, Y0, options) to [~,Ys] = ode15s(@PSA_pressurization, tbarp, Y0, options);
%
Y = zeros(noutp,4*n); % Added this!
for it = 1:noutp % time
    for i = 1:4*n % position
        Y(it,i) = Ys(it,i);
    end
end
%
% Cycle #1 adsorption
%
% Initial conditions for adsorption step
%
for i = 1:n
    Y0(0*n+i) = Y(noutp,0*n+i); % for Y1bar
    Y0(1*n+i) = Y(noutp,1*n+i); % for n1bar
    Y0(2*n+i) = Y(noutp,2*n+i); % for n2bar
    Y0(3*n+i) = Pa/(Pmax-Pmin); % for Pbar
end
%
[~,Ys] = ode15s(@PSA_adsorption, tbara, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_adsorption, tbara, Y0, options) to [t,Ys] = ode15s(@PSA_adsorption, tbara, Y0, options);
%
for it = 1:(nouta-1) % time
    for i = 1:4*n % position
        Y(it+noutp,i) = Ys(it+1,i);
    end
end
%
% Cycle #1 depressurization
%
% Initial conditions for depressurization step
%
for i = 1:n
    Y0(0*n+i) = Y(noutp+nouta-1,0*n+i); % for Y1bar
    Y0(1*n+i) = Y(noutp+nouta-1,1*n+i); % for n1bar
    Y0(2*n+i) = Y(noutp+nouta-1,2*n+i); % for n2bar
    Y0(3*n+i) = Y(noutp+nouta-1,3*n+i); % for Pbar
end
Y0(3*n+n) = Pd/(Pmax-Pmin);
%
[~,Ys] = ode15s(@PSA_depressurization, tbardep, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_depressurization, tbardep, Y0, options) to [~,Ys] = ode15s(@PSA_depressurization, tbardep, Y0, options);
%
for it = 1:(noutdep-1) % time
    for i = 1:4*n % position
        Y(it+noutp+nouta-1,i) = Ys(it+1,i);
    end
end
%
% Cycle #1 desorption
%
% Initial conditions for desorption step
%
for i = 1:n
    Y0(0*n+i) = Y(noutp+nouta+noutdep-2,0*n+i); % for Y1bar
    Y0(1*n+i) = Y(noutp+nouta+noutdep-2,1*n+i); % for n1bar
    Y0(2*n+i) = Y(noutp+nouta+noutdep-2,2*n+i); % for n2bar
    Y0(3*n+i) = Pd/(Pmax-Pmin); % for Pbar
end
%
[~,Ys] = ode15s(@PSA_desorption, tbard, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_desorption, tbard, Y0, options) to [~,Ys] = ode15s(@PSA_desorption, tbard, Y0, options);
%
for it = 1:(noutd-1) % time
    for i = 1:4*n
        Y(it+nouta+noutp+noutdep-2,i) = Ys(it+1,i);
    end
end
%
% For cycles after the first cycle
%
for ic = 1:(cycle-1)
    %
    % pressurization
    %
    % Initial conditions for pressurization step
    for i = 1:n
        Y0(0*n+i) = Y((noutcycle-1)*ic+1,0*n+i); % for Y1bar
        Y0(1*n+i) = Y((noutcycle-1)*ic+1,1*n+i); % for n1bar
        Y0(2*n+i) = Y((noutcycle-1)*ic+1,2*n+i); % for n2bar
        Y0(3*n+i) = Y((noutcycle-1)*ic+1,3*n+i); % for Pbar
    end
    Y0(3*n+1) = Pa/(Pmax-Pmin);
    %
    [~,Ys] = ode15s(@PSA_pressurization, tbarp, Y0, options); % Changed from [~,Ys] = ode15s(@PSA_pressurization, tbard, Y0, options) to [~,Ys] = ode15s(@PSA_pressurization, tbard, Y0, options);
    %
    for it = 1:(noutp-1) % time
        for i = 1:4*n
            Y(it+ic*(noutcycle-1)+1,i) = Ys(it+1,i);
        end
    end
    %
    % Adsorption
    %
    % Initial conditions for adsorption step
    %
    for i = 1:n
        Y0(0*n+i) = Y((noutcycle-1)*ic+noutp,0*n+i); % for Y1bar
        Y0(1*n+i) = Y((noutcycle-1)*ic+noutp,1*n+i); % for n1bar
        Y0(2*n+i) = Y((noutcycle-1)*ic+noutp,2*n+i); % for n2bar
        Y0(3*n+i) = Pa/(Pmax-Pmin); % for Pbar
    end
    [~,Ys] = ode15s(@PSA_adsorption, tbara, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_adsorption, tbard, Y0, options) to [~,Ys] = ode15s(@PSA_adsorption, tbard, Y0, options);
    for it = 1:(nouta-1) % time
        for i = 1:4*n
            Y(it+ic*(noutcycle-1)+noutp,i) = Ys(it+1,i);
        end
    end
    %
    % Depressurization
    %
    % Initial conditions for depressurization step
    %
    for i = 1:n
        Y0(0*n+i) = Y((noutcycle-1)*ic+noutp+nouta-1,0*n+i); % for Y1bar
        Y0(1*n+i) = Y((noutcycle-1)*ic+noutp+nouta-1,1*n+i); % for n1bar
        Y0(2*n+i) = Y((noutcycle-1)*ic+noutp+nouta-1,2*n+i); % for n2bar
        Y0(3*n+i) = Y((noutcycle-1)*ic+noutp+nouta-1,3*n+i); % for Pbar
    end
    Y0(3*n+n) = Pd/(Pmax-Pmin);
    %
    %
    [~,Ys] = ode15s(@PSA_depressurization, tbardep, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_depressurization, tbard, Y0, options) to [~,Ys] = ode15s(@PSA_depressurization, tbard, Y0, options);
    for it = 1:(noutdep-1) % time
        for i = 1:4*n
            Y(it+ic*(noutcycle-1)+nouta+noutp-1,i) = Ys(it+1,i);
        end
    end
    %
    % Desorption
    %
    % Initial conditions for desorption step
    %
    for i = 1:n
        Y0(0*n+i) = Y((noutcycle-1)*ic+nouta+noutp+noutdep-2,0*n+i); % for Y1bar
        Y0(1*n+i) = Y((noutcycle-1)*ic+nouta+noutp+noutdep-2,1*n+i); % for n1bar
        Y0(2*n+i) = Y((noutcycle-1)*ic+nouta+noutp+noutdep-2,2*n+i); % for n2bar
        Y0(3*n+i) = Pd/(Pmax-Pmin); % for Pbar
    end
    [~,Ys] = ode15s(@PSA_desorption, tbard, Y0, options); % Changed from [t,Ys] = ode15s(@PSA_desorption, tbard, Y0, options) to [~,Ys] = ode15s(@PSA_desorption, tbard, Y0, options);
    for it = 1:(noutd-1) % time
        for i = 1:4*n
            Y(it+ic*(noutcycle-1)+nouta+noutp+noutdep-2,i) = Ys(it+1,i);
        end
    end
    %
end
%
% One vector to three vectors
%
Y1bar = zeros(nout,n); % Added these six!
n1bar = zeros(nout,n);
n2bar = zeros(nout,n);
Pbar = zeros(nout,n);
y1bar = zeros(nout,n);
Mg = zeros(nout,n);
for it = 1:nout % time
    for i = 1:n % position
        Y1bar(it,i) = Y(it,0*n+i);
        n1bar(it,i) = Y(it,1*n+i);
        n2bar(it,i) = Y(it,2*n+i);
        Pbar(it,i) = Y(it,3*n+i);
        y1bar(it,i) = Y1bar(it,i)/(rhog0*L/(Qa*tt));
        Mg(it,i) = Mg1*y1bar(it,i)+Mg2*(1-y1bar(it,i));
    end
end
%
% Fix the boundary conditions
%
% ODE solver does take care of initial conditions at it = 1
%
for ic = 1:cycle
    for it = 2:(nouta+noutp-1) % time
        i = 1;
        y1bar(it+(noutcycle-1)*(ic-1),i) = y1a;
    end
    %
    for it = 2:(noutdep+noutd-1)
        i = n;
        y1bar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i) = y1d;
    end
end
%
% Isobaric, Isothermal
%
thetabar = zeros(nout,n); % Added these two!
rhogbar = zeros(nout,n);
for it = 1:nout
    for i = 1:n
        thetabar(it,i) = (T0-Tref)/(Tmax-Tmin);
        rhogbar(it,i) = rhog0*L/(Qa*tt);
    end
end
%
% Formulate nleqmbar(i), n2eqmbar(i), Qbar(i), usbar(i) with inital and boundary conditions
%
n1eqmbar = zeros(nout,n); % Added these two!
n2eqmbar = zeros(nout,n);
for it = 1:nout
    for i = 1:n
        py = (y1bar(it,i))*((Pmax-Pmin)*Pbar(it,i)/1013250);
        if py < 0.00013
            n1eqmbar(it,i) = 0;
        else
            n1eqmbar(it,i) = -(k1*rhob*L/Qa)*(3*0.069337*((y1bar(it,i))*((Pmax-Pmin)*Pbar(it,i)/1013250))^(-0.10825)/(1+0.069337*((y1bar(it,i))*((Pmax-Pmin)*Pbar(it,i)/1013250))^(-0.10825))-0.4633)/2/193.842*1000;
        end
        n2eqmbar(it,i) = 0;
    end
end
%
Qbar = zeros(nout,n); % Added these two!
usbar = zeros(nout,n);
for it = 1:nout % Changed from it = 1 to it = 1:nout
    for i = 1:n
        Qbar(it,i) = Qa/Qa;
        usbar(it,i) = Qbar(it,i)/rhogbar(it,i);
    end
end
%
for ic = 1:cycle % cycle
    for it = 2:(nouta+noutp-1) % time
        for i = 1
            Qbar(it+(noutcycle-1)*(ic-1),i) = Qa/Qa;
            usbar(it+(noutcycle-1)*(ic-1),i) = Qbar(it+(noutcycle-1)*(ic-1),i)/rhogbar(it+(noutcycle-1)*(ic-1),i);
        end
        for i = 2:n
            Qbar(it+(noutcycle-1)*(ic-1),i) = -((n1eqmbar(it+(noutcycle-1)*(ic-1),i)-n1bar(it+(noutcycle-1)*(ic-1),i))+(n2eqmbar(it+(noutcycle-1)*(ic-1),i)-n2bar(it+(noutcycle-1)*(ic-1),i)))*dzbar+Qbar(it+(noutcycle-1)*(ic-1),i-1);
            usbar(it+(noutcycle-1)*(ic-1),i) = Qbar(it+(noutcycle-1)*(ic-1),i)/rhogbar(it+(noutcycle-1)*(ic-1),i);
        end
    end
    for it = 2:(noutdep+noutd-1)
        for i = n
            Qbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i) = -Qd/Qa;
            usbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i) = Qbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i)/rhogbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i);
        end
        for i = n-1:-1:1
            Qbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i) = +((n1eqmbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i)-n1bar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i))+(n2eqmbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i)-n2bar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i)))*dzbar+Qbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i+1);
            usbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i) = Qbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i)/rhogbar(it+(noutcycle-1)*(ic-1)+nouta+noutp-2,i);
        end
    end
end
%
% Calculate O2 purity
%
O2 = 0;
Product = 0;
yy1bar = zeros(1,((td+tdep)/tstep)); % Added these three!
QQdbar = zeros(1,((td+tdep)/tstep));
c1bar = zeros(1,((td+tdep)/tstep));
for ic = 1:cycle
    for it = 1:((td+tdep)/tstep)
        yy1bar(it) = y1bar(it+(tc*(ic-1)+ta+tp)/tstep+1,1);
        QQdbar(it) = Qbar(it+(tc*(ic-1)+ta+tp)/tstep+1,1);
        c1bar(it) = yy1bar(it)*QQdbar(it);
    end
    t = tstep:tstep:((td+tdep)/tstep);
    O2 = O2 + trapz(t,c1bar);
    Product = Product + trapz(t,QQdbar);
end
Pu = O2/Product;
%
% Calculate O2 recovery
%
QQabar = zeros(1,((ta+tp)/tstep)); % Added this!
O2feed = 0;
for ic = 1:cycle
    for it = 1:((ta+tp)/tstep)
        QQabar(it) = Qbar(it+tc/tstep*(ic-1)+1,1);
    end
    t = tstep:tstep:((ta+tp)/tstep);
    O2feed = O2feed + trapz(t,QQabar)*y1a;
end

% Suggested this code block
% QQabar = zeros(cycle,((ta+tp)/tstep));
% O2feed = zeros(cycle,1);
% t = tstep:tstep:((ta+tp)/tstep);
% for ic = 1:cycle
%     for it = 1:((ta+tp)/tstep)
%         QQabar(ic,it) = Qbar(it+tc/tstep*(ic-1)+1,1);
%     end
%     %     t = tstep:tstep:((ta+tp)/tstep);
%     O2feed(ic) = O2feed(ic) + trapz(t,QQabar(ic,:))*y1a; % Changed from O2feed = O2feed + trapz(t,QQabar)*y1a to O2feed(i) = O2feed(i) + trapz(t,QQabar(i,:))*y1a;
% end

Re = -O2/O2feed;
% Calculate the amount of adsorbent (g)
m = A*L*rhob;
%
% Calculate O2 Productivity (mmol/sec/kg)
%
Pro = qa/m*1000*Re;
%
% Display selected output
%
fprintf('\n abstol = %8.1e reltol = %8.1e\n', abstol, reltol);
fprintf('\n ncall = %4d\n', ncall);
%
% Display important parameters
%
fprintf('\n process step = %d dYbardZbar differential method = %d',adsorptionORpurge, dYbardZbarDifferentialMethod);
fprintf('\n P0 = %g atm T0 = %g K yN20 = %g Q0 = %gmmol/cm2/sec', P0/1013250, T0, (1-y10), Q0);
fprintf('\n ta = %g sec ta = %g K yN2feed = %g Qa = %gmmol/cm2/sec Qa = %g mmol/sec', ta, ta, (1-y1a), Qa, Qa);
fprintf('\n dp = %g um L = %g cm D = %g cm', dp*10000, L, D);
fprintf('\n k1 = %g sec-1 k2 = %g sec-1 rhob = %g g/ml', k1, k2,rhob);
fprintf('\n Adsorbent mass = %g g', m);
fprintf('\n O2 purity = %g Oxygen Recovery = %g Oxygen productivity = %g mmol/sec/kg\n\n', Pu, Re, Pro);
%
% Plot numerical solutions vs independent variables (t,z)
%
z0 = 0.0;
zend = L;
z = linspace(z0,zend,n);
t0 = 0.0;
tend = tt;
t = linspace(t0,tend,nout);
%
choose = 2; % Added this!
if choose == 1 % start of choose == 1 : 2D subplots of y1, y2, n1, n2, T, P
    %
    % Subplots of y1, y2, n1, n2, T, P
    %
    %figure(1); % Comment out this!
    subplot(2,2,1)
    plot(t,y1bar(:,(n+1)/2),'-r',t,(1-y1bar(:,(n+1)/2)),'-b'); % axis tight
    % h = legend('y_{O2}','y_{N2}',1);
    % Replaced with below!
    h1 = legend('y_{IgG}','y_{Non-IgG}','Location','SW'); % Changed from h = legend('y_{O2}','y_{N2}','Loc','W'); to h1 = legend('y_{O2}','y_{N2}','Loc','SW');
    title('y_{IgG} and y_{Non-IgG} vs. time at L/2');
    xlabel('time (sec)');
    ylabel('y_{IgG} and y_{Non-IgG}')
    subplot(2,2,2)
    plot(t,(Qa/(k1*rhob*L)).*n1bar(:,(n+1)/2),'-r',t,(Qa/(k2*rhob*L).*n2bar(:,(n+1)/2)),'-b'); % axis tight
    % h = legend('n_{O2}','n_{N2}',1);
    % Replaced with below!
    h2 = legend('n_{IgG}','n_{Non-IgG}','Location','SW'); % Changed from h = legend('y_{O2}','y_{N2}','Loc','W'); to h2 = legend('y_{O2}','y_{N2}','Loc','W');
    title('n_{IgG} and n_{Non-IgG} vs. time at L/2');
    xlabel('time (sec)');
    ylabel('n_{IgG} and n_{Non-IgG} (mmol/g adsorbent')
    subplot(2,2,3)
    plot(t,Tref+(Tmax-Tmin).*thetabar(:,(n+1)/2)); % axis tight
    title('T vs. time at L/2');
    xlabel('time (sec)');
    ylabel('T (K)')
    subplot(2,2,4)
    plot(z,((Pmax-Pmin)/P0).*Pbar(nout,:)); % axis tight
    title('P vs. bed axial position at tfinal');
    xlabel('Bed axial position, z (cm)');
    ylabel('P (atm)')
end % end of choose == 1
%
if choose == 2 % start of choose == 2 : Individual 2D plots
    %
    figure(1)
    plot(t,y1bar(:,1),'-r>',t,(1-y1bar(:,1)),'-bs','LineWidth',2,'MarkerSize',3)
    legend('y_{IgG}','y_{Non-IgG}');
    xlabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    ylabel('Liquid phase mole fraction, y_{IgG} and y_{Non-IgG}','color','k','fontsize',10,'fontweight','b')
    title('Liquid phase mole fraction vs. time at the first cell','color','k','fontsize',12,'fontweight','b');
    axis([0 tt 0 1])
    grid on
    %
    figure(2)
    plot(t,y1bar(:,n),'-r>',t,(1-y1bar(:,n)),'-bs','LineWidth',2,'MarkerSize',3)
    legend('y_{IgG}','y_{Non-IgG}');
    xlabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    ylabel('Liquid phase mole fraction, y_{IgG} and y_{Non-IgG}','color','k','fontsize',10,'fontweight','b')
    title('Liquid phase mole fraction vs. time at the last cell','color','k','fontsize',12,'fontweight','b');
    axis([0 tt 0 1])
    grid on
    %
    figure(3)
    y1outbar = zeros((ta+tp)/tstep,1); % Added this! Maybe this: y1outbar = zeros(cycle,(ta+tp)/tstep);
    y1outbar(1) = y1bar(1,n);
    for ic = 1:cycle
        %         y1outbar = zeros((ta+tp)/tstep,1);
        for it = 1:((ta+tp)/tstep)
            y1outbar(it+(tc*(ic-1))/tstep+1) = y1bar(it+(tc*(ic-1))/tstep+1,n);
        end
        %         y1outbar = zeros((td+tdep)/tstep+(tc*(cycle-1)+ta+tp)/tstep+1);
        for it = 1:((td+tdep)/tstep)
            y1outbar(it+(tc*(ic-1)+ta+tp)/tstep+1) = y1bar(it+(tc*(ic-1)+ta+tp)/tstep+1,1);
        end
    end
    plot(t,y1outbar,'-r>','LineWidth',2,'MarkerSize',3)
    xlabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    ylabel('IgG mole fraction, y_{IgG} and y_{Non-IgG}','color','k','fontsize',10,'fontweight','b')
    title('IgG mole fraction vs. time at the outlet','color','k','fontsize',12,'fontweight','b');
    axis([0 tt 0 1])
    grid on
    %
    % Amount adsorbed on adsorbent, n1 and n2, vs. time at L/2
    %
    figure(4)
    plot(t,(Qa/(k1*rhob*L)).*n1bar(:,n),'-r>',t,(Qa/(k2*rhob*L)).*n2bar(:,n),'-bs','LineWidth',2,'MarkerSize',3);
    legend('n_{IgG}','n_{Non-IgG}');
    title('Amount adsorbed on adsorbent vs. time at the last cell','color','k','fontsize',12,'fontweight','b');
    xlabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b');
    ylabel('Amount adsorbed on adsorbent, n_{IgG} and n_{Non-IgG} (mmol IgG/g adsorbent)','color','k','fontsize',10,'fontweight','b')
    grid on
    %
    % Temperature, T, vs. time at L/2
    %
    figure(5)
    plot(t,(Tref+(Tmax-Tmin).*thetabar(:,(n+1)/2)),'LineWidth',2)
    xlabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    ylabel('Temperature, T (K)','color','k','fontsize',10,'fontweight','b')
    title('Temperature vs. time','color','k','fontsize',12,'fontweight','b');
    grid on
    %
    % Pressure, P, vs. time
    %
    figure(6)
    plot(t,((Pmax-Pmin)/P0).*Pbar(:,(n+1)/2),'LineWidth',2)
    xlabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    ylabel('Pressure, P (atm)','color','k','fontsize',10,'fontweight','b')
    title('Pressure vs. time','color','k','fontsize',12,'fontweight','b');
    grid on
    %
end % end of choose == 2
%
if choose == 3 % start of choose == 3 : Surface plots of y1, y2, n1, n2, T, P, rhog, Q, us
    %
    % Gas phase mole fraction, y1
    %
    figure(1)
    surf(z,t,y1bar,'edgecolor','none')
    colormap jet
    xlabel('Bed axial position, z (cm)','color','k','fontsize',10,'fontweight','b')
    ylabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    zlabel('Gas phase mole fraction,y_{IgG}','color','k','fontsize',10,'fontweight','b')
    axis tight
    %
    % Gas phase mole fraction, y2
    %
    figure(2)
    surf(z,t,(1-y1bar),'edgecolor','none')
    colormap jet
    xlabel('Bed axial position, z (cm)','color','k','fontsize',10,'fontweight','b')
    ylabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    zlabel('Gas phase mole fraction,y_{Non-IgG}','color','k','fontsize',10,'fontweight','b')
    axis tight
    %
    % Amount adsorbed on adsorbent, n1
    %
    figure(3)
    surf(z,t,(Qa/(k1*rhob*L)).*n1bar,'edgecolor','none')
    colormap jet
    xlabel('Bed axial position, z(cm)','color','k','fontsize',10,'fontweight','b')
    ylabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    zlabel('Amount adsorbed on adsorbent, n_{Non-IgG} (mmol Non-IgG/g adsorbent)','color','k','fontsize',10,'fontweight','b')
    axis tight
    %
    % Temperature, T
    %
    figure(4)
    surf(z,t,(Tmin+(Tmax-Tmin).*thetabar),'edgecolor','none')
    colormap jet
    xlabel('Bed axial position, z (cm)','color','k','fontsize',10,'fontweight','b')
    ylabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    zlabel('Temperature, T (K)','color','k','fontsize',10,'fontweight','b')
    axis tight
    %
    % Pressure, P
    %
    figure(5)
    surf(z,t,((Pmax-Pmin)/P0).*Pbar,'edgecolor','none')
    colormap jet
    xlabel('Bed axial position, z(cm)','color','k','fontsize',10,'fontweight','b')
    ylabel('Time, t (sec)','color','k','fontsize',10,'fontweight','b')
    zlabel('Pressure, P (atm)','color','k','fontsize',10,'fontweight','b')
    axis tight
    %
end % end of choose == 3
%
if choose == 4 % start of choose 4 : The contour plots
    %
    % Among nout points, take 5 time points spreading evenly between
    % [0,ta]
    %
    it = zeros(1,5); % Added this!
    for m = 1:5
        it(m) = 1+round((m-1)*(nout-1)/4); % Changed from it(m) = 1+(m-1)*(nout-1)/4; to it(m) = 1+ceil((m-1)*(nout-1)/4) or it(m) = 1+round((m-1)*(nout-1)/4);
    end
    % If nout = 101 ==> it = 1 26 51 76 101
    time = zeros(1,5); % Added this!
    for m = 1:5
        time(m) = t(it(m)); % Changed from time((m)) = t(it(m)); to time(m) = t(it(m));
    end
    %
    % Adsorption/Desorption profiles
    %
    % particle diameter in micrometer
    dp = dp * 10000;
    if adsorptionORpurge == 111
        % O2 initially present in the bed
        O2ingas = epsilon*A*L*(y10)*rhog0; % mmoles O2
        O2insolid = A*L*rhob*n10; % mmoles O2
        TotalO2inbed0_mmolesO2 = O2ingas + O2insolid; % mmoles O2
        % O2 remaining at every it, average over the whole bed i = 1:n
        ylacc = zeros(nout); % Added these initializations!
        rhogbaracc = zeros(nout);
        n1baracc = zeros(nout);
        y1ave = zeros(nout);
        rhogave = zeros(nout);
        n1ave = zeros(nout);
        O2remainingas = zeros(nout); % Added these six!
        O2remaininsolid = zeros(nout);
        TotalO2remaininbed_mmolesO2 = zeros(nout);
        frac_O2_desorbed = zeros(nout);
        y2acc = zeros(nout);
        n2baracc = zeros(nout);
        y2ave = zeros(nout);
        n2ave = zeros(nout);
        N2remainingas = zeros(nout);
        N2remaininsolid = zeros(nout);
        TotalN2remaininbed_mmolesN2 = zeros(nout);
        frac_N2_desorbed = zeros(nout);
        for it = 1:nout
            ylacc(it) = 0;
            rhogbaracc(it) = 0;
            n1baracc(it) = 0;
            for i = 1:n
                ylacc(it) = y1bar(it,i) + ylacc(it);
                rhogbaracc(it) = rhogbar(it,i) + rhogbaracc(it);
                n1baracc(it) = n1bar(it,i)+n1baracc(it);
            end
            y1ave(it) = ylacc(it)/n;
            rhogave(it) = (Qa*ta/L)*rhogbaracc(it)/n;
            n1ave(it) = (Qa/(k1*rhob*L))*n1baracc(it)/n;
        end
        for it = 1:nout
            O2remainingas(it) = epsilon*A*L*y1ave(it)*rhogave(it); % mmoles O2
            O2remaininsolid(it) = A*L*rhob*n1ave(it); % mmoles O2
            TotalO2remaininbed_mmolesO2(it) = O2remainingas(it)+O2remaininsolid(it); % mmoles O2
            frac_O2_desorbed(it) = (TotalO2inbed0_mmolesO2-TotalO2remaininbed_mmolesO2(it))/TotalO2inbed0_mmolesO2;
        end
        frac_O2_desorbedv = frac_O2_desorbed';
    end
    if adsorptionORpurge == 222
        % N2 initially present in the bed
        N2ingas = epsilon*A*L*(1-y10)*rhog0; % mmoles N2
        N2insolid = A*L*rhob*n20; % mmoles N2
        TotalN2inbed0_mmolesN2 = N2ingas + N2insolid; % mmoles N2
        % N2 remaining at every it, average over the whole bed i = 1:n
        for it = 1:nout
            y2acc(it) = 0;
            rhogbaracc(it) = 0;
            n2baracc(it) = 0;
            for i = 1:n
                y2acc(it) = (1-y1bar(it,i)) + y2acc(it);
                rhogbaracc(it) = rhogbar(it,i) + rhogbaracc(it);
                n2baracc(it) = n2bar(it,i)+n2baracc(it);
            end
            y2ave(it) = y2acc(it)/n;
            rhogave(it) = (Qa*ta/L)*rhogbaracc(it)/n;
            n2ave(it) = (Qa/(k2*rhob*L))*n2baracc(it)/n;
        end
        for it = 1:nout
            N2remainingas(it) = epsilon*A*L*y2ave(it)*rhogave(it); % mmoles N2
            N2remaininsolid(it) = A*L*rhob*n2ave(it); % mmoles N2
            TotalN2remaininbed_mmolesN2(it) = N2remainingas(it)+N2remaininsolid(it);
            % mmoles N2
            frac_N2_desorbed(it) = (TotalN2inbed0_mmolesN2-TotalN2remaininbed_mmolesN2(it))/TotalN2inbed0_mmolesN2;
        end
        frac_N2_desorbedv = frac_N2_desorbed';
    end
    %
    figure(1)
    plot(zbar,(y1bar(1,:)),zbar,(y1bar((nout-1)*0.2+1,:)),zbar,(y1bar((nout-1)*0.4+1,:)),zbar,(y1bar((nout-1)*0.6+1,:)),zbar,(y1bar((nout-1)*0.8+1,:)),zbar,(y1bar(nout,:)),'LineWidth',2,'MarkerSize',3); %axis tight
    % legend('tbar = 0.0','tbar = 0.2','tbar = 0.4','tbar = 0.6','tbar = 0.8','tbar = 1.0',1);
    legend('tbar = 0.0','tbar = 0.2','tbar = 0.4','tbar = 0.6','tbar = 0.8','tbar = 1.0','Location','S');
    text(.02,1.05,['\fontsize{11}d_p = ',num2str(dp),'\mum; \fontsize{11}L = ',num2str(L),' cm; \fontsize{11}D = ', num2str(D),'cm'])
    if adsorptionORpurge == 111
        title('Contour plot of IgG mole fraction','color','k','fontsize',12,'fontweight','b');
    end
    if adsorptionORpurge == 222
        title('Contour plot of IgG mole fraction','color','k','fontsize',12,'fontweight','b');
    end
    xlabel('Bed axial position,zbar','color','k','fontsize',12,'fontweight','b')
    ylabel('y_{IgG}','color','k','fontsize',12,'fontweight','b')
    axis([0 1 0 1.1])
    figure(2)
    plot(zbar,(1-y1bar(1,:)),zbar,(1-y1bar((nout-1)*0.2+1,:)),zbar,(1-y1bar((nout-1)*0.4+1,:)),zbar,(1-y1bar((nout-1)*0.6+1,:)),zbar,(1-y1bar((nout-1)*0.8+1,:)),zbar,(1-y1bar(nout,:)),'LineWidth',2,'MarkerSize',3); % axis tight
    % legend('tbar = 0.0','tbar = 0.2','tbar = 0.4','tbar = 0.6','tbar = 0.8','tbar = 1.0',1);
    % Replaced with below!
    legend('tbar = 0.0','tbar = 0.2','tbar = 0.4','tbar = 0.6','tbar = 0.8','tbar = 1.0','Location', 'N');
    text(.02,1.05,['\fontsize{11}d_p = ',num2str(dp),'\mum; \fontsize{11}L = ',num2str(L),' cm; \fontsize{11}D = ', num2str(D),'cm'])
    if adsorptionORpurge == 111
        title('Contour plot of Nitrogen molefraction','color','k','fontsize',12,'fontweight','b');
    end
    if adsorptionORpurge == 222
        title('Contour plot of Nitrogen molefraction','color','k','fontsize',12,'fontweight','b');
    end
    xlabel('Bed axial position,zbar','color','k','fontsize',12,'fontweight','b')
    ylabel('y_{Non-IgG}','color','k','fontsize',12,'fontweight','b')
    axis([0 1 0 1.1])
    figure(3)
    plot(zbar,((Qa/(k1*rhob*L)).*n1bar(1,:)),zbar,((Qa/(k1*rhob*L)).*n1bar((nout-1)*0.2+1,:)),zbar,((Qa/(k1*rhob*L)).*n1bar((nout-1)*0.4+1,:)),zbar,((Qa/(k1*rhob*L)).*n1bar((nout-1)*0.6+1,:)),zbar,((Qa/(k1*rhob*L)).*n1bar((nout-1)*0.8+1,:)),zbar,((Qa/(k1*rhob*L)).*n1bar(nout,:)),'LineWidth',2,'MarkerSize',3); %axis tight
    % legend('tbar = 0.0','tbar = 0.2','tbar = 0.4','tbar = 0.6','tbar = 0.8','tbar = 1.0',1);
    % Replaced with below!
    legend('tbar = 0.0','tbar = 0.2','tbar = 0.4','tbar = 0.6','tbar = 0.8','tbar = 1.0','Location','E');
    if adsorptionORpurge == 111
        title('Contour plot of adsorbed IgG','color','k','fontsize',12,'fontweight','b');
    end
    if adsorptionORpurge == 222
        title('Contour plot of adsorbed IgG','color','k','fontsize',12,'fontweight','b');
    end
    xlabel('Bed axial position,zbar','color','k','fontsize',12,'fontweight','b')
    ylabel('n_{O2} (mmol IgG/g adsorbent)','color','k','fontsize',12,'fontweight','b')
end