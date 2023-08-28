% Solve for y1, n1, n2, P, rhog, Q, us
%
function yt = PSA_adsorption(~,Y) % yt = PSA_adsorption(t,Y)
%
% parameters shared with the ODE routine
%
global Mg1 Mg2 De1 De2 taup dYbardZbarDifferentialMethod rhog0 T0 ncall n dzbar dtbar epsilon epsilonbar k1 k2 rhob dp L tc tt y1a Qin Qa Pa Pmax Pmin Tmax Tmin Tref MassAxial dtbara dtbarp dtbardep dtbard miu ta tp td tdep
%
% preallocating
%
Y1bar = zeros(1,n);
n1bar = zeros(1,n);
n2bar = zeros(1,n);
Pbar = zeros(1,n); % Pressure
y1bar = zeros(1,n);
Mg = zeros(1,n);
thetabar = zeros(1,n);
rhogbar = zeros(1,n);
n1eqmbar = zeros(1,n);
n2eqmbar = zeros(1,n);
Qbar = zeros(1,n);
usbar = zeros(1,n);
convecy1bar = zeros(1,n);
DL = zeros(1,n);
dY1bardtbar = zeros(1,n);
dn1bardtbar = zeros(1,n);
dn2bardtbar = zeros(1,n);
dPbardtbar = zeros(1,n); % I put this!
Yt = zeros(1,4*n);
Qin = Qa;
%
% One vector to four vectors
%
for i = 1:n
    Y1bar(i) = Y(0*n+i);
    n1bar(i) = Y(1*n+i);
    n2bar(i) = Y(2*n+i);
    Pbar(i) = Y(3*n+i);
end
for i = 1:n
    y1bar(i) = Y1bar(i)/(rhog0*L/(Qa*tt));
    Mg(i) = Mg1*y1bar(i)+Mg2*(1-y1bar(i));
end
%
% Boundary conditions after initial conditions
%
if (ncall ~= 0) % not initial condition
    i = 1; % Only use boundary conditions at the adsorption inlet
    y1bar(i) = y1a;
    Y1bar(i) = y1a*(rhog0*L/(Qa*tt));
end
% Formulate pre-parameters for ODEs
%
% Isobaric, Isothermal
%
for i = 1:n
    Pbar(i) = Pa/(Pmax-Pmin);
    thetabar(i) = (T0-Tref)/(Tmax-Tmin);
    rhogbar(i) = rhog0*L/(Qa*tt);
end
%
% Formulate n1eqmbar(i), n2eqmbar(i), Qbar(i), usbar(i) with initial and boundary conditions
%
for i = 1:n
    py = (y1bar(i))*((Pmax-Pmin)*Pbar(i)/1013250);
    if py < 0.00013
        n1eqmbar(i) = 0;
    else
        n1eqmbar(i) = -(k1*rhob*L/Qa)*(3*0.069337*((y1bar(i))*((Pmax-Pmin)*Pbar(i)/1013250))^(-0.10825)/(1+0.069337*((y1bar(i))*((Pmax-Pmin)*Pbar(i)/1013250))^(-0.10825))-0.4633)/2/193.842*1000;
    end
    n2eqmbar(i) = 0;
end
if (ncall == 0)
    for i = 1:n
        Qbar(i) = Qin/Qa;
        usbar(i) = Qbar(i)/rhogbar(i);
    end
end
if (ncall ~= 0)
    for i = 1
        Qbar(i) = Qin/Qa;
        usbar(i) = Qbar(i)/rhogbar(i);
    end
    for i = 2:n
        Qbar(i) = -((n1eqmbar(i)-n1bar(i))+(n2eqmbar(i)-n2bar(i)))*dzbar+Qbar(i-1);
        usbar(i) = Qbar(i)/rhogbar(i);
    end
end
%
% Formulate gas phase axial dispersion coefficient in mass balance, DL, cm2/sec
%
if MassAxial == 881
    for i = 1:n
        DL(i) = 0;
    end
end
if MassAxial == 882
    DM = (epsilon/taup)*(1/(1/De1+1/De2));
    for i = 1:n
        DL(i) = 0.7*DM +0.5*dp*(abs(Qbar(i))*Qa)/(rhogbar(i)*(Qa*tt/L)*epsilonbar);
    end
end
%
%
% Start of 1st order (2 points) upwind difference
if dYbardZbarDifferentialMethod == 333
    %
    % Formulate convective terms for i = 3:n-2
    %
    for i = 3:n-2
        % Le = left edge of cell
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        % Re = right edge of cell
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = 0;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = 0;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    % Formulate convective terms for i = 1, 2, n-1, n
    %
    for i = 1
        Qbar_eps_Le = (Qbar(i)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = 1;
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = 1;
        end
        superbeey1barLe = 0;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = 0;
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = 0;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = 2
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = 0;
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = 0;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = 0;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = n-1
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = 0;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = 0;
        end
        superbeey1barRe = 0;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = n
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = 0;
        end
        superbeey1barLe = 0;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = 1;
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = 1;
        end
        superbeey1barRe = 0;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
end % end of dYbardZbarDifferentialMethod == 333
%
%
% Start of 2nd order (3 points) upwind difference
if dYbardZbarDifferentialMethod == 444
    %
    % Formulate convective terms for i = 3:n-2
    %
    for i = 3:n-2
        % Le = left edge of cell
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        % Re = right edge of cell
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = 1;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = 1;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    % Formulate convective terms for i = 1, 2, n-1, n
    %
    for i = 1
        Qbar_eps_Le = (Qbar(i)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = 1;
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = 1;
        end
        superbeey1barLe = 1;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = 0;
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = 1;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = 2
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = 0;
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = 1;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = 1;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = n-1
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = 1;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = 0;
        end
        superbeey1barRe = 1;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = n
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = 0;
        end
        superbeey1barLe = 1;
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = 1;
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = 1;
        end
        superbeey1barRe = 1;
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
end % end of dYbardZbarDifferentialMethod == 444
%
% Start of Superbee Flux Limiter
if dYbardZbarDifferentialMethod == 777
    %
    % Formulate convective terms for i = 3:n-2
    %
    for i = 3:n-2
        % Le = left edge of cell
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        % Re = right edge of cell
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = max(0,max(min(1,2*ry1barLe),min(2,ry1barLe)));
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = max(0,max(min(1,2*ry1barRe),min(2,ry1barRe)));
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    % Formulate convective terms for i = 1, 2, n-1, n
    %
    for i = 1
        Qbar_eps_Le = (Qbar(i)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = 1;
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = 1;
        end
        superbeey1barLe = max(0,max(min(1,2*ry1barLe),min(2,ry1barLe)));
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = 0;
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = max(0,max(min(1,2*ry1barRe),min(2,ry1barRe)));
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = 2
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = 0;
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = max(0,max(min(1,2*ry1barLe),min(2,ry1barLe)));
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = (y1bar(i+2)-y1bar(i+1))/(y1bar(i+1)-y1bar(i));
        end
        superbeey1barRe = max(0,max(min(1,2*ry1barRe),min(2,ry1barRe)));
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = n-1
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i+1))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = (y1bar(i+1)-y1bar(i))/(y1bar(i)-y1bar(i-1));
        end
        superbeey1barLe = max(0,max(min(1,2*ry1barLe),min(2,ry1barLe)));
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = (y1bar(i)-y1bar(i-1))/(y1bar(i+1)-y1bar(i));
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = 0;
        end
        superbeey1barRe = max(0,max(min(1,2*ry1barRe),min(2,ry1barRe)));
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i+1))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i+1)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
    for i = n
        Qbar_eps_Le = (Qbar(i-1)+Qbar(i))/(2*epsilon);
        Qbar_eps_Re = (Qbar(i)+Qbar(i))/(2*epsilon);
        % for left edge i-1/2
        if Qbar_eps_Le >= 0
            flowdirectionLe = 1;
            ry1barLe = (y1bar(i-1)-y1bar(i-2))/(y1bar(i)-y1bar(i-1));
        elseif Qbar_eps_Le <= 0
            flowdirectionLe = -1;
            ry1barLe = 0;
        end
        superbeey1barLe = max(0,max(min(1,2*ry1barLe),min(2,ry1barLe)));
        fluxy1barLe = 0.5*Qbar_eps_Le*((1+flowdirectionLe)*y1bar(i-1)+(1-flowdirectionLe)*y1bar(i))+0.5*abs(Qbar_eps_Le)*(1-abs(Qbar_eps_Le*dtbar/dzbar))*superbeey1barLe*(y1bar(i)-y1bar(i-1));
        % for right edge i+1/2
        if Qbar_eps_Re >= 0
            flowdirectionRe = 1;
            ry1barRe = 1;
        elseif Qbar_eps_Re <= 0
            flowdirectionRe = -1;
            ry1barRe = 1;
        end
        superbeey1barRe = max(0,max(min(1,2*ry1barRe),min(2,ry1barRe)));
        fluxy1barRe = 0.5*Qbar_eps_Re*((1+flowdirectionRe)*y1bar(i)+(1-flowdirectionRe)*y1bar(i))+0.5*abs(Qbar_eps_Re)*(1-abs(Qbar_eps_Re*dtbar/dzbar))*superbeey1barRe*(y1bar(i)-y1bar(i));
        convecy1bar(i) = (fluxy1barLe-fluxy1barRe)/dzbar;
    end
    %
end % end of dYbardZbarDifferentialMethod == 777
%
%
% Formulate ODEs - can be used for adsorption and purge - both directions
%
for i = 1
    dY1bardtbar(i) = convecy1bar(i)-(n1eqmbar(i)-n1bar(i))/epsilon+(epsilonbar*DL(i)*tt/(epsilon*L^2))*(Y1bar(i+2)-2*Y1bar(i+1)+Y1bar(i))/(dzbar)^2;
end
for i = 2:n-1
    dY1bardtbar(i) = convecy1bar(i)-(n1eqmbar(i)-n1bar(i))/epsilon+(epsilonbar*DL(i)*tt/(epsilon*L^2))*(Y1bar(i+1)-2*Y1bar(i)+Y1bar(i-1))/(dzbar)^2;
end
for i = n
    dY1bardtbar(i) = convecy1bar(i)-(n1eqmbar(i)-n1bar(i))/epsilon+(epsilonbar*DL(i)*tt/(epsilon*L^2))*(Y1bar(n)-2*Y1bar(n-1)+Y1bar(n-2))/(dzbar)^2;
end
%
if MassAxial == 882 % Turn on mass axial dispersion
    for i = n
        dY1bardtbar(i) = convecy1bar(i)-(n1eqmbar(i)-n1bar(i))/epsilon + 0;
    end
end
%
for i = 1:n
    dn1bardtbar(i) = k1*tt*(n1eqmbar(i)-n1bar(i));
    dn2bardtbar(i) = k2*tt*(n2eqmbar(i)-n2bar(i));
    dPbardtbar(i) = 0;
end
%
% Four vectors into one vector
%
for i = 1:n
    Yt(0*n+i) = dY1bardtbar(i);
    Yt(1*n+i) = dn1bardtbar(i);
    Yt(2*n+i) = dn2bardtbar(i);
    Yt(3*n+i) = dPbardtbar(i);
end
Yt(0*n+1) = 0; % because Y1bar(1) = y1a*rhog(1)
%
yt = Yt';
% Increment calls to psa_1
ncall = ncall+1;
%