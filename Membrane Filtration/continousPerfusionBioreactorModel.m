% Perfusion Bioreactor Model
function dcdt = continousPerfusionBioreactorModel(t,c)
c1target = 100; % Viable cells density, [E6 Cells/mL]
PR = 2;
CSPR = 0.05e-6;
%bleedRateConstant = 1e-2;
% Kp = -0.2;
% Ti = 0.5;
mugmax = 3*0.481360322;
mudmax = 0.05;
mudmin = 0.1*mudmax;
%mulmax = 0.20;
%mulmin = 0.20;
%Kl2 = 1.00;
Ks4 = 0.10;
Ks5 = 0.25*Ks4;
Ks6 = 0.25*Ks4;
Ks10 = 0.50*Ks4;
%Ks10app = Ks10+(Ks10/Ks4)*c(4);
Ki7 = 50;
Ki8 = 10;
Kd7 = 25;
Kd8 = 5;
Km4 = 0.005;
Km5 = 0.00025*Km4;
Km6 = 0.025*Km4;
Km10 = 0.5*Km4;
Kdi5 = 0.000025;
Kdi6 = 0.025;
Kdi10 = 0.025;
%Ku7 = 0.005;
% if c(4) < c4to7
%    Ku7 = 0.005;
% else
%    Ku7 = 0.00;
% end
if t <= 4
    Ku7 = 0;
else
    Ku7 = 0.0005;
end
Ku8 = 0.005;
% KlO2 = 1.75;
% KlN2 = 0.75;
% klCO2 = 1.75;
q4max = 0.25*0.068597286;
q4 = q4max*(c(4)/(Ks4+c(4)));
q5f = 0.011911535; % Glutamine formation
q5cmax = 0.05*q5f;
q5c = q5cmax*(c(5)/(Ks5+c(5))); % Glutamine consumption for cell growth
%q5 = 0.4*0.011911535*(Ks5/(Ks5+c(5))); % q5 = Y5/11*q11
q6max = 0.010453245;
q6 = q6max*(c(6)/(Ks6+c(6)));
q7 = 0.06821422;
%q7 = 0.06821422*(Ki7/(Ki7+c(7))); % q7 = Y7/4*q4
%q8 = 0.07413363;
%q8 = 0.07413363*(Ki8/(Ki8+c(8))); % q8 = Y8/5,6&11*q5,6&11 = Y8/5*q5+Y8/6*q6+Y8/11*q11
q9 = 0.00578808; % q9 = Y9/1*q1 = Y9/1*mug; q9 = q9max
q10max = 0.0125*0.005437439;
q10 = q10max*(c(10)/(Ks10+c(10)));
%q10 = q10max*(c(10)/(Ks10app+c(10)));
%q10 = q10max*(c(10)/(Ks10+c(10)))*((1-q4/q4max)^2);
%q10 = q10max*(c(10)/(Ks10app+c(10)))*((1-q4/q4max)^2);
q11 = 1.00;
c1feed = 0.00;
c2feed = 0.00;
c3feed = c1feed+c2feed;
%c4feed = 9.030312*1000/180.156;
c4feed = 9.030312;
c5feed = 0.00;
c6feed = 0.25*c4feed;
c7feed = 0.00;
c8feed = 0.00;
c9feed = 0.00;
c10feed = 2.5269;
% c11feed = 21; % N2
% c12feed = 78; % O2
% c13feed = 1; % CO2
% c11sat = 100;
% c12sat = 100;
% c13sat = 100;

% if t <= 6
%     Fharvest = CSPR*c(1)*c(11); % [mL/d];
% else
%     Fharvest = PR*c(11); % [mL/d];
% end
%
% Fmedia = Fharvest + Fbleed; % [mL/d]

%         if t <= 6
%             Fmedia = CSPR*c(1)*c(11); % [mL/d];
%         else
%             Fmedia = PR*c(11); % [mL/d];
%         end

% if t == 0
%     ebleed = 0;
% else
%     ebleed = c1target - c(1);
% end

% % t = tinit + n*deltat;
% ebleed = zeros(length(t),1);
% MVbleed = zeros(length(t),1);
% Fbleed = zeros(length(t),1);
% for i = 1:length(t)
%     ebleed(i) = c1target - c(i,1);
%     MVbleed(i) = Kp*(ebleed(i) + 1/Ti*(ebleed(i) - ebleed(i-1)));
%     Fbleed(i) = max(0, min(Fmedia(i), Fbleed(i-1) + MVbleed(i)));
% end

% ebleedti = c1target - c(1);
% MVbleedti = Kp*ebleedti + Kp/Ti*(ebleedti - ebleedti-1);
% Fbleed = max(0, min(Fmedia, Fbleedti-1 + MVbleedti));
% ebleedti = Xvtarget - Xvti;
% dbleedti = Kp*ebleedti + Kp/Ti*(ebleedti - ebleedti-1);
% Fbti = max(0, min (Ff, Fbti-1 + dbleedti));
% where ebleedti is the deviation at time i of the cell concentration (Xvti) from the target setpoint (Xvtarget),
% MVbleedti is the controller output (with Kp and Ti as proportional and integral terms) that will be used to adjust
% the bleeding rate imposed at time i-1 (Fbti-1) such as to maintain the target setpoint (Xvtarget). For this study,
% the PI control parameters have been hand tuned and set to Kp = −0.2 and Ti = 0.5.

% % time = [];
% c1Error = zeros(length(tSpan),1);
% controllerOutput = zeros(length(tSpan),1);
% Fbleed = zeros(length(tSpan),1);
% for j = 1:length(tSpan)
%     c1Error(j) = c1target - c(j,1); % c(1) error
%     controllerOutput(j) = Kp*(c1Error(j) + 1/Ti*(c1Error(j) - c1Error(j-1)));
%     Fbleed(j) = max(0, min(Fmedia(j), Fbleed(j-1) + controllerOutput(j)));
% end
% This is time event, try to avoid it!
% if t <= 10
%     Fbleed = 0;
% else
%     Fbleed = bleedRateConstant*c(1); %Fmedia/1000; % Need to still cork on this!
% end

%         % This is state event, use this!
%         if c(1) <= c1target
%             Fbleed = 0;
%         else
%             % solve((Fmedia/c(11))*c1feed-(Fbleed/c(11))*c(1)+r1-(c(1)/c(11))*dcdt(11) == 0, Fbleed);
%             Fbleed = (c(11)/c(1))*((Fmedia/c(11))*c1feed+r1-(c(1)/c(11))*dcdt(11));
%         end

%         Fharvest = Fmedia - Fbleed; % [mL/d]

% V = 2.5; %2.5/1000; % 2.5 L = 2.5/1000 m3
% alphaLac = 0.55;
% alphaAmm = 0.15;
% alphaCO2 = 1-(alphaLac+alphaAmm);
% betta = 1e-3;

mug = mugmax*(c(4)/(Ks4+c(4)))*(c(5)/(Ks5+c(5)))*(c(6)/(Ks6+c(6)))*(c(10)/(Ks10+c(10)))*(Ki7/(Ki7+c(7)))*(Ki8/(Ki8+c(8)));
%mug = mugmax*(c(4)/(Ks4+c(4)))*(c(10)/(Ks10+c(10)))*(Ki7/(Ki7+c(7)))*(Ki8/(Ki8+c(8)));
%mug = mugmax*((c(4)/(Ks4+c(4)))+(c(10)/(Ks10+c(10))))*(Ki7/(Ki7+c(7)))*(Ki8/(Ki8+c(8)));
mud = mudmin*(Kd7/(Kd7+c(7)))*(Kd8/(Kd8+c(8)))+mudmax*(c(7)/(Kd7+c(7)))*(c(8)/(Kd8+c(8)));
%mud = mudmin+mudmax*(Kd7/(Kd7+c(7)))*(Kd8/(Kd8+c(8)));
%mul = mulmin*(Kl2/(Kl2+c(2)))+mulmax*(c(2)/(Kl2+c(2)));
r1 = (mug-mud)*c(1); % mug*c(1)-mud*c(1)

% if t <= 6
%     Fmedia = CSPR*c(1)*c(11); % [mL/d];
% else
%     Fmedia = PR*c(11); % [mL/d];
% end
% % This is a state event, use this!
% if c(1) <= c1target
%     Fbleed = 0;
% else
%     % solve((Fmedia/c(11))*c1feed-(Fbleed/c(11))*c(1)+r1-(c(1)/c(11))*dcdt(11) == 0, Fbleed);
%     Fbleed = (c(10)/c(1))*((Fmedia/c(11))*c1feed+r1);%-(c(1)/c(11))*r11);
% end
% Fharvest = Fmedia - Fbleed; % [mL/d]

Fmedia = CSPR*c(1)*c(10); % [mL/d]
Fharvest = CSPR*c(1)*c(10); % [mL/d];
Fbleed = Fmedia-Fharvest;

r2 = mud*c(1);
r3 = r1+r2;
r4 = -q4*c(1)-Km4*c(1);
%r5 = q5*c(1)-Km5*c(1)-Kdi5*c(5);
r5 = q5f*c(8)*c(10)-q5c*c(1)-Km5*c(1)-Kdi5*c(5);
r6 = -q6*c(1)-Km6*c(1)-Kdi6*c(6);
%r7 = q7*c(1)-Ku7*c(1);
r7 = q7*c(4)-Ku7*c(1);
%r8 = q8*c(1)+Kdi5*c(5)+Kdi6*c(6)+Kdi10*c(10);
%r8 = q8*c(1);
r8 = Kdi5*c(5)+Kdi6*c(6)+Kdi10*c(10)-Ku8*c(8);
r9 = q9*c(1);
r10 = -q10*c(1)-Km10*c(1)-Kdi10*c(10);
r11 = q11*(Fmedia-Fharvest-Fbleed);
r12 = 0;
dcdt = zeros(12,1);
%         dcdt(1) = (Fmedia/c(11))*c1feed-(Fharvest/c(11))*c(1)*0-(Fbleed/c(11))*c(1)+r1-(c(1)/c(11))*dcdt(11);
%         dcdt(2) = (Fmedia/c(11))*c2feed-(Fharvest/c(11))*c(2)*0-(Fbleed/c(11))*c(2)+r2-(c(2)/c(11))*dcdt(11);
%         dcdt(3) = (Fmedia/c(11))*c3feed-(Fharvest/c(11))*c(3)*0-(Fbleed/c(11))*c(3)+r3-(c(3)/c(11))*dcdt(11);
dcdt(1) = (Fmedia/c(11))*c1feed-(Fbleed/c(11))*c(1)+r1-(c(1)/c(11))*dcdt(11);
dcdt(2) = (Fmedia/c(11))*c2feed-(Fbleed/c(11))*c(2)+r2-(c(2)/c(11))*dcdt(11);
dcdt(3) = (Fmedia/c(11))*c3feed-(Fbleed/c(11))*c(3)+r3-(c(3)/c(11))*dcdt(11);
dcdt(4) = (Fmedia/c(11))*c4feed-(Fharvest/c(11))*c(4)-(Fbleed/c(11))*c(4)+r4-(c(4)/c(11))*dcdt(11);
dcdt(5) = (Fmedia/c(11))*c5feed-(Fharvest/c(11))*c(5)-(Fbleed/c(11))*c(5)+r5-(c(5)/c(11))*dcdt(11);
dcdt(6) = (Fmedia/c(11))*c6feed-(Fharvest/c(11))*c(6)-(Fbleed/c(11))*c(6)+r6-(c(6)/c(11))*dcdt(11);
dcdt(7) = (Fmedia/c(11))*c7feed-(Fharvest/c(11))*c(7)-(Fbleed/c(11))*c(7)+r7-(c(7)/c(11))*dcdt(11);
dcdt(8) = (Fmedia/c(11))*c8feed-(Fharvest/c(11))*c(8)-(Fbleed/c(11))*c(8)+r8-(c(8)/c(11))*dcdt(11);
dcdt(9) = (Fmedia/c(11))*c9feed-(Fharvest/c(11))*c(9)-(Fbleed/c(11))*c(9)+r9-(c(9)/c(11))*dcdt(11);
dcdt(10) = (Fmedia/c(11))*c10feed-(Fharvest/c(11))*c(10)-(Fbleed/c(11))*c(10)+r10-(c(10)/c(11))*dcdt(11);
dcdt(11) = r11;
dcdt(12) = ((1.0/c(3))*dcdt(1)-(c(1)/c(3)^2)*dcdt(3))*100+r12;
% dcdt(13) = (Fperfusion/c(11))*c13feed-(Foutlet/c(11))*c(13)-(Fbleed/c(11))*c(13)+r13-(c(13)/c(11))*dcdt(11); % N2
% dcdt(14) = (Fperfusion/c(11))*c14feed-(Foutlet/c(11))*c(14)-(Fbleed/c(11))*c(14)+r14-(c(14)/c(11))*dcdt(11); % O2
% dcdt(15) = (Fperfusion/c(11))*c15feed-(Foutlet/c(11))*c(15)-(Fbleed/c(11))*c(15)+r15-(c(15)/c(11))*dcdt(11); % CO2
% dcdt(16) = alphaLac*c(7)+alphaAmm*c(8)+alphaCO2*c(15)+r16;
% dcdt(17) = betta*c(17)+r17;
end

% fileName = 'output.xlsx';
% sheetName 'sheet1';
% xlswrite(fileName,C,sheetName,'c(1)');