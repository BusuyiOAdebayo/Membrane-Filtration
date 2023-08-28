% SPTFF Model Optimizer
% 3 Biochemical Components/species
close all; clc

tic

De_A = 1e-4;
De_B = 1e-6;
%De_C = 1e-6;

theta0 = [De_A, De_B];%, De_C]; % X0

% theta0 = [1e-4;1e-6;1e-6]; % X0

% For data in number matrix format
% t = SPTFFData(1:end,1); % XDATA
% c = SPTFFData(1:end,2:end); 

% For data in Table format
t = SPTFFData.Time; % XDATA
c = SPTFFData.cTFFA; 

% Below didn't work: Unrecognized table variable name 'ForFitting'.
% t = SPTFFData.ForFitting.Time; % XDATA
% c = SPTFFData.ForFitting.cTFFA; 

LB = (rand(3,1))';

for i = 1:3 
    LB(i) = 1e-8; %[];
    %UB(i) = 50; %[];
end
UB = []; % No upper bound!

% Solver-Based Solution Approach

% For lsqnonlin
% options = optimoptions('lsqcurvefit', 'OptimalityTolerance', 1e-10, 'FunctionTolerance', 1e-10);
%[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@SPTFFModelSimulator,theta0,t,c);
[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@SPTFFModelSimulator,theta0,t,c,LB,UB);
%[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(@SPTFFModelSimulator,theta0,t,c,LB,UB,options);

% For lsqnonlin
%[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqnonlin(@SPTFFModelSimulator,theta0);
%[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqnonlin(@SPTFFModelSimulator,theta0,LB,UB);
%[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqnonlin(@SPTFFModelSimulator,theta0,LB,UB,options);

fprintf(1,'\tRate Constants:\n')
for k1 = 1:length(theta)
    fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta(k1))
end

tfit = (linspace(min(t), max(t)))';
Cfit = SPTFFModelSimulator(theta, tfit);