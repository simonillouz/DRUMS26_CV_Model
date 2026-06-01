% Copyright (C) 2019 Mette S. Olufsen.
% All Rights Reserved.
% 
% Note: This a trial version built to develop and test a research version.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software") for research and
% development purposes subject to the following conditions:
% 
% The above copyright notice and the README.txt file shall be included in all copies of 
% any portion of this data.
% 
% Whenever reasonable and possible in publications and presentations when this model or 
% data is used in whole or part, please include an acknowledgement to the
% manuscript listed below and cite the official website for this code: https://wp.math.ncsu.edu/cdg/
% The software can befound under supplemental material  
%
% Explicit written permission should be obtained from M.S. Olufsen prior to posting or 
% distributing any portion of this software. Please email msolufse@ncsu.edu to obtain
% written permission.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
% IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
function  DriverBasic(PID,opt)
% Solves the five compartment model described in detail in the manuscript
% by Williams, Brady, Gilmore, Gremaud, Tran, Ottesen, Mehlsen and Olufsen,
% in J Math Biol, Cardiovascular dynamics during head-up tilt assessed via
%pulsatile and non-pulsatilemodels, https://doi.org/10.1007/s00285-019-01386-9
close all;
global ODE_TOL % absolute and relative tolerance for ODE solver


% [mdata,data] = LoadData(datn);

data = Patient(PID);

%Load nominal parameter values (pars) and initial conditions for the ODEs
[pars parsG times Init] = load_global(PID, data);

pars = exp(pars); 

% Load optimized parameters if opt = 1                
% if opt == 1
%   s = strcat('../Optimization/Opt',num2str(datn),'.mat');
%   load(s);
%   pars = exp(parsLM);
% end

% Display parameter values
% if opt == 1
%     disp('Estimated parameters')
% else
%     disp('Nominal parameters')
% end

% pars;

% Define resistances
Rpul  = pars(1);
Rmval = pars(2);
Raval = pars(3);
Rsys   = pars(4);
Rtval  = pars(5);
Rpval  = pars(6);

% Define compliances
Cpv = pars(7);
Csa	= pars(8);
Csv	= pars(9);
Cpa = pars(10);

% Heart parameters
ElvM = pars(11);
Elvm = pars(12);
ErvM = pars(13);
Ervm = pars(14);

TC = pars(15); 
TR = pars(16);

%Initialize empty vectors for solutions for pressure
PpvSa  = []; 
PlvSa  = []; 
PsaSa  = []; 
PsvSa  = [];
PrvSa  = [];
PpaSa  = [];

%Initialize empty vectors for solutions for volume
VpvSa  = []; 
VlvSa  = []; 
VsaSa  = []; 
VsvSa  = [];
VrvSa  = [];
VpaSa  = [];

%Initialize empty vectors for solutions for flow
QpulSa  = [];
QmvalSa = [];
QavalSa = []; 
QtvalSa = []; 
QpvalSa = []; 
QsysSa  = []; 

options=odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);  %sets how accurate the ODE solver is

T  = times.T;
k1 = 1; % index of first time step in first period
k2 = round(T/times.dt)+k1; %index of last time step in first period
% Ts = length(times.td(k1:k2));
for i = 1:times.NC % go through loop NC times
    clear Vpv Vlv Vsa Vsv Vrv Vpa 
    clear Ppv Plv Psa Psv Prv Ppa Elv Erv
    clear Qpul Qmval Qaval Qtval Qpval Qsys
    
    tdc = times.td(k1:k2); %current period

    sol = ode15s(@modelBasic,tdc,Init,options,pars, parsG, tdc(1)); %the ODE solver calls modelBasic and enters in all the following values
    sols= deval(sol,tdc); %takes the solutions at each time step in the period
 
    
    %assigns each row as a temporary vector to store the solutions
    %for the current time period
    Vpv  = sols(1,:)'; 
    Vlv  = sols(2,:)'; 
    Vsa  = sols(3,:)'; 
    Vsv  = sols(4,:)'; 
    Vrv  = sols(5,:)';
    Vpa  = sols(6,:)';
   
    Ppv = (Vpv-parsG.VpvU)/Cpv;
    Psa = (Vsa-parsG.VsaU)/Csa;
    Psv = ((Vsv-parsG.VsvU)/Csv)';
    Ppa = (Vpa-parsG.VpaU)/Cpa;
    
    %Determines the elasticity of the left ventricle at each period timestep
    for j = 1:length(tdc)
        Elv(j) = Elastance(tdc(j)-tdc(1),ElvM, Elvm, TC, TR);
    
        %Pressure of the left ventricle
        Plv(j) = Elv(j).*(Vlv(j) - parsG.VlvU); 

        % flow into LH (mitral valve)
        if Ppv(j) > Plv(j)
            Qmval(j) = (Ppv(j)-Plv(j))/Rmval;
        else
            Qmval(j) = 0;
            
        end
        % flow out of LH (aortic valve)
        if Plv(j) > Psa(j)
            Qaval(j) = (Plv(j)-Psa(j))/Raval;
        else
            Qaval(j) = 0;
        end 
    end
    Psa = Psa';

    %Determines the elasticity of the right ventricle at each period timestep
    for j = 1:length(tdc)
        Erv(j) = Elastance(tdc(j)-tdc(1),ErvM, Ervm, TC, TR);
    
        %Pressure of the right ventricle
        Prv(j) = Erv(j).*(Vrv(j) - parsG.VrvU); 

       % flow into RH (tricuspid valve)
        if Psv(j) > Prv(j)
            Qtval(j) = (Psv(j)-Prv(j))/Rtval;
        else
            Qtval(j) = 0;
        end
        % flow out of RH (pulmonary valve)
        if Prv(j) > Ppa(j)
            Qpval(j) = (Prv(j)-Ppa(j))/Rpval;
        else
            Qpval(j) = 0;
        end
    end
    
    %Flows defined by Ohm's Law
    Qpul = ((Ppa - Ppv)/Rpul)';
    Qsys = (Psa - Psv)/Rsys;

    % Adds to pressure solution vectors every iteration of the loop
    PpvSa  = [PpvSa  Ppv(1:end-1)'];
    PlvSa  = [PlvSa  Plv(1:end-1)];
    PsaSa  = [PsaSa  Psa(1:end-1)];
    PsvSa  = [PsvSa  Psv(1:end-1)];
    PrvSa  = [PrvSa  Prv(1:end-1)]; 
    PpaSa =  [PpaSa  Ppa(1:end-1)'];

    
    % Adds to volume solution vectors every iteration of the loop
    VpvSa  = [VpvSa  Vpv(1:end-1)'];
    VlvSa  = [VlvSa  Vlv(1:end-1)'];
    VsaSa  = [VsaSa  Vsa(1:end-1)'];
    VsvSa  = [VsvSa  Vsv(1:end-1)'];
    VrvSa  = [VrvSa  Vrv(1:end-1)']; 
    VpaSa =  [VpaSa  Vpa(1:end-1)'];
    
    % Adds to flow solution vectors every iteration of the loop
    QpulSa  = [QpulSa   Qpul(1:end-1)];
    QmvalSa = [QmvalSa  Qmval(1:end-1)];
    QavalSa = [QavalSa  Qaval(1:end-1)];
    QtvalSa = [QtvalSa  Qtval(1:end-1)];
    QpvalSa = [QpvalSa  Qpval(1:end-1)];
    QsysSa  = [QsysSa   Qsys(1:end-1)];
    
    Init = [Vpv(end) Vlv(end) Vsa(end) Vsv(end) Vrv(end) Vpa(end)];
    
    % if i < 20
    %   T  = (data.per(i+1)-data.per(i));
    % 
    %   k1 = k2; %sets last index of this loop as first index for next loop
    %   k2 = round(T/times.dt)+k1;
    %   if k2 > length(data.td)
    %       k2 = length(data.td);
    %   end
    % end 
    k1 = k2;
    k2 = round(times.T/times.dt)+k1;
end

%After for loop is done, saves last value for each of these vectors
PpvSa  = [PpvSa  Ppv(end)];
PlvSa  = [PlvSa  Plv(end)];
PsaSa  = [PsaSa  Psa(end)];
PsvSa  = [PsvSa  Psv(end)];
PrvSa  = [PrvSa  Prv(end)]; 
PpaSa =  [PpaSa  Ppa(end)];

VpvSa  = [VpvSa  Vpv(end)];
VlvSa  = [VlvSa  Vlv(end)];
VsaSa  = [VsaSa  Vsa(end)];
VsvSa  = [VsvSa  Vsv(end)];
VrvSa  = [VrvSa  Vrv(end)]; 
VpaSa =  [VpaSa  Vpa(end)];

QpulSa  = [QpulSa   Qpul(end)];
QmvalSa = [QmvalSa  Qmval(end)];
QavalSa = [QavalSa  Qaval(end)];
QtvalSa = [QtvalSa  Qtval(end)];
QpvalSa = [QpvalSa  Qpval(end)];
QsysSa  = [QsysSa   Qsys(end)];

figure(20);clf;
h=plot(times.td,PpaSa,'b');
set(gca,'fontsize',24);
set(h,'linewidth',2);
xlabel('time (s)');
ylabel('Ppa (mmHg)');
legend('data','model');
grid on;
% if opt == 1
%   spf2 = strcat('modelPau_opt',num2str(datn),'.png');
% else
%   spf2 = strcat('modelPau_nom',num2str(datn),'.png');
% end
% print(spf2,'-dpng');

N = length(times.td);
N = N-min(6*length(tdc),length(times.td))+1;
% figure(30);clf;
% h=plot(data.td(N:end),PpaSa(N:end),data.td(N:end),PsaSa(N:end),...
%        data.td(N:end),PsvSa(N:end),data.td(N:end),PpvSa(N:end),...
%        data.td(N:end),PlvSa(N:end), data.td(N:end), PrvSa(N:end));
% set(gca,'fontsize',24);
% set(h,'linewidth',2);
% xlabel('time (s)');
% ylabel('pressure (mmHg)');
% xlim([data.td(N) data.td(end)]);
% legend('Ppa','Psa','Psv','Ppv','Plv', 'Prv');
% grid on;
% if opt == 1
%   spf3 = strcat('modelAllPres_opt',num2str(datn),'.png');
% else
%   spf3 = strcat('modelAllPres_nom',num2str(datn),'.png');
% end
% print(spf3,'-dpng');
% 
% ss = strcat('Res',num2str(datn),'.mat');
% save(ss,'data','PpaSa');
