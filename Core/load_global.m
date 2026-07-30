% This function initializes the parameters for the model and sets initial
% values for the variables.
function [pars, parsG, times, Init, low, hi] = load_global(PID,data)

Qtot = data.Qtot;

% Flows (related to subject)
Qmval = Qtot;
Qaval = Qtot;
Qsys  = Qtot;
Qtval = Qtot;
Qpval = Qtot;
Qpul  = Qtot;

PsaS = data.PsaS;
PsaD = data.PsaD;
PsaM = PsaS/3 + 2*PsaD/3;

PpaS = data.PpaS;
PpaD = data.PpaD;
PpaM = PpaS/3 + 2*PpaD/3;

Psv = data.Psv;
Ppv = data.Ppv;

%Resistances (Ohm's law)
Rpul   = (PpaM - Ppv)/Qpul;
if Rpul < 0
    PpaM = data.PrvS*0.98;
    Rpul = (PpaM - Ppv)/Qpul;
end

Rmval  = (Ppv - 0.99*Ppv)/Qmval;   %0.001;   % Mitral valve resistance
Raval  = (PsaS*1.01 - PsaS)/Qaval; %0.001;   % Arterial valve resistance
Rsys   = (PsaM - Psv)/Qsys;
Rtval  = (Psv - 0.99*Psv)/Qtval;   %0.001;
Rpval  = (PpaS*1.01 - PpaS)/Qpval; %0.001;

% Total Volumes (Beneken and deWit)
% Peskin & Hoppensteadt pg 8 (21 on pdf)
VpvT = data.TotalVol*0.125;    % Upper body arterial volume 
VsaT = data.TotalVol*0.13;     % Lower body arterial volume 

RV = data.VrvM/data.TotalVol;
LV = data.VlvM/data.TotalVol;
VsvT = (0.715 - RV - LV) * data.TotalVol;
VpaT = data.TotalVol*0.03;   % Lower body venous volume  

% disp([data.TotalVol VpvT+VsaT+VsvT+VpaT+data.VlvM+data.VrvM])
% disp([VpvT VsaT VsvT VpaT])
% pause;

% Stressed Volumes - used Beneken and deWit for percentages 
% --> stressed volume over total volume
VpvS = VpvT*.11;  % 54/514   = 0.1051 
VsaS = VsaT*.27;  % 160/585  = 0.2735  
VsvS = VsvT*.075; % 219/2916 = 0.0751 
VpaS = VpaT*.58;  % 69/119   = 0.5798

% Unstressed Volumes
parsG.VpvU = VpvT - VpvS;
parsG.VsaU = VsaT - VsaS;
parsG.VsvU = VsvT - VsvS;
parsG.VpaU = VpaT - VpaS;
parsG.VlvU = data.VlvU;
parsG.VrvU = data.VrvU;

% Compliances, stressed volume percentages from Beneken are weighted averages
Cpv = VpvS/Ppv;   % pulmonary veins compliance
Csa = VsaS/PsaS;  % systemic arteries compliance  
Csv = VsvS/Psv;   % Lower body arterial compliance
Cpa = VpaS/PpaS;  % Lower body venous compliance

VlvM = data.VlvM;     % Max volume left ventricle (from literature)
Vlvm = data.Vlvm;     % Minimum volume left ventricle (from literature)

VrvM = data.VrvM;     % Max volume right ventricle (from literature)
Vrvm = data.Vrvm;     % Minimum volume right ventricle (from literature)

Elvm  = Ppv/(VlvM -parsG.VlvU); % Minimum cardiac elastance - left heart
ElvM  = PsaM/(Vlvm-parsG.VlvU); % Maximum cardiac elastance - left heart

Ervm = Psv/(VrvM -parsG.VrvU);  % minimum cardiac elastance - right heart
ErvM = PpaM/(Vrvm-parsG.VrvU);  % max cardiac elastance - right heart

TC  = data.TC;    % Time for maximum contractility relative to the lenght of the heart beat
TR  = data.TR;    % Time for relaxation relative to the lenght of the heart beat

%disp('Check')
%disp([VlvM-Vlvm,VrvM-Vrvm,data.SV]);
%pause;

% times
times.dt = 0.01;
times.T = round(1/data.HR*60/times.dt)*times.dt;

% Parameter vector
pars = [Rpul Rmval Raval Rsys Rtval Rpval ...  %1-6
        Cpv Csa Csv Cpa ...                    %7-10
        ElvM Elvm ErvM Ervm TC TR]';           %11-16

%[Init] = initial_volumes(pars, parsG, data);
%disp([Init(1) VpvT]);
%disp([Init(2) VlvM]);
%disp([Init(3) VsaT]);
%disp([Init(4) VsaT]);

Init = [VpvT VlvM VsaT VsvT VrvM VpaT];
%disp([VlvM Vlvm VrvM Vrvm VpvT VsaT VsvT VpaT ])

pars = log(pars);
low = pars - log(4);
%low(10) = pars(10) - log(8); %pushes lower bound even lower
%low(i) = low(i) - log(?)
hi  = pars + log(4);
