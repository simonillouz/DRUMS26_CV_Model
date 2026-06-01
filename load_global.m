% This function initializes the parameters for the model and sets initial
% values for the variables.
function [x0, parsG, times, Init, low, hi] = load_global(PID, data)

Qtot = data.Qtot;

% Flows (related to subject)
Qmval = Qtot;
Qaval = Qtot;
Qsys = Qtot;
Qtval = Qtot;
Qpval = Qtot;
Qpul = Qtot;

PsaS = data.PsaS;
PsaD = data.PsaD;
PsaM = PsaS/3 + 2*PsaD/3;

PpaS = data.PpaS;
PpaD = data.PpaD;
PpaM = PpaS/3 + 2*PpaD/3;

Psv = data.Psv;
Ppv = data.Ppv;

%Resistances (Ohm's law)
Rpul  = (PpaD - Ppv)/Qpul;
Rmval = 0.001;   % Mitral valve resistance
Raval = 0.001;   % Arterial valve resistance

Rsys   = (PsaD - Psv)/Qsys;
Rtval  = 0.001;
Rpval  = 0.001;

% Total Volumes (Beneken and deWit)

VtotA= 0.20*data.TotalVol; % Arterial volume (20%)
VtotV= 0.80*data.TotalVol; % Venous volume (80%)

% Peskin & Hoppensteadt pg 8 (21 on pdf)
VpvT = VtotV*0.80;   % Upper body arterial volume 
VsaT = VtotA*0.22;   % Lower body arterial volume 
VsvT = VtotV*0.78;   % Upper body venous volume 
VpaT = VtotA*0.20;   % Lower body venous volume  
                        
% Stressed Volumes - used Beneken and deWit for percentages
VpvS = VpvT*.075;     
VsaS = VsaT*.27;    
VsvS = VsvT*.105;
VpaS = VpaT*.58;

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
Vlvm = data.Vlvm;      % Minimum volume left ventricle (from literature)

VrvM = data.VrvM;     % Max volume right ventricle (from literature)
Vrvm = data.Vrvm;      % Minimum volume right ventricle (from literature)

Elvm  = Psv/VlvM;   % Minimum cardiac elastance - left heart
ElvM  = PsaS/Vlvm;  % Maximum cardiac elastance - left heart

Ervm = Psv/VrvM; % minimum cardiac elastance - right heart
ErvM = PsaS/Vrvm; % max cardiac elastance - right heart

TC  = data.TC;    % Time for maximum contractility relative to the lenght of the heart beat
TR  = data.TR;    % Time for relaxation relative to the lenght of the heart beat

% times
times.dt = 0.01;
times.T = round(data.HR/60/times.dt)*times.dt;
times.NC = 3;
times.td = 0:times.dt:times.T*times.NC;

Init = [VpvT VlvM VsaT VsvT VrvM VpaT]; % Steady state initial conditions

% Parameter vector
x0 = [Rpul Rmval Raval Rsys Rtval Rpval ...  %1-6
      Cpv Csa Csv Cpa ...          %7-10
      ElvM Elvm ErvM Ervm TC TR]';  %11-16

x0 = log(x0);

% MA432 Required for optimization
low = x0 - log(4);
hi  = x0 + log(4);
