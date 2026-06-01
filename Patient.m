function data = Patient(PID)

ODE_TOL  = 1e-8;
DIFF_INC = 1e-4;

H      = 165; % Subject height cm
W      = 65; % Subject weight kg
Gender = 1; % Male 1, female 2

BSA = sqrt(H*W/3600);  % Body surface area

if Gender == 2
  data.TotalVol = (3.47*BSA - 1.954)*1000; % Female
elseif Gender == 1
  data.TotalVol = (3.29*BSA - 1.229)*1000; % Male
end
data.Qtot = data.TotalVol/60; % Assume blood circulates in one min

data.PsaS = 120;
data.PsaD = 80;

data.PpaS = 25;
data.PpaD = 10;

data.Psv = 2;
data.Ppv = 5;

data.VlvM = 120; % left heart max volume
data.Vlvm = 60;  

data.VrvM = 120;
data.Vrvm = 60;

% unstressed volumes (healthy patients)
data.VlvU = 10;
data.VrvU = 10;

data.TC  = 0.12;    % Time for maximum contractility relative to the lenght of the heart beat
data.TR  = 0.18;    % Time for relaxation relative to the lenght of the heart beat

data.HR = 60; % beats per min