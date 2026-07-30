function data = Patient(PID, teich)
% PID: instert patient id
% teich: 1 for on, 0 for off

data.ODE_TOL  = 1e-5;
data.DIFF_INC = sqrt(data.ODE_TOL);

% Loading Patient Data
patients = readtable("AllPatsMRI_2_REU.xlsx", 'VariableNamesRange', 'A1:D1', 'DataRange', 'A3:D345'); % loading patient ID, height, weight, and sex
rhcdata  = readtable("AllPatsMRI_2_REU.xlsx", 'VariableNamesRange', 'E1:N1', 'DataRange', 'E3:N345'); % loading patient data from RHC
ttedata  = readtable("AllPatsMRI_2_REU.xlsx", 'VariableNamesRange', 'Q1:R1', 'DataRange', 'Q3:R345'); % loading patient data from TTE
mridata  = readtable("AllPatsMRI_2_REU.xlsx", 'VariableNamesRange', 'AI1:AM1', "DataRange", 'AI3:AM345'); % loading patient data from MRI

data.H   = patients{PID, 'H'}; % Subject height cm
data.W   = patients{PID, 'W'};  % Subject weight kg
data.Sex = patients{PID, "Sex"};   % Male 1, female 2

Hgt = data.H/100; % Calculating height in m for below calculations

% Calculations for total volume - (Nader, Hidalgo & Bloch, 1961)
if data.Sex == 2
  TotalVol = 0.414 * Hgt^3 + 0.0328 * data.W - 0.030; % Female
elseif data.Sex == 1
  TotalVol = 0.417 * Hgt^3 + 0.0450 * data.W - 0.030; % Male
end

data.TotalVol = TotalVol*1000;

% Heart rate
data.HR = rhcdata{PID, "HR"};       % beats per min

% Cardiac output
data.CO   = rhcdata{PID, "CO"}; % patient cardiac output, L/min
data.Qtot = data.CO*1000/60; % Assume blood circulates in one min

% Stroke volume
data.SV = data.CO/data.HR*1000;

% Pressures
data.PsaS = rhcdata{PID, "PsaS"};
data.PsaD = rhcdata{PID, "PsaD"};

data.PpaS = rhcdata{PID, "PpaS"};
data.PpaD = rhcdata{PID, "PpaD"};

data.PrvS = rhcdata{PID, 'PrvS'};
data.PrvD = rhcdata{PID, 'PrvD'};

data.PpcW = rhcdata{PID, 'PpcW'}; 

% TTE data
data.IDlvS = ttedata{PID, 'IDLVS'}; % measured in cm
data.IDlvD = ttedata{PID, 'IDLVD'}; % measured in cm

data.Psv = 1.025*data.PrvD;      % not given in dataset, here it is calculated
data.Ppv = rhcdata{PID, 'PpcW'}; % using PpcW, which is a good measure of Ppv, this is a mean value

data.VlvM = mridata{PID, "VlvD"}; % using diastolic volume and assuming that it is end-diastolic volume
data.VlvmT = mridata{PID, "VlvS"}; % for teichholz calculations, need to make sure that we are using actual data points
data.Vlvm = data.VlvM-data.SV; 

data.VrvM = mridata{PID, "VrvD"}; % using diastolic volume and assuming that it is end-diastolic volume
data.VrvmT = mridata{PID, "VrvS"}; % for teichholz calculations
data.Vrvm = data.VrvM-data.SV; 

if teich == 1
    a = 7.3098; % from estimated parameters
    b = 14.6489;
    data.VlvM = a*data.IDlvD^2+b;
    data.Vlvm = a*data.IDlvS^2+b;
end

% If data.CO = Nan
if isnan(data.CO)
    data.SV = data.VlvM - data.VlvmT; % use LV Volume
    data.CO = data.HR*data.SV;
    data.Qtot = data.CO*1000/60;
    data.Vlvm = data.VlvM-data.SV;
    data.Vrvm = data.VrvM-data.SV;
end

% Check to make sure that volume is positive
if data.Vlvm < 0
    data.Vlvm = data.VlvmT;
end
if data.Vrvm < 0
    data.Vrvm = data.VrvmT;
end

% Unstressed volumes
data.VlvU = 0.10 * data.Vlvm; % scaled values
data.VrvU = 0.10 * data.Vrvm; % scaled values

data.TC  = 0.12; % Time for maximum contractility relative to the lenght of the heart beat
data.TR  = 0.18; % Time for relaxation relative to the lenght of the heart beat

data.EF = mridata{PID, 'EF'};

