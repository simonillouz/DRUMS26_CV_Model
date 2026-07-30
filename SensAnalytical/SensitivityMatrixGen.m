%% 1. Define State Variables
syms Vpv Vlv Vsa Vsv Vrv Vpa
X = [Vpv; Vlv; Vsa; Vsv; Vrv; Vpa];

%% 2. Define the 16 Parameters (Columns of the sensitivity matrix)
syms Rpul Rmval Raval Rsys Rtval Rpval Cpv Csa Csv Cpa ElvM Elvm ErvM Ervm TC TR
pars_sym = [Rpul, Rmval, Raval, Rsys, Rtval, Rpval, Cpv, Csa, Csv, Cpa, ElvM, Elvm, ErvM, Ervm, TC, TR];

% Define the Unstressed Volumes (grouped together to mirror parsG structure)
syms VpvU VlvU VsaU VsvU VrvU VpaU
parsG_vector = [VpvU, VlvU, VsaU, VsvU, VrvU, VpaU];

% Time-dependent heart elastances are treated as localized inputs during calculation
syms Elv Erv 

%% 3. Recreate the Cardiovascular Model Equations (Matching JacobianDisplay.m)
Plv = Elv * (Vlv - VlvU);
Prv = Erv * (Vrv - VrvU);
Ppv = (Vpv - VpvU) / Cpv;
Psa = (Vsa - VsaU) / Csa;
Psv = (Vsv - VsvU) / Csv;
Ppa = (Vpa - VpaU) / Cpa;

% Valve flows using symbolic piecewise conditions (mimics max(..., 0))
Qmval = piecewise(Ppv > Plv, (Ppv - Plv) / Rmval, 0);
Qaval = piecewise(Plv > Psa, (Plv - Psa) / Raval, 0);
Qtval = piecewise(Psv > Prv, (Psv - Prv) / Rtval, 0);
Qpval = piecewise(Prv > Ppa, (Prv - Ppa) / Rpval, 0);

% Linear structural vessel flows
Qpul = (Ppa - Ppv) / Rpul;
Qsys = (Psa - Psv) / Rsys;

% Differential Equations System (The state derivatives dX/dt)
dVpv = Qpul - Qmval;
dVlv = Qmval - Qaval;
dVsa = Qaval - Qsys;
dVsv = Qsys - Qtval;
dVrv = Qtval - Qpval;
dVpa = Qpval - Qpul;

% Differential equation vector (Rows of the sensitivity matrix)
F = [dVpv; dVlv; dVsa; dVsv; dVrv; dVpa];

%% 4. Compute the Core Parameter Sensitivity Matrix
% 1. Derivatives with respect to the first 10 structural circuit parameters
S_circuit = jacobian(F, pars_sym(1:10)); 

% 2. Derivatives with respect to Elv and Erv themselves (The Chain Rule Links)
S_Elv = jacobian(F, Elv);
S_Erv = jacobian(F, Erv);

% Combine them into a 6 row x 12 column matrix
% Columns 1-10: Structural parameters
% Column 11: dF/dElv
% Column 12: dF/dErv
S_combined = [S_circuit, S_Elv, S_Erv];

disp(size(S_combined))

%% 5. Compile and Generate the Function File
% Note: S_combined is now 6x12. We exclude the 6 elastance/time parameters from vars
% because we will compute their chain rule components dynamically in the driver.
matlabFunction(S_combined, 'File', 'get_SensitivityMatrix', 'Vars', {X, pars_sym(1:10), parsG_vector, Elv, Erv});

disp('==================================================================')
disp('Success! An updated "get_SensitivityMatrix.m" (6x12) has been generated.')
disp('==================================================================')
