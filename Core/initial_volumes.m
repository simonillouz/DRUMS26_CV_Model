% Init function will take all parameters defined in load_global and patient.m and output intial volume conditions

function [Init] = initial_volumes(pars, parsG, data)

% Define Parameters 
Vtot = data.TotalVol;


% Compliances (Calculated as 1/E for LV and RV)
Cpv = pars(7);
Csa	= pars(8);
Csv	= pars(9);
Cpa = pars(10);


Elvm = pars(12);

Ervm = pars(14);

Clv = 1/(Elvm);
Crv = 1/(Ervm);

% Unstressed Volumes

VpvU = parsG.VpvU;
VlvU = parsG.VlvU;
VsaU = parsG.VsaU;
VsvU = parsG.VsvU;
VrvU = parsG.VrvU;
VpaU = parsG.VpaU;

% NEW: Unstressed Pressures (Placeholders - Please update these!)
PlvU = 0;
PpaU = 0;
PpvU = 0;
PrvU = 0;
PsaU = 0;
PsvU = 0;

% % Calculate Denominator (D) 
% % The shared denominator across all the new equations
% D = Clv + Cpa + Cpv + Crv + Csa + Csv;
% 
% % --- 3. Calculate Specific Volumes ---
% Vlv = -(Clv*Cpa*PlvU - Clv*Cpa*PpaU + Clv*Cpv*PlvU - Clv*Cpv*PpvU + Clv*Crv*PlvU - Clv*Crv*PrvU + Clv*Csa*PlvU - Clv*Csa*PsaU + Clv*Csv*PlvU - Clv*Csv*PsvU + Clv*VpaU + Clv*VpvU + Clv*VrvU + Clv*VsaU + Clv*VsvU - Clv*Vtot - Cpa*VlvU - Cpv*VlvU - Crv*VlvU - Csa*VlvU - Csv*VlvU) / D;
% 
% Vpa = (Clv*Cpa*PlvU - Clv*Cpa*PpaU - Cpa*Cpv*PpaU + Cpa*Cpv*PpvU - Cpa*Crv*PpaU + Cpa*Crv*PrvU - Cpa*Csa*PpaU + Cpa*Csa*PsaU - Cpa*Csv*PpaU + Cpa*Csv*PsvU + Clv*VpaU - Cpa*VlvU - Cpa*VpvU - Cpa*VrvU - Cpa*VsaU - Cpa*VsvU + Cpa*Vtot + Cpv*VpaU + Crv*VpaU + Csa*VpaU + Csv*VpaU) / D;
% 
% Vrv = (Clv*Crv*PlvU - Clv*Crv*PrvU + Cpa*Crv*PpaU - Cpa*Crv*PrvU + Cpv*Crv*PpvU - Cpv*Crv*PrvU - Crv*Csa*PrvU + Crv*Csa*PsaU - Crv*Csv*PrvU + Crv*Csv*PsvU + Clv*VrvU + Cpa*VrvU + Cpv*VrvU - Crv*VlvU - Crv*VpaU - Crv*VpvU - Crv*VsaU - Crv*VsvU + Crv*Vtot + Csa*VrvU + Csv*VrvU) / D;
% 
% Vsa = (Clv*Csa*PlvU - Clv*Csa*PsaU + Cpa*Csa*PpaU - Cpa*Csa*PsaU + Cpv*Csa*PpvU - Cpv*Csa*PsaU + Crv*Csa*PrvU - Crv*Csa*PsaU - Csa*Csv*PsaU + Csa*Csv*PsvU + Clv*VsaU + Cpa*VsaU + Cpv*VsaU + Crv*VsaU - Csa*VlvU - Csa*VpaU - Csa*VpvU - Csa*VrvU - Csa*VsvU + Csa*Vtot + Csv*VsaU) / D;
% 
% Vsv = (Clv*Csv*PlvU - Clv*Csv*PsvU + Cpa*Csv*PpaU - Cpa*Csv*PsvU + Cpv*Csv*PpvU - Cpv*Csv*PsvU + Crv*Csv*PrvU - Crv*Csv*PsvU + Csa*Csv*PsaU - Csa*Csv*PsvU + Clv*VsvU + Cpa*VsvU + Cpv*VsvU + Crv*VsvU + Csa*VsvU - Csv*VlvU - Csv*VpaU - Csv*VpvU - Csv*VrvU - Csv*VsaU + Csv*Vtot) / D;

% Shared denominator
D = Clv + Cpa + Cpv + Crv + Csa + Csv;

% Solution Equations
Vlv = -(Clv*VpaU + Clv*VpvU + Clv*VrvU + Clv*VsaU + Clv*VsvU - Clv*Vtot - Cpa*VlvU - Cpv*VlvU - Crv*VlvU - Csa*VlvU - Csv*VlvU) / D;

Vpv = (Clv*VpvU + Cpa*VpvU - Cpv*VlvU - Cpv*VpaU - Cpv*VrvU - Cpv*VsaU - Cpv*VsvU + Cpv*Vtot + Crv*VpvU + Csa*VpvU + Csv*VpvU) / D;

Vrv = (Clv*VrvU + Cpa*VrvU + Cpv*VrvU - Crv*VlvU - Crv*VpaU - Crv*VpvU - Crv*VsaU - Crv*VsvU + Crv*Vtot + Csa*VrvU + Csv*VrvU) / D;

Vsa = (Clv*VsaU + Cpa*VsaU + Cpv*VsaU + Crv*VsaU - Csa*VlvU - Csa*VpaU - Csa*VpvU - Csa*VrvU - Csa*VsvU + Csa*Vtot + Csv*VsaU) / D;

Vsv = (Clv*VsvU + Cpa*VsvU + Cpv*VsvU + Crv*VsvU + Csa*VsvU - Csv*VlvU - Csv*VpaU - Csv*VpvU - Csv*VrvU - Csv*VsaU + Csv*Vtot) / D;

% Calculate the remaining volume for the Pulmonary Vein
Vpa = Vtot - Vlv - Vpv - Vrv - Vsa - Vsv;

Init = [Vpv Vlv Vsa Vsv Vrv Vpa];


% % Display the Results
% fprintf('Vlv = %.4f\n', Vlv);
% fprintf('Vpa = %.4f\n', Vpa);
% fprintf('Vpv = %.4f\n', Vpv);
% fprintf('Vrv = %.4f\n', Vrv);
% fprintf('Vsa = %.4f\n', Vsa);
% fprintf('Vsv = %.4f\n', Vsv);

end 
