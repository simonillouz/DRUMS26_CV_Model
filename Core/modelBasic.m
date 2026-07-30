function xdot = modelBasic(t,y,pars, parsG, tst)
% Right hand side of the differential equations (from original code)

% Define volumes
Vpv = y(1);
Vlv = y(2);
Vsa = y(3);
Vsv = y(4);
Vrv = y(5);
Vpa = y(6);

% Define unstressed volumes
VpvU = parsG.VpvU;
VlvU = parsG.VlvU;
VsaU = parsG.VsaU;
VsvU = parsG.VsvU;
VrvU = parsG.VrvU;
VpaU = parsG.VpaU;

% Define resistances
Rpul  = pars(1);
Rmval = pars(2);
Raval = pars(3);

Rsys   = pars(4);
Rtval  = pars(5);
Rpval  = pars(6);

% Define compliances - elastance for right heart and left heart
Cpv = pars(7);
Csa	= pars(8);
Csv	= pars(9);
Cpa = pars(10);

% Define elastances (heart parameters)
ElvM = pars(11);
Elvm = pars(12);
ErvM = pars(13);
Ervm = pars(14);

% Define heart timing
TC = pars(15); 
TR = pars(16);

% Left heart dynamics
Elv  = Elastance(t-tst,ElvM, Elvm, TC, TR);
Plv  = Elv*(Vlv-VlvU);

% Right heart dynamics
Erv = Elastance(t-tst, ErvM, Ervm, TC, TR);
Prv = Erv*(Vrv-VrvU);

Ppv = (Vpv-VpvU)/Cpv;
Psa = (Vsa-VsaU)/Csa;
Psv = (Vsv-VsvU)/Csv;
Ppa = (Vpa-VpaU)/Cpa;

% flow into LH (mitral valve)
if Ppv > Plv
    Qmval = (Ppv-Plv)/Rmval;
else
    Qmval = 0;
    
end
% flow out of LH (aortic valve)
if Plv > Psa
    Qaval = (Plv-Psa)/Raval;
else
    Qaval = 0;
end

% flow into RH (tricuspid valve)
if Psv > Prv
    Qtval = (Psv-Prv)/Rtval;
else
    Qtval = 0;
end
% flow out of RH (pulmonary valve)
if Prv > Ppa
    Qpval = (Prv-Ppa)/Rpval;
else
    Qpval = 0;
end

%Calculate flows
Qpul = (Ppa - Ppv)/Rpul;
Qsys = (Psa - Psv)/Rsys;

% Differential Equations
dVpv = Qpul - Qmval;
dVlv = Qmval - Qaval;
dVsa = Qaval - Qsys;
dVsv = Qsys - Qtval;
dVrv = Qtval - Qpval;
dVpa = Qpval - Qpul;

% Right hand side of differential equations
xdot = [dVpv; dVlv; dVsa; dVsv; dVrv; dVpa];
