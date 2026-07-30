function RunDB_sens(func_name, PID, N_samples, Target_output)
% Function name: GSA or LSA
% PID: takes in patient ID number
% N_samples: the number of samples needed for GSA (if doing LSA, user can
    % set as 0)
% Target_output: the state of interest ('Ppv','Plv','Psa','Psv','Prv','Ppa,
    % 'Vpv','Vlv','Vsa','Vsv','Vrv','Vpa', 'Qpul','Qmval','Qaval','Qtval',
    % 'Qpval','Qsys', 'VlvM','Vlvm','VrvM','Vrvm', 'CO','SV')
    % (if doing LSA, user can set at 0)

switch func_name
    case 'GSA'
        cd GSA
        run_cardio_GSA(PID, N_samples, Target_output)
    case 'LSA'
        cd SensAnalytical/
        SensitivityAnalyze2(PID)
end
end