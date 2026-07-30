function DriverBasic_opt(PID, teich)
% Name: DriverBasic_opt
% Inputs: takes in patient ID and flags 1 for using calculated volumes 
% (LV), 0 for measured volumes (LV)
% Requires:
% Function

close all

cd Core
data = Patient(PID, teich); % Read data from patient (PID)

% Load nominal parameter values and initial conditions, upper and lower
% bounds for the parameter estimation
[pars, parsG, times, Init, low, hi] = load_global(PID,data);

% Save all global variables to the structure data
% Note no globals for parallelization
data.parsG    = parsG;
data.Init     = Init;
data.dt       = times.dt;
data.T        = times.T;
data.CC       = 5;
data.ALLPARS  = pars;
data.INDMAP   = [1 4 8:14]; % I was able to add 6 or 7, adding both make 6 hit lower bound
%data.INDMAP   = [1 4 8 9 11:14]; % I was able to add 6 or 7, adding both make 6 hit lower bound
data.hi       = hi;
data.low      = low;

optx   = pars(data.INDMAP);
opthi  = data.hi(data.INDMAP);
optlow = data.low(data.INDMAP);
sprintf('optx %d \n opthi %d \n optlow %d \n',optx,opthi,optlow);
if teich == 0
    folderName = fullfile('../MultiOptPatients', sprintf('Patient_%d', PID));
elseif teich == 1
    % folderName = fullfile('../MultiOptPatientsTeich', sprintf('Patient_%d', PID));
    folderName = fullfile('../MultiOptPatientsCalc', sprintf('Patient_%d', PID)); % so that we have a different file name
end

if ~exist(folderName, 'dir')
    mkdir(folderName);
end
    
data.numMSOpt = 15;
maxiter  = 40; 
mode     = 2;
nu0      = 2e-1;

% Matrix with random values from -1 to 1, to ensure that the random numbers
% are different for each simulation, we use "shuffle" as the seed
rng('shuffle');
rt           = rand(data.numMSOpt+1, length(optx)) .* 2 - 1;
rt(1,:)      = 0; %first pass is nominal with 0 noise
parsNoiseMat = (1 + 0.05 .* rt) .* optx(:)';  % (numMSOpt+1) x nParams

parfor i = 1:data.numMSOpt+1
    pars_noise_i = parsNoiseMat(i, :)';
    [xopt_i, histout_i, ~,~,~,~,~] = ...
        newlsq_v2(pars_noise_i, 'opt_wrap', 1e-4, maxiter, mode, nu0, ...
                  opthi, optlow, data);
    nomPars              = data.ALLPARS;
    nomPars(data.INDMAP) = pars_noise_i;
    optPars              = data.ALLPARS;
    optPars(data.INDMAP) = xopt_i;
    
    s = struct("histout",histout_i,"nomPars",nomPars,"optPars",optPars,"data",data);
    if teich == 0
        save(sprintf('../MultiOptPatients/Patient_%d/Opt_Run%d.mat',PID,i-1),"-fromstruct",s);
    elseif teich == 1
        save(sprintf('../MultiOptPatientsTeich/Patient_%d/Opt_Run%d.mat',PID,i-1),"-fromstruct",s);
    end
    
end

cd ..

end % Function