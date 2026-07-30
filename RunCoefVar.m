function RunCoefVar(PIDs)

cd MultiOptPatients/
missing  = [];
for i = 1:343
    foldName = strcat('Patient_', num2str(i));
    if ~exist(foldName,"dir")
        missing = [missing i];
    else
        if numel(dir(foldName)) <= 2
            missing = [missing i];
        end
    end
end
cd ../Core %to RunModel/Core

for i = PIDs
    if ~(ismember(i, missing))
        coefVariation(i)
        cd Core
    end
end

cd .. %to RunModel where this file lives