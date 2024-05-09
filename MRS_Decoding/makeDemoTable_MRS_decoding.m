% Make demographics table for MRS subjects
% KWK - 20210421

% Current dir
options.curDur = pwd;

% Load in subject number list
datafull = load('mrsDecode_20210421.mat','subj_date');
fileNames = {datafull.subj_date};
options.subj_number = fileNames{1}(:,1);
options.date_number = fileNames{1}(:,2);
options.partGroup = zeros([length(options.subj_number) 1]);
options.partGroup(options.subj_number < 2000000) = 1;   % Controls
options.partGroup(options.subj_number >= 2000000 &...
    options.subj_number < 6000000) = 2;   % Relatives
options.partGroup(options.subj_number >= 6000000) = 3;   % Patients
options.subj_group_def = 1;

% Call MPS phcp_demographics function
cd /home/shaw-raid1/matlab_tools/mpsCode/
options.demoData = phcp_demographics(options);
cd(options.curDur);

% Make table of demo data for the 3 groups








