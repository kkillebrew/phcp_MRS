%% get the paths
%%% add our custom matlab code:
addpath(genpath('/home/shaw-raid1/matlab_tools/MRS_scripts/'))

%%% add romain's code from GitHub
addpath(genpath('/home/shaw-raid1/matlab_tools/matspec/'))

%%% add our own custom functions
addpath(genpath('/home/shaw-raid1/matlab_tools/matspec_custom_phcp/'),'-begin')

%%% add only spect processing
% addpath(genpath('/home/range4-raid1/gosia/matlab/from_Romain_151126/matspec/spect_processing'))

%%% add Gosia's older copy of the toolbox
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/anat_processing/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/spect_processing/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/readdicom/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/lcmodel/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/plotting/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/tools/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/lana_fit/',path);
% path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/hsvd_processing/',path);


if datenum(version('-date')) > datenum('February 11, 2014')
    error(['This code relies on SPM8, and will NOT function with versions '...
        'of Matlab beyond 2014a!!'])
end
%% choose between OCC & PFC...
which_ROI = questdlg('Which ROI do you want to analyze?','Which ROI?','OCC','PFC','OCC');
if strcmp(which_ROI,'OCC')
    ROI_idx = 1;
elseif strcmp(which_ROI,'PFC')
    ROI_idx = 2;
end
%% set processing parameters

% setting up processing parameters
par=processing_spec;
par.mean_line_broadening=4;

% params for relaxation, water correction, from Gosia 20190610
par.gmT1 = 2132;
par.gmT2 = 50; 
par.wmT1 = 1220;
par.wmT2 = 55;
par.csfT1 = 4425;
par.csfT2 = 141;
par.TR = 5000;
par.TE = 8;
par.water_fraction = [0.80 0.71 0.97];

% show figures at all? 1 = yes, 0 = no
par.figure = 1;

% close figures after saving them? 1 = yes, 0 = no
close_all = 1;

%% for subjects with an issue with w1, use a different water reference
subj_with_w1_issues = {'P6004071_20180801/MR-SE028-OCC_steam_eja_w1/'};

alternative_h2o_ref = {'P6004071_20180801/MR-SE030-OCC_steam_eja_w4_RES/'};

%% Throw out known bad data points for particular scans

% OCC first, ROI_idx = 1
subj_with_data_to_toss{1} = {'P6010622_20180908','P2110417_20181029','P1010213_20180211','P6004213_20190116','P6004663_20180911','P6010465_20180825'};

%TRs_to_toss{1} = {[22:24 28:30],[78],[59],[24],[86],[58]};
TRs_to_toss{1} = {[2,4,8,23,29,32,36,49,56,62,79],[78],[59],[24],[86],[58]};

% PFC 2nd, ROI_idx = 2
subj_with_data_to_toss{2} = {'P6011098_20190507'};

TRs_to_toss{2} = {[42]};

%% find metab spectra
d=get_subdir_regex('/home/shaw-raid1/data/MRS/dicom_data','P\d\d\d\d\d\d\d*',...
    ['MR-SE\d\d\d-' which_ROI '_steam_eja_metab']);
toss_me = [];
for iD = 1:numel(d)
    if isempty(strfind(d{iD},'RES'))
        toss_me = [toss_me iD];
    end
end
d(toss_me) = [];

%% find water spectra
dw=get_subdir_regex('/home/shaw-raid1/data/MRS/dicom_data','P\d\d\d\d\d\d\d*',...
    ['MR-SE\d\d\d-' which_ROI '_steam_eja_w1']);
toss_me = [];
issue_list = [];
for iD = 1:numel(dw)
    if ~isempty(strfind(dw{iD},'noOVS'))
        toss_me = [toss_me iD];
    end
    find_issue_idx = []; % clear this
    for iIssue = 1:numel(subj_with_w1_issues)
        find_issue_idx = strfind(dw{iD},subj_with_w1_issues{iIssue});
        if ~isempty(find_issue_idx)
            dw{iD} = fullfile(dw{iD}(1:find_issue_idx-1),...
                alternative_h2o_ref{iIssue});
        end
    end
    if ~isempty(find_issue_idx)
        issue_list(iD) = 1;
    else issue_list(iD) = 0;
    end
end
dw(toss_me) = [];
issue_list(toss_me) = [];
issue_idx = find(issue_list);

%% reading in the data and organizing 
f=explore_spectro_data(char(d));

warndlg('Crudely tossing any scan that has the same name and date as scan after it, assuming this is an unusable data set that was repeated and should be tossed...');
toss_me = [];
for iF = 1:(numel(f)-1) % can't include last one...
    if strcmp(f(iF).sujet_name , f(iF+1).sujet_name)
        warning(['tossing suspected duplicate metab scan for ' f(iF).sujet_name])
        toss_me = [toss_me iF];
    end
end
f(toss_me) = [];

f=concatenate_fid(f);


fw=explore_spectro_data(char(dw));

toss_me = [];
for iF = 1:(numel(fw)-1) % can't include last one...
    if strcmp(fw(iF).sujet_name , fw(iF+1).sujet_name)
        warning(['tossing suspected duplicate w1 scan for ' fw(iF).sujet_name])
        toss_me = [toss_me iF];
    end
    if sum(issue_idx == iF)
        fw(iF).fid = fw(iF).fid(:,1);
        % only use the first of the w4 series
    end
end
fw(toss_me) = [];

% check that subjects in f & fw match...
doesnt_match = [];
for iSubj = 1:min([numel(f) numel(fw)])
    if ~strcmp(f(iSubj).sujet_name, fw(iSubj).sujet_name)
        doesnt_match = [doesnt_match ; iSubj];
    end
end
if ~isempty(doesnt_match)
    error(['Some scan names don''t match, start by looking at ' ...
        f(doesnt_match(1)).sujet_name]);
end

% initial processing
% eddy current compensation
f_ecc=eddycorrection(f,fw);
fw_ecc=eddycorrection(fw,fw);

% remove first point and add extra 0 at the end (to have the same number of
% points as initially
fo=remove_first_points_fillup(f_ecc,1);

%% toss known bad TRs
% mps 20190111 - adding the ability to toss TRs...
% to check for bad TRs: mps_check_bad_avgs(fo, subj_idx, par)
% where subj_idx is the # representing where in the list of subjects this person is (the Nth subject)
% use the data cursor to select the bad data in the plot, the z value is
% the TR # to exclude. par is optional, include if you want to run line
% broadening! Can also be used for other (e.g., processed) data sets, such
% as fc3_3, rather than fo. For example: mps_check_bad_avgs(fc3_3, 83, par)

n_toss = 0;
for iSubj = 1:numel(fo)
    find_toss_idx = find(strcmp(fo(iSubj).sujet_name,subj_with_data_to_toss{ROI_idx})); % find the subject in the toss list
    if ~isempty(find_toss_idx)
        fo(iSubj).fid(:,TRs_to_toss{ROI_idx}{find_toss_idx}) = []; % find which TRs to toss for this person, and then toss them
        fo(iSubj).Nex = fo(iSubj).Nex - numel(TRs_to_toss{ROI_idx}{find_toss_idx}); % update these numbers... I hope this is right!!
        fo(iSubj).Number_of_spec = fo(iSubj).Number_of_spec - numel(TRs_to_toss{ROI_idx}{find_toss_idx});
        n_toss = n_toss +1;
        warning(['Tossing TRs ' num2str(TRs_to_toss{ROI_idx}{find_toss_idx}) ...
            ' for ' fo(iSubj).sujet_name ', as requested...']);
    end
end
if n_toss > 0
    warndlg('Crudely tossing some TRs, as specified!')
end

%% freq & phase correction - option #1 = Corr
% par.method='correlation';
% par.correlation_bound = [];
% fc=processing_spec(fo,par);
% 
% % par.same_fig = 1; % turn this on to plot on the same figure;
% % par.arg_give_me_a_new_one = 1; % choose figure number
% for iSubj = 1:numel(fo)
% plot_spectrum(fo(iSubj),par) % before processing
% get_ax = axis;
% axis([0.5 2.5 get_ax(3) get_ax(4)]);
% title([fo(iSubj).sujet_name ' Before'])
% set(gcf,'color','w')
% end
% 
% for iSubj = 1:numel(fo)
% plot_spectrum(fc(iSubj),par) % after processing
% get_ax = axis;
% axis([0.5 2.5 get_ax(3) get_ax(4)]);
% title([fc(iSubj).sujet_name ' 1) Corr.'])
% set(gcf,'color','w')
% end

% % to know if frequency and phase correction were successful, look for
% % variability in the peak freq of NAA (~2ppm - for freq) and variability in
% % the shape (e.g., weird dips) on the flanks of the NAA peak (~1.95ppm -
% % for phase).

%% option #2 = corr w/ bounds
% 
% par.correlation_bound = [1.8 3.5];
% par.method='correlation';
% fc2=processing_spec(fo,par);
% 
% for iSubj = 1:numel(fo)
% plot_spectrum(fc2(iSubj),par) % after processing
% get_ax = axis;
% axis([0 6 get_ax(3) get_ax(4)]);
% title([fo(iSubj).sujet_name ' 2) Corr. w/ bound'])
% set(gcf,'color','w')
% end

%% option #3 = max + abs

par.method='max';
par.correct_freq_mod='abs';
fc3_1=processing_spec(fo,par);
if close_all
    close all
end
par.correct_freq_mod='abs';
fc3_2=processing_spec(fc3_1,par);
if close_all
    close all
end
par.correct_freq_mod='real';
fc3_3=processing_spec(fc3_2,par);
if close_all
    close all
end

plot_dir = fullfile('/home/shaw-raid1/data/MRS/saved_plots',which_ROI);
plot_dir_contents = dir([plot_dir '/*.fig']);
if ~isempty(plot_dir_contents)
    clear ask_wipe
    ask_wipe = ['Saved plots folder is not empty (' plot_dir '), do you want to delete all contents and re-write, or quit?'];
    wipe_dir = questdlg(ask_wipe,'Do you want to delete?','Delete','Quit','Quit');
    if strcmp(wipe_dir,'Delete')
        warning(['Deleting contents of ' plot_dir ', as requested...']);
        eval(['! rm -f ' plot_dir '/*.fig']);
    elseif strcmp(wipe_dir,'Quit')
        return
    end
end

if par.figure % if we are plotting figures...
    for iSubj = 1:numel(fo)
        plot_spectrum(fo(iSubj),par) % before processing
        get_ax = axis;
        axis([0 6 get_ax(3) get_ax(4)]);
        title([fo(iSubj).sujet_name ' uncorrected'])
        set(gcf,'color','w')
        savefig(fullfile(plot_dir,[fo(iSubj).sujet_name '_' which_ROI '_uncorrected.fig']))
        if close_all
            close(gcf);
        end
        
        plot_spectrum(fc3_3(iSubj),par) % after processing
        get_ax = axis;
        axis([0 6 get_ax(3) get_ax(4)]);
        title([fo(iSubj).sujet_name ' 3) Max + Abs. x2 & Real'])
        set(gcf,'color','w')
        savefig(fullfile(plot_dir,[fo(iSubj).sujet_name '_' which_ROI '_corrected.fig']))
        if close_all
            close(gcf);
        end
        
        
        SD_NAA(iSubj) = mps_SD_NAA(fc3_3(iSubj),par);
    end
end

%% plot all results overlapping
% for iSubj = 1:numel(fo)
% figure; hold on
% scale_ppm = mps_scale_ppm(fo(iSubj)); % no difference if I use fc, so just calculate this once
% plot_orig = mean(real(fftshift(fft(fo(iSubj).fid),1)),2); % stole from correct_freq_and_phase_by_correlation.m
% plot(scale_ppm,plot_orig,'k-')
% plot_opt1 = mean(real(fftshift(fft(fc(iSubj).fid),1)),2);
% plot(scale_ppm,plot_opt1,'r-')
% plot_opt2 = mean(real(fftshift(fft(fc2(iSubj).fid),1)),2);
% plot(scale_ppm,plot_opt2,'b-')
% plot_opt3 = mean(real(fftshift(fft(fc3_3(iSubj).fid),1)),2);
% plot(scale_ppm,plot_opt3,'g-')
% 
% set(gca,'Xdir','reverse','fontsize',12);
% get_ax = axis;
% axis([0 6 get_ax(3) get_ax(4)])
% legend('orig','corr','corr w/ bound','max + abs. x2 & real')
% set(gcf,'color','w')
% xlabel('freq. (ppm)')
% title(fo(iSubj).sujet_name)
% end
%% to obtain water linewidth
fw_ecc=get_water_width(fw_ecc, par.figure);
fw_ecc.water_width_no_cor % gosia does not trust the correction, so don't use it
% (already have done eddy current compensation) - will be a bit smaller
% than the linewidth of water we measured at the scanner
%% preparing for LCModel analysis:

pp.root_dir=['/home/shaw-raid1/data/MRS/processed_data/' which_ROI];
% it uses fc.sujet_name and fc.SerDescr for naming RAW files
pp.gessfile=5; % different options for setting file names based on dicom info

% remove exisiting proc files
proc_dir_contents = dir(fullfile(pp.root_dir,'*.RAW'));
if ~isempty(proc_dir_contents)
    ask_wipe = ['Processed data folder is not empty (' pp.root_dir '), do you want to delete all contents and re-write, or quit?'];
    wipe_dir = questdlg(ask_wipe,'Do you want to delete?','Delete','Quit','Quit');
    if strcmp(wipe_dir,'Delete')
        warning(['Deleting contents of ' pp.root_dir ', as requested...']);
        eval(['! rm ' pp.root_dir '/*.PLOTIN']);
        eval(['! rm ' pp.root_dir '/*.RAW']);
        eval(['! rm ' pp.root_dir '/*.H2O']);
    elseif strcmp(wipe_dir,'Quit')
        return
    end
end

fc = fc3_3; % use Abs + max, Abs + max, Real, as above...
write_fid_to_lcRAW(fc,pp,fw_ecc);

% remove existing lcm files
lcm_fit_dir = fullfile('/home/shaw-raid1/data/MRS/processed_data/', which_ROI,'lcm_fit');
lcm_fit_contents = dir(lcm_fit_dir);
if ~isempty(lcm_fit_contents)
    clear ask_wipe
    ask_wipe = ['LCM folder is not empty (' lcm_fit_dir '), do you want to delete all contents and re-write, or quit?'];
    wipe_dir = questdlg(ask_wipe,'Do you want to delete?','Delete','Quit','Quit');
    if strcmp(wipe_dir,'Delete')
        warning(['Deleting contents of ' lcm_fit_dir ', as requested...']);
        eval(['! rm -f ' lcm_fit_dir '/*.CONTROL']);
        eval(['! rm -f ' lcm_fit_dir '/*.COORD']);
        eval(['! rm -f ' lcm_fit_dir '/*.PRINT']);
        eval(['! rm -f ' lcm_fit_dir '/*.PS']);
    elseif strcmp(wipe_dir,'Quit')
        return
    end
end
%% run the LCModel analysis

[foo, server_name] = system('hostname');
if strcmp(server_name(1:16),'kea.cmrr.umn.edu') && datenum(version('-date')) <= datenum('February 12, 2009')
    % processing_LCmodel('steam_phcp','lcm_fit') % 20190208 mps testing new
    % control file params to fit lipids
    processing_LCmodel('steam_phcp_lipids','lcm_fit')
    % n.b. matlab will only run on kea if you select a 32 bit installation,
    % such as 2009a...
    
else
    error('This command needs to be run on the server kea.cmrr.umn.edu -- ssh over there and try again. Be sure to use a 32-bit version of Matlab, such as 2009a...');
end
%% convert the data files to .pdf and .csv
% 1. generate pdf file for review of the fits
% be sure to select the folder you wrote the LCmodel data to in the step
% above, e.g., lcm_fit

if ~exist('proc_dir','var')
    proc_dir = '/home/shaw-raid1/data/MRS/processed_data';
end

[foo, server_name] = system('hostname');
if strcmp(server_name(1:20),'range4.cmrr.umn.edu')
    
    dir_res=get_subdir_regex;
    concat_ps_file(dir_res)
else
    error('This command needs to be run on the server range4.cmrr.umn.edu -- ssh over there and try again...');
end

% 2. generate spreadsheet

% fix file permissions for the files you just wrote here...
eval(['! chmod g+w -R ' dir_res{1}]);

% check if we know which ROI to use
if numel(dir_res) > 1
    error('You selected more than 1 lcm folder...')
end
if ~exist('which_ROI','var')
    occ_idx = regexp(dir_res{1},'OCC');
    pfc_idx = regexp(dir_res{1},'PFC');
    if isempty(occ_idx) && isempty(pfc_idx)
        error('which_ROI is not defined, and was not detected automatically from the lcm folder name!')
    elseif ~isempty(occ_idx) && ~isempty(pfc_idx)
        error('which_ROI is not defined, and the lcm folder name has BOTH OCC & PFC in it??')
    elseif ~isempty(occ_idx) && isempty(pfc_idx)
        which_ROI = 'OCC';
        warning(['automatically detected which ROI (' which_ROI ') based on lcm folder name']);
    elseif isempty(occ_idx) && ~isempty(pfc_idx)
        which_ROI = 'PFC';
        warning(['automatically detected which ROI (' which_ROI ') based on lcm folder name']);
    end
end

c=get_result(dir_res);
write_name = fullfile(proc_dir,...
    [datestr(now,'yyyymmdd') '_phcp_' which_ROI '_' num2str(numel(c.suj)) 'subj_unscaled.csv']);
if ~exist(write_name,'file')
    write_conc_res_to_csv(c,write_name)
else
    error(['The following file already exists, to create a new version, you must delete the existing one! - ' write_name]);
end

% scale to tCr
c_cr = correct_result_ration(c,'','Cr_PCr',8);
write_name = fullfile(proc_dir,...
    [ datestr(now,'yyyymmdd') '_phcp_' which_ROI '_' num2str(numel(c.suj)) 'subj_Cr_scaled.csv']);
if ~exist(write_name,'file')
    write_conc_res_to_csv(c_cr,write_name)
else
    error(['The following file already exists, to create a new version, you must delete the existing one! - ' write_name]);
end

% scale to Water
% first get all tissue fractions
tissue_fract_data = mps_run_all_MRS_tissue_correct;

not_missing_subjs = ~isnan( mean( tissue_fract_data.tissue_fract(...
    :, ROI_idx, :), 3));

cw = c;
if sum(not_missing_subjs) ~= numel(c.suj)
    error('number of subjects in tissue_fract_data does not match data in c!');
else
    cw.fgray = tissue_fract_data.tissue_fract(not_missing_subjs, ROI_idx, 1); % GM first
    cw.fwhite = tissue_fract_data.tissue_fract(not_missing_subjs, ROI_idx, 2); % WM second
    cw.fcsf = tissue_fract_data.tissue_fract(not_missing_subjs, ROI_idx, 3); % CSF third
end

corr_factor = water_content_attenuation(cw, par);
c_wscale = correct_result(cw, corr_factor);

write_name = fullfile(proc_dir,...
    [ datestr(now,'yyyymmdd') '_phcp_' which_ROI '_' num2str(numel(c.suj)) 'subj_H2O_scaled.csv']);
if ~exist(write_name,'file')
    write_conc_res_to_csv(c_wscale,write_name)
else
    error(['The following file already exists, to create a new version, you must delete the existing one! - ' write_name]);
end

% rename .pdf and .ps

write_name = fullfile(proc_dir,...
    [datestr(now,'yyyymmdd') '_phcp_' which_ROI '_' num2str(numel(c.suj)) 'subj.pdf']);
run_cmd = ['! mv ' proc_dir '/' which_ROI '_lcm_fit.pdf ' write_name ' ;'];
eval(run_cmd);

write_name = fullfile(proc_dir,...
    [datestr(now,'yyyymmdd') '_phcp_' which_ROI '_' num2str(numel(c.suj)) 'subj.ps']);
run_cmd = ['! mv ' proc_dir '/' which_ROI '_lcm_fit.ps ' write_name ' ;'];
eval(run_cmd);
