function output = single_page_phcpNoiseQuant_func(options)
% useage: output = single_page_phcpNoiseQuant_func(options)
%
% mps 20190205
%% set options
if ~exist('options','var')      
    options = [];
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 0; % 1 = yes, 0 = no
end
if ~isfield(options,'add_paths')
    options.add_paths = 'default'; % default, only_spect_proc, gosia_old
end
% only doing OCC
options.which_ROI = 'OCC';
ROI_idx = 1;

% close figures after saving them? 1 = yes, 0 = no
if ~isfield(options,'close_all')
    options.close_all = 1;
end

% set frequency & phase correction options
if ~isfield(options,'freq_phase_corr')
    options.freq_phase_corr = '(max + abs) x2 + real';
    % valid options are:
    %   (max + abs) x2 + real
    %   corr
    %   corr w/ bounds
end

% set option for which part of the process you want to run
if ~isfield(options,'which_step')
    options.which_step = questdlg('Which processing step do you want to perform?',...
        'Which step?','pre-LCM','LCM','post-LCM','pre-LCM');
end
if ~sum(strcmp(lower(options.which_step), ...
        lower({'pre-LCM','LCM','post-LCM'}) ))
    error(['Unrecognized value for options.which_step: ' options.which_step]);
end

% set option for how to exclude bad data -- manually vs. automatic



output = [];

%% get the paths
default_strs = {'default', 'standard', 'regular', 'defaults'};
spect_proc_strs = {'only_spect_proc', 'only spect proc', 'spect_proc',...
    'spect proc', 'spec_proc', 'spec proc', 'spect_proc_only', 'spect proc only'};
gosia_old_strs = {'gosia_old','gosias_old','gosia old','gosias old',...
    'old_gosia','old gosia','old_gosias','old gosias'};

if sum(strcmp(lower(options.add_paths), default_strs))
    %%% add our custom matlab code:
    addpath(genpath('/home/shaw-raid1/matlab_tools/MRS_scripts/'))
    
    %%% add romain's code from GitHub
    addpath(genpath('/home/shaw-raid1/matlab_tools/matspec/'))
    
    %%% add our own custom functions
    addpath(genpath('/home/shaw-raid1/matlab_tools/matspec_custom_phcp/'),'-begin')
    
elseif sum(strcmp(lower(options.add_paths), spect_proc_strs))
    %%% add only spect processing
    addpath(genpath('/home/range4-raid1/gosia/matlab/from_Romain_151126/matspec/spect_processing'))
    
elseif sum(strcmp(lower(options.add_paths), gosia_old_strs))
    %%% add Gosia's older copy of the toolbox
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/anat_processing/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/spect_processing/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/readdicom/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/lcmodel/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/plotting/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/tools/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/lana_fit/',path);
    path('/home/orochi4-raid2/gosia/matlab/from_Romain_151126/matspec/hsvd_processing/',path);
else
    error(['Value for options.add_path not recognized: ' options.add_paths]);
end

%% set processing parameters
if ~isfield(options,'par')
    par = processing_spec;
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
    
end

if ~isfield(par,'figure')
    % show figures at all? 1 = yes, 0 = no
    if options.displayFigs
        par.figure = 1;
    else
        par.figure = 0;
    end
end

%% now split into 3 processing steps

if strcmp(lower(options.which_step), 'pre-lcm')
%% check matlab version
if datenum(version('-date')) > datenum('February 11, 2014')
    error(['This code relies on SPM8, and will NOT function with versions '...
        'of Matlab beyond 2014a!!'])
end

%% for subjects with an issue with w1, use a different water reference
% subj_with_w1_issues = {'P6004071_20180801/MR-SE028-OCC_steam_eja_w1/'};
% 
% alternative_h2o_ref = {'P6004071_20180801/MR-SE030-OCC_steam_eja_w4_RES/'};
% already dealt with this 1 subj during dicom linking

%% Throw out known bad data points for particular scans

% % OCC first, ROI_idx = 1
% subj_with_data_to_toss{1} = {'P6010622_20180908','P2110417_20181029',...
%     'P1010213_20180211','P6004213_20190116','P6004663_20180911',...
%     'P6010465_20180825','P1010309_20191221'};
% 
% TRs_to_toss{1} = {[2,4,8,23,29,32,36,49,56,62,79],[78],...
%     [59],[24],[86],...
%     [58],[34 35]};
% 
% % PFC 2nd, ROI_idx = 2
% subj_with_data_to_toss{2} = {'P6011098_20190507'};
% 
% TRs_to_toss{2} = {[42]};

% already tossed data during dicom linking

%% find metab spectra
d1=get_subdir_regex('/home/shaw-raid1/data/MRS/processed_data/noise_quant/dicom_shortcuts',...
    'P1\d\d\d\d\d\d.*chunk.*');
d2=get_subdir_regex('/home/shaw-raid1/data/MRS/processed_data/noise_quant/dicom_shortcuts',...
    'P6\d\d\d\d\d\d.*chunk.*');
d = cat(2,d1,d2); % use only patients & controls

%% find water spectra
dw1=get_subdir_regex('/home/shaw-raid1/data/MRS/processed_data/noise_quant/dicom_shortcuts',...
    'P1\d\d\d\d\d\d.*_water');
dw2=get_subdir_regex('/home/shaw-raid1/data/MRS/processed_data/noise_quant/dicom_shortcuts',...
    'P6\d\d\d\d\d\d.*_water');
dw = cat(2,dw1,dw2); % use only patients & controls

dw = repmat(dw,17,1);

dw = dw(:)'; % reshape like this to use same water ref for each of 17 chunks

%% reading in the data and organizing 
% metab data
f=explore_spectro_data(char(d));

% f=concatenate_fid(f); don't want to concat, that's the whole point...

% water data
fw=explore_spectro_data(char(dw));

for iSubj = 1:min([numel(f) numel(fw)])
    f(iSubj).sujet_name = f(iSubj).ser_dir(1:17);
    fw(iSubj).sujet_name = fw(iSubj).ser_dir(1:17);
end


% check that subjects in f & fw match...
output.subj_date = [];
doesnt_match = [];
for iSubj = 1:min([numel(f) numel(fw)])
    if ~strcmp(f(iSubj).sujet_name, fw(iSubj).sujet_name)
        doesnt_match = [doesnt_match ; iSubj];
    end
    
    output.subj_date{iSubj,1} = f(iSubj).sujet_name;
end
if ~isempty(doesnt_match)
    error(['Some scan names don''t match, start by looking at ' ...
        f(doesnt_match(1)).sujet_name]);
end

% initial processing
% eddy current compensation
f_ecc = eddycorrection(f, fw);
fw_ecc = eddycorrection(fw, fw);

% remove first point and add extra 0 at the end (to have the same number of
% points as initially
fo = remove_first_points_fillup(f_ecc, 1);

%% toss known bad TRs
% mps 20190111 - adding the ability to toss TRs...
% to check for bad TRs: mps_check_bad_avgs(fo, subj_idx, par)
% where subj_idx is the # representing where in the list of subjects this person is (the Nth subject)
% use the data cursor to select the bad data in the plot, the z value is
% the TR # to exclude. par is optional, include if you want to run line
% broadening! Can also be used for other (e.g., processed) data sets, such
% as fc, rather than fo. For example: mps_check_bad_avgs(fc, 83, par)

% n_toss = 0;
% for iSubj = 1:numel(fo)
%     find_toss_idx = find(strcmp(fo(iSubj).sujet_name,subj_with_data_to_toss{ROI_idx})); % find the subject in the toss list
%     if ~isempty(find_toss_idx)
%         fo(iSubj).fid(:,TRs_to_toss{ROI_idx}{find_toss_idx}) = []; % find which TRs to toss for this person, and then toss them
%         fo(iSubj).Nex = fo(iSubj).Nex - numel(TRs_to_toss{ROI_idx}{find_toss_idx}); % update these numbers... I hope this is right!!
%         fo(iSubj).Number_of_spec = fo(iSubj).Number_of_spec - numel(TRs_to_toss{ROI_idx}{find_toss_idx});
%         n_toss = n_toss +1;
%         warning(['Tossing TRs ' num2str(TRs_to_toss{ROI_idx}{find_toss_idx}) ...
%             ' for ' fo(iSubj).sujet_name ', as requested...']);
%     end
% end
% if n_toss > 0
%     warndlg('Crudely tossing some TRs, as specified!')
% end

% previously excluded mps 20200909

%% freq & phase correction 
corr_str = {'corr','correl','correlation','correlate'};
corr_bound_str = {'corr with bound','correlation with bound','corr with bounds',...
    'correlation with bounds','corr. with bound','corr. with bounds',...
    'corr w/ bound','corr w/ bounds','correlation w/ bound','correlation w/ bounds',...
    'corr_bound','correlation_bound'};
max_abs_real_str = {'(max + abs) x2 + real', 'max + abs, real','max abs real'};

if sum(strcmp(lower(options.freq_phase_corr) , corr_str))
%%% option #1 = Corr
    par.method='correlation';
    par.correlation_bound = [];
    fc = processing_spec(fo,par);

    % par.same_fig = 1; % turn this on to plot on the same figure;
    % par.arg_give_me_a_new_one = 1; % choose figure number
    for iSubj = 1:numel(fo)
    plot_spectrum(fo(iSubj),par) % before processing
    get_ax = axis;
    axis([0.5 2.5 get_ax(3) get_ax(4)]);
    title([fo(iSubj).sujet_name ' Before'])
    set(gcf,'color','w')
    end

    for iSubj = 1:numel(fo)
    plot_spectrum(fc(iSubj),par) % after processing
    get_ax = axis;
    axis([0.5 2.5 get_ax(3) get_ax(4)]);
    title([fc(iSubj).sujet_name ' 1) Corr.'])
    set(gcf,'color','w')
    end

% % to know if frequency and phase correction were successful, look for
% % variability in the peak freq of NAA (~2ppm - for freq) and variability in
% % the shape (e.g., weird dips) on the flanks of the NAA peak (~1.95ppm -
% % for phase).

elseif sum(strcmp(lower(options.freq_phase_corr) , corr_bound_str))
%%% option #2 = corr w/ bounds
    if ~isfield(par,'correlation_bound')
        par.correlation_bound = [1.8 3.5];
        warning(['Using default correlation bounds for freq. & phase correction: '...
            num2str(par.correlation_bound) ]);
    end
    par.method='correlation';
    fc = processing_spec(fo,par);

    for iSubj = 1:numel(fo)
        plot_spectrum(fc(iSubj),par) % after processing
        get_ax = axis;
        axis([0 6 get_ax(3) get_ax(4)]);
        title([fo(iSubj).sujet_name ' 2) Corr. w/ bound'])
        set(gcf,'color','w')
    end

elseif sum(strcmp(lower(options.freq_phase_corr) , max_abs_real_str))
%%% option #3 = (max + abs) x2 & real

    par.method='max';
    par.correct_freq_mod='abs';
    fc_1 = processing_spec(fo,par);
    if options.close_all
        close all
    end
    par.correct_freq_mod='abs';
    fc_2 = processing_spec(fc_1,par);
    if options.close_all
        close all
    end
    par.correct_freq_mod='real';
    fc = processing_spec(fc_2,par);
    if options.close_all
        close all
    end
    
else
    error(['Unknown value for options.freq_phase_corr: ' options.freq_phase_corr]);
end

%% show & save some figures for freq. & phase correction
% plot_dir = fullfile('/home/shaw-raid1/data/MRS/saved_plots',options.which_ROI);
% plot_dir_contents = dir([plot_dir '/*.fig']);
% if ~isempty(plot_dir_contents)
%     clear ask_wipe
%     ask_wipe = ['Saved plots folder is not empty (' plot_dir ...
%         '), do you want to delete all contents and re-write, or quit?'];
%     wipe_dir = questdlg(ask_wipe,'Do you want to delete?','Delete','Quit','Quit');
%     if strcmp(wipe_dir,'Delete')
%         warning(['Deleting contents of ' plot_dir ', as requested...']);
%         eval(['! rm -f ' plot_dir '/*.fig']);
%     elseif strcmp(wipe_dir,'Quit')
%         return
%     end
% end

% if options.displayFigs & par.figure % if we are plotting figures...
    for iSubj = 1:numel(fo)
%         plot_spectrum(fo(iSubj),par) % before processing
%         get_ax = axis;
%         axis([0 6 get_ax(3) get_ax(4)]);
%         title([fo(iSubj).sujet_name ' uncorrected'])
%         set(gcf,'color','w')
%         savefig(fullfile(plot_dir,[fo(iSubj).sujet_name '_' options.which_ROI '_uncorrected.fig']))
%         if options.close_all
%             close(gcf);
%         end
%         
%         plot_spectrum(fc(iSubj),par) % after processing
%         get_ax = axis;
%         axis([0 6 get_ax(3) get_ax(4)]);
%         title([fo(iSubj).sujet_name ' 3) Max + Abs. x2 & Real'])
%         set(gcf,'color','w')
%         savefig(fullfile(plot_dir,[fo(iSubj).sujet_name '_' options.which_ROI '_corrected.fig']))
%         if options.close_all
%             close(gcf);
%         end
%         
        
        output.SD_NAA(iSubj) = mps_SD_NAA(fc(iSubj),par);
    end
% end

%% to obtain water linewidth
fw_ecc = get_water_width(fw_ecc, par.figure);
output.water_width_no_cor = [];
for iS = 1:numel(fw_ecc)
    output.water_width_no_cor(iS,1) = fw_ecc(iS).water_width_no_cor;
end
% gosia does not trust the correction, so don't use it
% (already have done eddy current compensation) - will be a bit smaller
% than the linewidth of water we measured at the scanner

%% preparing for LCModel analysis:

pp.root_dir=['/home/shaw-raid1/data/MRS/noise_quant/processed_data'];
% it uses fc.sujet_name and fc.SerDescr for naming RAW files
pp.gessfile=4; % different options for setting file names based on dicom info
% mps 20200909 unlike regular PHCP =5, use =4 here to use ser_dir

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

write_fid_to_lcRAW(fc,pp,fw_ecc);

% remove existing lcm files
lcm_fit_dir = fullfile(pp.root_dir,'lcm_fit');
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

output.fo = fo;
output.fc = fc;

%% run the LCModel analysis
elseif strcmp(lower(options.which_step), 'lcm')

[foo, server_name] = system('hostname');
if ( strcmp(server_name(1:end-1),'kea.cmrr.umn.edu') || ...
        strcmp(server_name(1:end-1),'lcmodel.cmrr.umn.edu') ) && ...
        datenum(version('-date')) <= datenum('February 12, 2009')

    processing_LCmodel('steam_phcp_lipids','lcm_fit')
    % n.b. matlab will only run on kea if you select a 32 bit installation,
    % such as 2009a...
    
else
    error('This command needs to be run on the server lcmodel.cmrr.umn.edu -- ssh over there and try again. You probably should use a 32-bit version of Matlab, such as 2009a...');
end

%% post-LCModel processing
elseif strcmp(lower(options.which_step), 'post-lcm')

%% check matlab version & server
if datenum(version('-date')) > datenum('February 11, 2014')
    error(['This code relies on SPM8, and will NOT function with versions '...
        'of Matlab beyond 2014a!!'])
end

[foo, server_name] = system('hostname');
if strcmp(server_name(1:19),'range4.cmrr.umn.edu')
    
    dir_res=get_subdir_regex;
    concat_ps_file(dir_res)
else
    error('This command needs to be run on the server range4.cmrr.umn.edu -- ssh over there and try again...');
end

%% convert the data files to .pdf and .csv
% 1. generate pdf file for review of the fits
% be sure to select the folder you wrote the LCmodel data to in the step
% above, e.g., lcm_fit

if ~exist('proc_dir','var')
    proc_dir = '/home/shaw-raid1/data/MRS/noise_quant/processed_data';
end

% dir_res=get_subdir_regex;
% concat_ps_file(dir_res)

% 2. generate spreadsheet

% fix file permissions for the files you just wrote here...
eval(['! chmod g+w -R ' dir_res{1}]);

% check if we know which ROI to use
if numel(dir_res) > 1
    error('You selected more than 1 lcm folder...')
end


c=get_result(dir_res);
write_name = fullfile(proc_dir,...
    [datestr(now,'yyyymmdd') '_phcp_' options.which_ROI '_' num2str(numel(c.suj)) 'subj_unscaled.csv']);
if ~exist(write_name,'file')
    write_conc_res_to_csv(c,write_name)
else
    error(['The following file already exists, to create a new version, you must delete the existing one! - ' write_name]);
end

% scale to tCr
c_cr = correct_result_ration(c,'','Cr_PCr',8);
write_name = fullfile(proc_dir,...
    [ datestr(now,'yyyymmdd') '_phcp_' options.which_ROI '_' num2str(numel(c.suj)) 'subj_Cr_scaled.csv']);
if ~exist(write_name,'file')
    write_conc_res_to_csv(c_cr,write_name)
else
    error(['The following file already exists, to create a new version, you must delete the existing one! - ' write_name]);
end

% scale to Water
% first get all tissue fractions
tissue_fract_data = mps_run_all_MRS_tissue_correct;

cw = c;
cw.fgray = [];
cw.fwhite = [];
cw.fcsf = [];

for iSubj = 1:numel(cw.suj)
    subj_date = [str2num(cw.suj{iSubj}(2:8)) datenum(cw.suj{iSubj}(10:17),'yyyymmdd')];
    subj_idx = find(tissue_fract_data.subj_date(:,1) == subj_date(1) & ...
        tissue_fract_data.subj_date(:,2) == subj_date(2));
    cw.fgray(iSubj) = tissue_fract_data.tissue_fract(subj_idx, ROI_idx, 1); % GM first
    cw.fwhite(iSubj) = tissue_fract_data.tissue_fract(subj_idx, ROI_idx, 2); % WM second
    cw.fcsf(iSubj) = tissue_fract_data.tissue_fract(subj_idx, ROI_idx, 3); % CSF third
end

corr_factor = water_content_attenuation(cw, par);
c_wscale = correct_result(cw, corr_factor);

write_name = fullfile(proc_dir,...
    [ datestr(now,'yyyymmdd') '_phcp_' options.which_ROI '_' num2str(numel(c.suj)) 'subj_H2O_scaled.csv']);
if ~exist(write_name,'file')
    write_conc_res_to_csv(c_wscale,write_name)
else
    error(['The following file already exists, to create a new version, you must delete the existing one! - ' write_name]);
end

% rename .pdf and .ps

write_name = fullfile(proc_dir,...
    [datestr(now,'yyyymmdd') '_phcp_' options.which_ROI '_' num2str(numel(c.suj)) 'subj.pdf']);
run_cmd = ['! mv ' proc_dir '/' options.which_ROI '_lcm_fit.pdf ' write_name ' ;'];
eval(run_cmd);

write_name = fullfile(proc_dir,...
    [datestr(now,'yyyymmdd') '_phcp_' options.which_ROI '_' num2str(numel(c.suj)) 'subj.ps']);
run_cmd = ['! mv ' proc_dir '/' options.which_ROI '_lcm_fit.ps ' write_name ' ;'];
eval(run_cmd);

end % options.which_step

%% output
options.par = par;
output.options = options;
output.date_run = datestr(now);

end