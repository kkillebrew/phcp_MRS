function [output] = look_at_mrs_group_differences(options)
% usage: [output] = look_at_mrs_group_differences(options)
%
% mps 20181123
%% opt
addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode')); % add mpsCode to path

if ~exist('options','var')      
    options = [];
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 1; % 1 = yes, 0 = no
end
if ~isfield(options,'mrs_file')
    error('No options.mrs_file provided!')
    % e.g., options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20210422_phcp_OCC_149subj_H2O_scaled.csv';
end
if ~isfield(options,'mrs_n_col')
    %     KWK - using the actual number of columns in csv file (68) (after editing the csv file, by deleting the extra column
%     labels ('met' and 'cormat') and getting rid of the random commas that
%     presumably created).
%     options.mrs_n_col = 68; % These numbers (507 and 439) don't seem to work w/ the output that
    options.mrs_n_col = 504; % 507 = if using Gosia's notes, 439 = LCM default
end
if ~isfield(options,'mrs_header_lines')
    options.mrs_header_lines = 6; % 6 = if using Gosia's notes, 8 = LCM default
%     options.mrs_header_lines = 4; % Worked for KWK along w/ above noted changes.
end
if ~isfield(options,'mrs_overwrite_mat')
    options.mrs_overwrite_mat = 0; % by default, load exisiting .mat file, 1 to overwrite
end
mrs_opt.target_file = options.mrs_file;
mrs_opt.n_col = options.mrs_n_col;
mrs_opt.header_lines = options.mrs_header_lines;
mrs_opt.overwrite_mat = options.mrs_overwrite_mat;
if ~isfield(options,'mrs_struct')
    options.mrs_struct = read_in_LCModel_results(mrs_opt);
end

if ~isfield(options,'toss_subj_num_date')
    options.toss_subj_num_date = []; % this subjects looks like they have bad MRS data on 20190331
    %options.toss_subj_num_date = [6004213 datenum('20190116','yyyymmdd')]; % this subjects looks like they have bad MRS data on 20190331
end
if ~isfield(options,'which_metab')
%      options.which_metab = {'Glu','GABA'};
    options.which_metab = {'MacY','Asc','Asp','PCho','GPC','Cr','PCr','GABA','Glc','Gln','Glu','GSH','Ins','Lac','NAA','NAAG','PE','sIns','Tau'};   % all metabs
    warning('options.which_metab not specified, assuming you want to look at only Glu and GABA...')
end
if ~isfield(options,'avg_repeats')
    options.avg_repeats = 1;
end
if ~isfield(options,'toss_CRLB')
    options.toss_CRLB = 0; % set value to 0 to turn off, 1 for on
end
if ~isfield(options,'toss_mrs_quality')
    options.toss_mrs_quality = 1;
end
if ~isfield(options,'quality_file')
    options.quality_file = []; % e.g., options.quality_file = '/home/shaw-raid1/data/MRS/processed_data/data_quality/data_coil_parameters_performance_OCC_20200420.csv';
end
if options.toss_mrs_quality
    if ~isempty(options.quality_file)
        qual_opt.target_file = options.quality_file;
        load_quality = read_in_gosia_data_quality( qual_opt );
        quality_data = load_quality.csv_data;
    else
        error('you set options.toss_mrs_quality = 1, but options.quality_file is empty or missing, so I quit!')
    end    
end
if ~isfield(options,'outlier_n_SDs')
    options.outlier_n_SDs = 3;
end
if ~isfield(options,'corr_metabs')
    options.corr_metabs = 0;
end
if ~isfield (options, 'subj_group_def')
    options.subj_group_def = 1; % 1 = controls, relatives, probands; 
    % 2 = controls, SZ, BP
    % 3 = SZ, schizoaffective (SCA), BP; 
    % 4 = healthy (con+rel), SZ+SCA, bipolar,
    % 5 = controls, probands, relatives (flip order of P & R)
end
if ~isfield(options,'corr_symp')
    options.corr_symp = 0; % correlate symptoms
end
if ~isfield(options,'which_symp')
    options.which_symp = 'SGITotal';
    % valid options currently are SGITotal, BPRSDisorganization, BPRSTotal, SAPSGlobalPositiveSymptoms
end
if ~isfield(options,'plot_groups_separate')
    options.plot_groups_separate = 1;
end
if ~isfield(options,'show_stars')
    options.show_stars = 0;
end
if ~isfield(options,'combine_Cr_Cho')
    options.combine_Cr_Cho = 1;
    
    cr_idx = strcmp(options.which_metab , 'Cr');
    pcr_idx = strcmp(options.which_metab , 'PCr');
    gpc_idx = strcmp(options.which_metab , 'GPC');
    pcho_idx = strcmp(options.which_metab , 'PCho');
    
    if sum(cr_idx) && sum(pcr_idx) 
        warning('Combining Cr and PCr as tCr for decoding...');
        options.which_metab(logical(cr_idx + pcr_idx)) = [];
        options.which_metab{end+1} = 'tCr';
        options.mrs_struct.tCr = options.mrs_struct.Cr + options.mrs_struct.PCr;
    end    
    if sum(gpc_idx) && sum(pcho_idx)
        warning('Combining GPC and PCho as tCho for decoding...');
        options.which_metab(logical(gpc_idx + pcho_idx)) = [];
        options.which_metab{end+1} = 'tCho';
        options.mrs_struct.tCho = options.mrs_struct.GPC + options.mrs_struct.PCho;
    end
end

%% get mrs file info
for iFile = 1:numel(options.mrs_struct.row_name)-2 % skip last 2, mean and sd/mean
    name_idx = regexp(options.mrs_struct.row_name{iFile},'P\d\d\d\d\d\d\d');
    options.mrs_struct.subj_number(iFile,1) = str2num(options.mrs_struct.row_name{iFile}...
        (name_idx+1:name_idx+7));
    date_idx = regexp(options.mrs_struct.row_name{iFile},'\d\d\d\d\d\d\d\d');
    options.mrs_struct.date_number(iFile,1) = datenum(options.mrs_struct.row_name{iFile}...
        (date_idx:date_idx+7),'yyyymmdd');
end

%% toss subj
MRS_subj_date = [options.mrs_struct.subj_number options.mrs_struct.date_number];

toss_idx = [];
if ~isempty(options.toss_subj_num_date)
    for iSubj = 1:size(MRS_subj_date,1)
        for iToss = 1:size(options.toss_subj_num_date,1)
            if sum(MRS_subj_date(iSubj,:) == options.toss_subj_num_date(iToss,:)) == 2
                % both subject num and date are the same, so toss
                toss_idx = [toss_idx ; iSubj];
            end
        end
    end
    if ~isempty(toss_idx)
        warning(['Tossing ' num2str(numel(toss_idx)) ' data sets, as requested...']);
        MRS_subj_date(toss_idx,:) = [];
    end
end

find_out_mrs = [];

if options.toss_mrs_quality

    mrs_qual_metrics = {'SNR','lw_H2O','lw_tCr'};
    high_is_good = [1 0 0];
    
    for iQ = 1:numel(mrs_qual_metrics)
        if ~sum(ismember(quality_data.Properties.VariableNames, ...
                mrs_qual_metrics{iQ}))
            error(['You asked to toss data based on MRS quality metrics, but '...
                'the following metric is not included in your MRS quality data set: '...
                mrs_qual_metrics{iQ} ' -- Are you loading data quality '...
                'metrics from Gosia??']);
        end
        
        use_qual = [];
        for iS = 1:size(MRS_subj_date,1)
            subj_idx = find(quality_data.subject_num == MRS_subj_date(iS,1) & ...
                quality_data.date_num == MRS_subj_date(iS,2));
            if ~isempty(subj_idx)
                use_qual(iS,1) = quality_data.(mrs_qual_metrics{iQ})(subj_idx);
            else
                use_qual(iS,1) = NaN;
            end
        end
                
        if high_is_good(iQ)
            toss_me = find( (use_qual - nanmean(use_qual,1) )...
                < -options.outlier_n_SDs * nanstd(use_qual,0,1) );
            
        else
            toss_me = find( (use_qual - nanmean(use_qual,1) )...
                > options.outlier_n_SDs * nanstd(use_qual,0,1) );
        end
        output.tossed_subj_date{iQ} = MRS_subj_date(toss_me,:);
        find_out_mrs = [find_out_mrs ; toss_me];
    end
    
    if ~isempty(find_out_mrs)
        MRS_subj_date(find_out_mrs,:) = [];
        warning(['tossing ' num2str(numel(unique(find_out_mrs))) ...
            ' subject(s) for poor MRS data quality.']);
    end
end
toss_idx = [toss_idx ; find_out_mrs];

%% if we're using symptom levels, find them BEFORE tossing repeated data
if options.corr_symp
    if strcmp(options.which_symp,'SGITotal') || strcmp(...
            options.which_symp,'SAPSGlobalPositiveSymptoms')
        date_cutoff = inf;
    elseif strcmp(options.which_symp,'BPRSTotal') || strcmp(...
            options.which_symp,'BPRSDisorganization') 
        date_cutoff = 30;
    end
    
    addpath(genpath('/home/shaw-raid1/matlab_tools/COP_scripts'));
    psych_opts.target_file = '/home/shaw-raid1/data/7T/demographics/PHCP7TfMRIPsych.csv';
    psych_data = read_in_psych_data(psych_opts);
    
    all_symp = psych_data.(options.which_symp);
    
    corr_symp = nan(size(MRS_subj_date(:,1)));
    
    for iS = 1:numel(corr_symp)
        symp_idx = find( strcmp(['P' num2str(MRS_subj_date(iS,1))] , ...
            psych_data.Record_ID ));
        if isempty(symp_idx)
            error('cannot find subject number in psych table!');
        end
        date_diff = MRS_subj_date(iS,2) - cell2mat(...
            psych_data.datenumber(symp_idx) );
        [~, sort_idx] = sort(abs(date_diff),'ascend'); % find the most recent symptom dates
        
        for iSort = 1:numel(sort_idx)
            if (date_diff(sort_idx(iSort)) <= date_cutoff) & (~isnan(all_symp...
                    (symp_idx(sort_idx(iSort)))) ) % use most recent date that is below cuttoff and not missing
                corr_symp(iS) = all_symp(symp_idx(sort_idx(iSort)));
                break
            end
        end
    end
end
%% average repeated data
if options.avg_repeats
    [~, unique_idx] = unique(MRS_subj_date(:,1),'rows','stable');
    [~, repeat_scan_idx] = setxor(MRS_subj_date(:,:),...
        MRS_subj_date(unique_idx,:),'rows','stable'); % find not intersection
    
    output.retest.date.mean = mean(MRS_subj_date(repeat_scan_idx,2) - ...
        MRS_subj_date(repeat_scan_idx-1,2));
    output.retest.date.median = median(MRS_subj_date(repeat_scan_idx,2) - ...
        MRS_subj_date(repeat_scan_idx-1,2));
    output.retest.date.sd = std(MRS_subj_date(repeat_scan_idx,2) - ...
        MRS_subj_date(repeat_scan_idx-1,2));
    output.retest.date.range = [min(MRS_subj_date(repeat_scan_idx,2) - ...
        MRS_subj_date(repeat_scan_idx-1,2)) ...
        max(MRS_subj_date(repeat_scan_idx,2) - ...
        MRS_subj_date(repeat_scan_idx-1,2))];
    
    output.retest.group_Ns.control = sum(MRS_subj_date(repeat_scan_idx,1) ...
        < 2000000);
    output.retest.group_Ns.relative = sum(MRS_subj_date(repeat_scan_idx,1) ...
        >= 2000000 & MRS_subj_date(repeat_scan_idx,1) < 6000000);
    output.retest.group_Ns.proband = sum(MRS_subj_date(repeat_scan_idx,1) ...
        >= 6000000);
    
    for iRep = 1:numel(repeat_scan_idx)
        if MRS_subj_date(repeat_scan_idx(iRep),1) ~= MRS_subj_date(repeat_scan_idx(iRep)-1,1)
            % make sure that the repeated scans are immediately after the
            % 1st scan for that subject in the list
           error(['Session # ' num2str(repeat_scan_idx(iRep)) ' is a repeat,'...
               'but is not the same subj# as the scan before it??']);
        end
        if repeat_scan_idx(iRep) > 2
            if MRS_subj_date(repeat_scan_idx(iRep),1) == MRS_subj_date(repeat_scan_idx(iRep)-2,1)
                % make sure there are no 3x repeats..
                error(['Session # ' num2str(repeat_scan_idx(iRep)) ' is a 3x repeat??']);
            end
        end
    end
    % average repeats with the scan before in the list (should be the
    % 1st scan)
    MRS_subj_date_with_repeats = MRS_subj_date;
    warning('averaging repeated subjects, as requested...')
    MRS_subj_date(repeat_scan_idx,:) = [];
end
%% Subject group indices (from refit_COP_dataRK) added by MPS 13 NOV 2019
demog_opts.target_file = '/home/shaw-raid1/data/7T/demographics/PHCP7TfMRIDemo.csv';
demog_data = read_in_demog_data(demog_opts);

subj_number = MRS_subj_date(:,1);

dx_list = nan(numel(subj_number),1);
missing_dx_list = [];
for iSubj = 1:numel(subj_number)
    dx_idx = strcmp(['P' num2str(subj_number(iSubj))],demog_data.Record_ID);
    if isempty(dx_idx) % if this subject isn't in the demographics .csv file
        missing_dx_list = [missing_dx_list ; subj_number(iSubj)];
        continue
    end
    dx_list(iSubj) = demog_data.Dx_code(dx_idx);
end

%Legend for dx codes: 0=none; 1=MDD; 2=SZ; 3=SZaff; 4=BP1; 5=BP2;
%6=Panic; 7=DeprNOS; 8=PsychNOS; 9=ADHD

if options.subj_group_def == 1; %use controls, probands, and relatives as the different groups
    
    use_colors = {'g','b','r'};
    g1_idx = subj_number < 2000000;
    g1_label = 'Control';
    g1_short = 'Ctrl';
    g2_idx = subj_number >= 2000000 & subj_number < 6000000;
    g2_label = 'Relative';
    g2_short = 'Rel';
    g3_idx = subj_number > 6000000;
    g3_label = 'Proband';
    g3_short = 'Psy';
    
elseif options.subj_group_def == 2; % look at controls, SZ, BP
    use_colors = {'g',[0.75 0 0.75],[255 204 0]./255};
    g1_idx = (dx_list == 0) & subj_number < 2000000; % no dx, and control subject ID
    g1_label = 'Control';
    g1_short = 'C';
    g2_idx = (dx_list == 2) & subj_number > 6000000;
    g2_label = 'SZ';
    g2_short = 'SZ';
    g3_idx = ( (dx_list == 4 | dx_list == 5) ) & subj_number > 6000000;
    g3_label = 'BP';
    g3_short = 'BP';
    
elseif options.subj_group_def == 3; %look at SZ, SCA, BP
    use_colors = {[0.75 0 0.75],[255 102 0]./255 , [255 204 0]./255};
    g1_idx = (dx_list == 2) & subj_number > 6000000;
    g1_label = 'SZ';
    g1_short = 'SZ';
    g2_idx = (dx_list == 3) & subj_number > 6000000;
    g2_label = 'SCA';
    g2_short = 'SCA';
    g3_idx = ( (dx_list == 4 | dx_list == 5) ) & subj_number > 6000000;
    g3_label = 'BP';
    g3_short = 'BP';
    
elseif options.subj_group_def == 4; %look at controls, SCZ+SCA, BP
    use_colors = {'c','m',[255 204 0]./255};
    g1_idx = (dx_list == 0) & subj_number < 6000000;
    g1_label = 'Ctrl + Rel';
    g1_short = 'C+R';
    g2_idx = (dx_list == 2 | dx_list == 3) & subj_number > 6000000;
    g2_label = 'SZ+SCA';
    g2_short = 'SZ+A';
    g3_idx = ( (dx_list == 4 | dx_list == 5) ) & subj_number > 6000000;
    g3_label = 'BP';
    g3_short = 'BP';
    
elseif options.subj_group_def == 5; % flip the order for probands and relatives...
    
    use_colors = {'g','r','b'};
    g1_idx = subj_number < 2000000;
    g1_label = 'Control';
    g1_short = 'C';
    g2_idx = subj_number > 6000000;
    g2_label = 'Proband';
    g2_short = 'P';
    g3_idx = subj_number >= 2000000 & subj_number < 6000000;
    g3_label = 'Relative';
    g3_short = 'R';

else
    error(['Unknown value for options.subj_group_def = ' num2str(options.subj_group_def)]);
end
g1_idx_bin = g1_idx;
g2_idx_bin = g2_idx;
g3_idx_bin = g3_idx;
g1_idx = find(g1_idx); % keep binary indices, but also make numeric
g2_idx = find(g2_idx);
g3_idx = find(g3_idx);

%% look at group differences in whichever metabolites were selected...
all_p_to_FDR = [];

for iM = 1:numel(options.which_metab)
    metab = options.mrs_struct.(options.which_metab{iM})(1:end-2); % toss last 2 values - mean and sd/mean    
    metab(toss_idx) = []; % toss excluded subj
    
    checkSD = 0;
    eval(['checkSD = isfield(options.mrs_struct,''SD' options.which_metab{iM} ''');']); % toss last 2 values - mean and sd/mean

    if checkSD
        eval(['metab_CRLB = options.mrs_struct.SD' options.which_metab{iM} '(1:end-2);']); % toss last 2 values - mean and sd/mean
        metab_CRLB(toss_idx) = []; % toss excluded subj


        if options.toss_CRLB % remove values that are unreliable (CRLB too high)
            warning(['tossing ' num2str(sum(metab_CRLB > options.toss_CRLB)) ...
                ' ' options.which_metab{iM} ' values for CRLB > ' num2str(...
                options.toss_CRLB)]);
            metab( metab_CRLB > options.toss_CRLB ) = NaN;
        end
        this_metab.CRLB = metab_CRLB;
    end
    
    clear this_metab
    if options.avg_repeats
        
        this_metab.stats.corr_retest.type = 'Spearman';
        
        use_metab_corr_idx = ~isnan(metab(repeat_scan_idx-1)) & ~isnan(...
            metab(repeat_scan_idx));
        
        repeat_scan_idx = repeat_scan_idx( use_metab_corr_idx ); % exclude NaNs
        [this_metab.stats.corr_retest.r, this_metab.stats.corr_retest.p] = ...
            corr(metab(repeat_scan_idx-1),metab(repeat_scan_idx),...
            'type',this_metab.stats.corr_retest.type);
        this_metab.stats.corr_retest.df = numel(repeat_scan_idx)-2;
        
        this_metab.stats.ICC_retest = ICC(1, 'k', [metab(repeat_scan_idx-1) ...
            metab(repeat_scan_idx)]);
        
        if options.displayFigs
            figure; hold on
            
            [poly_fit] = polyfit(metab(repeat_scan_idx-1), ...
                metab(repeat_scan_idx), 1);
                        
            fit_x = [min(metab(repeat_scan_idx-1)) max(metab(repeat_scan_idx-1))];
            fit_y = poly_fit(1).*fit_x + poly_fit(2);
            y_range = [min(metab(repeat_scan_idx)) max(metab(repeat_scan_idx))];
            
            g1_with_rep = find( MRS_subj_date_with_repeats(:,1) < 2000000 );
            g2_with_rep = find( MRS_subj_date_with_repeats(:,1) >= 2000000 & ...
                MRS_subj_date_with_repeats(:,1) < 2000000 );
            g3_with_rep = find( MRS_subj_date_with_repeats(:,1) >= 6000000 );
            rep_g1_idx = intersect(repeat_scan_idx, g1_with_rep);
            rep_g2_idx = intersect(repeat_scan_idx, g2_with_rep);
            rep_g3_idx = intersect(repeat_scan_idx, g3_with_rep);

            if options.plot_groups_separate
                plot_colors = {[0.33 1 0.33],[0.33 0.33 1],[1 0.33 0.33]};
            else
                plot_colors = {'w','w','w'};
            end
            
            use_legend = [];
            idx_leg = 0;
            if ~isempty(rep_g1_idx)
                plot(-1, -1 ,'ko',...
                    'MarkerFaceColor',plot_colors{1},'linewidth',2,'MarkerSize',8)
                idx_leg = idx_leg+1;
                use_legend{idx_leg} = [g1_short ', n = ' num2str(numel(rep_g1_idx))];
            end
            if ~isempty(rep_g2_idx)
                plot(-1, -1 ,'ko',...
                    'MarkerFaceColor',plot_colors{2},'linewidth',2,'MarkerSize',8)
                idx_leg = idx_leg+1;
                use_legend{idx_leg} = [g2_short ', n = ' num2str(numel(rep_g2_idx))];
            end
            if ~isempty(rep_g3_idx)
                plot(-1, -1 ,'ko',...
                    'MarkerFaceColor',plot_colors{3},'linewidth',2,'MarkerSize',8)
                idx_leg = idx_leg+1;
                use_legend{idx_leg} = [g3_short ', n = ' num2str(numel(rep_g3_idx))];
            end
            plot(fit_x,fit_y,'k-','linewidth',2)

            plot(metab(rep_g1_idx-1), metab(rep_g1_idx), 'ko', ...
                'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', plot_colors{1})
            plot(metab(rep_g2_idx-1), metab(rep_g2_idx), 'ko', ...
                'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', plot_colors{2})
            plot(metab(rep_g3_idx-1), metab(rep_g3_idx), 'ko', ...
                'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', plot_colors{3})

            ylabel([options.which_metab{iM} ' scan 2 (mM)'],'Color','k');
            xlabel([options.which_metab{iM} ' scan 1 (mM)'],'Color','k');
            set(gcf,'color','w','POS',[357   275   560   560])
            box off
            set(gca,'XColor','k','YColor','k','Fontsize',18)

            both_range = [ min([metab(repeat_scan_idx-1) ; metab(repeat_scan_idx)]) ...
                max([metab(repeat_scan_idx-1) ; metab(repeat_scan_idx)]) ];
            both_diff = both_range(2) - both_range(1);
            axis([both_range(1)-both_diff*0.1 both_range(2)+both_diff*0.1 ...
                both_range(1)-both_diff*0.1 both_range(2)+both_diff*0.1])
            text(fit_x(1),both_range(2),...
                ['n = ' num2str(numel(repeat_scan_idx))],'fontsize',18)
            text(fit_x(1),both_range(2) - both_diff*0.1,...
                ['ICC = ' num2str(round(100*this_metab.stats.ICC_retest)/100)],'fontsize',18)
            use_tick = get(gca,'YTick');
            set(gca,'XTick',use_tick,'DataAspectRatio',[1 1 1]);
            
            if options.plot_groups_separate
                legend(use_legend)
            end
            
            pause(1)
        end
        
        avg_metab = -1*ones(size(metab));        
        avg_metab(unique_idx) = metab(unique_idx);
        avg_metab(repeat_scan_idx-1) = nanmean([metab(repeat_scan_idx-1) ...
            metab(repeat_scan_idx)],2); % avg repeated subj
        metab = avg_metab(~(avg_metab < 0)); % do it this way, to deal with NaNs for tossed values with CRLB > limit
    end
    
    this_metab.data = metab;
    this_metab.g1_data = metab(g1_idx);
    this_metab.g1_mean = mean(this_metab.g1_data);
    this_metab.g1_median = median(this_metab.g1_data);
    this_metab.g1_sd = std(this_metab.g1_data);
    
    this_metab.g2_data = metab(g2_idx);
    this_metab.g2_mean = mean(this_metab.g2_data);
    this_metab.g2_median = median(this_metab.g2_data);
    this_metab.g2_sd = std(this_metab.g2_data);
    
    this_metab.g3_data = metab(g3_idx);
    this_metab.g3_mean = mean(this_metab.g3_data);
    this_metab.g3_median = median(this_metab.g3_data);
    this_metab.g3_sd = std(this_metab.g3_data);
    
    [h p ci this_metab.stats.ttest2_g1_g2] = ttest2(this_metab.g1_data, ...
        this_metab.g2_data);
    this_metab.stats.ttest2_g1_g2.p = p;
    
    showKW = 'off';
    
    group_data = [this_metab.g1_data ; this_metab.g2_data ; ...
        this_metab.g3_data];
    group_idx = [ones(numel(this_metab.g1_data),1) ; ...
        2*ones(numel(this_metab.g2_data),1) ; ...
        3*ones(numel(this_metab.g3_data),1)];
    
    [p, this_metab.stats.kruskallwallis.table, ...
        this_metab.stats.kruskallwallis.stats] = kruskalwallis(...
        group_data, group_idx, showKW);
    this_metab.stats.kruskallwallis.p = p;
    all_p_to_FDR = [all_p_to_FDR p]; % placeholder, come back at the end and correct
    this_metab.stats.kruskallwallis.p_bonf = p*numel(options.which_metab);
    if this_metab.stats.kruskallwallis.p_bonf > 1
        this_metab.stats.kruskallwallis.p_bonf = 1;
    end
    
    [p, this_metab.stats.levene] = vartestn(group_data, group_idx, ...
        'display', showKW, 'testtype', 'LeveneAbsolute');
    this_metab.stats.levene.p = p;
          
    % correlate metab with symptoms, if requested
    if options.corr_symp
        
        if iM == 1 % first time through, need to toss / avg sympt data
            corr_symp(toss_idx) = NaN;
            
            if options.avg_repeats
                avg_symp = -1*ones(size(corr_symp));
                avg_symp(unique_idx) = corr_symp(unique_idx);
                avg_symp(repeat_scan_idx-1) = nanmean([corr_symp(repeat_scan_idx-1) ...
                    corr_symp(repeat_scan_idx)],2); % avg repeated subj
                corr_symp = avg_symp(~(avg_metab < 0)); % do it this way, to deal with NaNs for tossed values with CRLB > limit
            end
        end
        
        this_metab.corr_symp.name = options.which_symp;
        this_metab.corr_symp.type = 'Spearman';
        
        corr_symp_idx = find(~isnan(corr_symp) & ~isnan(metab));
        [this_metab.corr_symp.r, this_metab.corr_symp.p] = corr(...
            metab(corr_symp_idx), corr_symp(corr_symp_idx), 'type',...
            this_metab.corr_symp.type);
        this_metab.corr_symp.df = numel(corr_symp_idx)-2;
        
        if options.displayFigs
            figure; hold on
            pause(1) % let this finish plotting...
            
            [poly_fit] = polyfit(metab(corr_symp_idx), ...
                corr_symp(corr_symp_idx), 1);
                        
            fit_x = [min(metab(corr_symp_idx)) max(metab(corr_symp_idx))];
            fit_y = poly_fit(1).*fit_x + poly_fit(2);
            y_range = [min(corr_symp(corr_symp_idx)) max(corr_symp(corr_symp_idx))];
            plot(fit_x,fit_y,'k-','linewidth',2)

            plot(metab(corr_symp_idx), corr_symp(corr_symp_idx), 'ko', ...
                'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w')
            
            text(fit_x(1),y_range(2) - (y_range(2) - y_range(1))*0,...
                ['n = ' num2str(numel(corr_symp_idx))],'fontsize',18)
            text(fit_x(1),y_range(2) - (y_range(2) - y_range(1))*0.1,...
                ['r = ' num2str(round(100*this_metab.corr_symp.r)/100)],'fontsize',18)
            text(fit_x(1),y_range(2) - (y_range(2) - y_range(1))*0.2,...
                ['p = ' num2str(round(100*this_metab.corr_symp.p)/100)],'fontsize',18)

            ylabel([options.which_symp ' score'],'Color','k');
            if strcmp(options.which_metab{iM},'MacY')
                use_x_label = 'MacM (arb. units)';
            else
                use_x_label = [options.which_metab{iM} ' (mM)'];
            end
            xlabel(use_x_label,'Color','k');
            set(gcf,'color','w','POS',[357   275   560   560])
            box off
            set(gca,'XColor','k','YColor','k','Fontsize',18)
        end
    end
    
    output.(options.which_metab{iM}) = this_metab;
        
%% plot, if requested 
    if options.displayFigs        
        figure; hold on
        pause(1) % let this finish plotting...
        
        plot([0 0],repmat(this_metab.g1_median,[1 2]),'-','color',[0.5 0.5 0.5],...
            'linewidth',4)
        plot([0 0],repmat(this_metab.g1_median,[1 2]),'-','color',[0.5 0.5 0.5],...
            'linewidth',2)
        plot([0 0],repmat(this_metab.g1_median,[1 2]),'--','color',[0.5 0.5 0.5],...
            'linewidth',2)
        h = boxplot(group_data, group_idx);
        set(h,'linewidth',2)
        for iH = 1:size(h,2)
            set(h(6,iH), 'linewidth', 4)
            set(h(:,iH), 'color', use_colors{iH})
            set(h(:,iH), 'MarkerEdgeColor', use_colors{iH})
        end
        
        hp = plotSpread({group_data},'binWidth',0.1,...
                'distributionColors',{[0.7 0.7 0.7]},...
                'distributionIdx',group_idx,'spreadWidth',0.33);
        set(hp{1},'MarkerSize',14)
        
        if options.show_stars
            ax = axis;
            plot([1 3],repmat(1.03*max(group_data),[1 2]),'k-','linewidth',2)
            plot([2], 1.45*max(group_data) ,'kp','MarkerFaceColor','k',...
                'linewidth',2)
            axis([ax(1) ax(2) ax(3) 1.06*max(group_data)])
        end
            
        % Plot the stats on the graph below the title - KWK 20201216
        title({sprintf('%s%s','\bf ',options.which_metab{iM}),...
            sprintf('%s%s%s%.4f%s%.4f','\rm ',this_metab.stats.kruskallwallis.table{1,5},' = ',this_metab.stats.kruskallwallis.table{2,5},...
            '; p = ',this_metab.stats.kruskallwallis.table{2,6})})
        set(gca,'XTick',[1 2 3],'XTickLabel',{[g1_short ', n = ' num2str(sum(...
            ~isnan(this_metab.g1_data)))], [g2_short ', n = ' num2str(sum(...
            ~isnan(this_metab.g2_data)))], [g3_short ', n = ' num2str(sum(...
            ~isnan(this_metab.g3_data)))]})
        if strcmp(options.which_metab{iM},'MacY')
            use_y_label = 'MacM (arb. units)';
        else
            use_y_label = [options.which_metab{iM} ' (mM)'];
        end
        ylabel(use_y_label,'Color','k');
        set(gcf,'color','w','POS',[95   185   750   420])
        box off
        h_leg = legend('Median','25-75%','Range','Data');
        set(h_leg,'Location','northeast','textcolor','k')
        set(gca,'XColor','k','YColor','k','Fontsize',18)
        
        % Plot the stats on the graph - KWK 20201216
        
    end
end
% FDR correction
[sort_p, sort_idx] = sort(all_p_to_FDR);

FDR_p = sort_p .* [numel(options.which_metab):-1:1];
FDR_p(FDR_p > 1) = 1;

for iM = 1:numel(options.which_metab)
    output.(options.which_metab{iM}).stats.kruskallwallis.p_FDR = FDR_p(sort_idx == iM);
end

%% check if metabolites correlate, if requested
if options.corr_metabs
    for iM1 = 1:(numel(options.which_metab)-1)
        for iM2 = (iM1+1):numel(options.which_metab) % don't repeat
            corr_M1 = output.(options.which_metab{iM1}).data;
            corr_M2 = output.(options.which_metab{iM2}).data;
                        
            use_metab_corr_idx = ~isnan(corr_M1) & ~isnan(corr_M2);
            
            [r, p] = corr(corr_M1(use_metab_corr_idx), ...
                corr_M2(use_metab_corr_idx),...
                'type','Spearman');
            df = sum(use_metab_corr_idx)-2;
            
            eval(['output.corr_metab.' options.which_metab{iM1} '_' ...
                options.which_metab{iM2} '.r = r;']);
            eval(['output.corr_metab.' options.which_metab{iM1} '_' ...
                options.which_metab{iM2} '.p = p;']);
            eval(['output.corr_metab.' options.which_metab{iM1} '_' ...
                options.which_metab{iM2} '.df = df;']);
            
            if options.displayFigs
                figure; hold on
                
                [poly_fit] = polyfit(corr_M1(use_metab_corr_idx), ...
                    corr_M2(use_metab_corr_idx), 1);
                
                fit_x = [min(corr_M1(use_metab_corr_idx)) max(corr_M1(use_metab_corr_idx))];
                fit_y = poly_fit(1).*fit_x + poly_fit(2);
                y_range = [min(corr_M2(use_metab_corr_idx)) max(corr_M2(use_metab_corr_idx))];
                plot(fit_x,fit_y,'k-','linewidth',2)
                
                plot(corr_M1(use_metab_corr_idx), corr_M2(use_metab_corr_idx), 'ko', ...
                    'linewidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w')
                
                text(fit_x(1),y_range(2) - (y_range(2) - y_range(1))*0,...
                    ['n = ' num2str(sum(use_metab_corr_idx))],'fontsize',18)
                text(fit_x(1),y_range(2) - (y_range(2) - y_range(1))*0.1,...
                    ['r = ' num2str(round(100*r)/100)],'fontsize',18)
                text(fit_x(1),y_range(2) - (y_range(2) - y_range(1))*0.2,...
                    ['p = ' num2str(round(100*p)/100)],'fontsize',18)
                
                ylabel([options.which_metab{iM2} ' (mM)'],'Color','k');
                xlabel([options.which_metab{iM1} ' (mM)'],'Color','k');
                set(gcf,'color','w','POS',[357   275   560   560])
                box off
                set(gca,'XColor','k','YColor','k','Fontsize',18)
            end
        end
    end
end
%% check E/I, if we're looking at Glu & GABA...
if sum(strcmp('Glu',options.which_metab)) && sum(strcmp('GABA',options.which_metab))
    EI.g1_data = output.Glu.g1_data ./ output.GABA.g1_data;
    EI.g3_data = output.Glu.g3_data ./ output.GABA.g3_data;
    EI.g2_data = output.Glu.g2_data ./ output.GABA.g2_data;
    
    [h p ci EI.stats.ttest2_g1_g2] = ttest2(EI.g1_data, ...
        EI.g3_data);
    EI.stats.ttest2_g1_g2.p = p;
    
    showKW = 'off';
    
    group_data = [EI.g1_data ; EI.g2_data ; ...
        EI.g3_data];
    group_idx = [ones(numel(EI.g1_data),1) ; ...
        2*ones(numel(EI.g2_data),1) ; ...
        3*ones(numel(EI.g3_data),1)];
    
    [p EI.stats.kruskallwallis.table,...
        EI.stats.kruskallwallis.stats] = kruskalwallis(...
        group_data, group_idx, showKW);
    EI.stats.kruskallwallis.p = p;
    
    output.EI = EI;
    
%% plot, if requested
    if options.displayFigs
        figure; hold on
        pause(1) % let this finish plotting...
        
        plot([0 0],repmat(this_metab.g1_median,[1 2]),'-','color',[0.5 0.5 0.5],...
            'linewidth',4)
        plot([0 0],repmat(this_metab.g1_median,[1 2]),'-','color',[0.5 0.5 0.5],...
            'linewidth',2)
        plot([0 0],repmat(this_metab.g1_median,[1 2]),'--','color',[0.5 0.5 0.5],...
            'linewidth',2)
        h = boxplot(group_data, group_idx);
        set(h,'linewidth',2)
        for iH = 1:size(h,2)
            set(h(6,iH), 'linewidth', 4)
            set(h(:,iH), 'color', use_colors{iH})
            set(h(:,iH), 'MarkerEdgeColor', use_colors{iH})
        end
        
        hp = plotSpread({group_data},'binWidth',0.1,...
            'distributionColors',{[0.7 0.7 0.7]},...
            'distributionIdx',group_idx,'spreadWidth',0.33);
        set(hp{1},'MarkerSize',14)
        
        title({'Glu / GABA',...
            sprintf('%s%s%s%.4f%s%.4f','\rm ',output.EI.stats.kruskallwallis.table{1,5},' = ',output.EI.stats.kruskallwallis.table{2,5},...
            '; p = ',output.EI.stats.kruskallwallis.table{2,6})})
        set(gca,'XTick',[1 2 3],'XTickLabel',{[g1_short ', n = ' num2str(numel(...
            EI.g1_data))], [g2_short ', n = ' num2str(numel(...
            EI.g2_data))], [g3_short ', n = ' num2str(numel(...
            EI.g3_data))]})
        ylabel(['Glu / GABA ratio'],'Color','k');
        set(gcf,'color','w','POS',[95   185   750   420])
        box off
        h_leg = legend('Median','25-75%','Range','Data');
        set(h_leg,'Location','Best','textcolor','k')
        set(gca,'XColor','k','YColor','k','Fontsize',18)
    end
end

%% out
output.subj_date = MRS_subj_date;
output.options = options;
output.date_run = datestr(now);
end