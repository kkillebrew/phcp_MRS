function [output] = decode_mrs( options )
% usage: [output] = decode_mrs( options )
%
%
%
% mps 20201004
%
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
    % options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20190913_phcp_PFC_83subj_H2O_scaled.csv';
end
if ~isfield(options,'mrs_n_col')
    options.mrs_n_col = 503; % 68 if using tissue data without correlations, 507 = if using Gosia's notes, 439 = LCM default
end
if ~isfield(options,'mrs_header_lines')
    options.mrs_header_lines = 6; % Works for KWK
%     options.mrs_header_lines = 6; % 6 = if using Gosia's notes, 8 = LCM default
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
%     options.which_metab = {'Glu','Gln','GABA','NAA','GSH','NAAG'};   % 6 aprior metabolites w/ support in the lit
%         options.which_metab = {'MacY','Asc','GPC','Cr','PCr','GABA','Gln','Glu','GSH','Ins','NAA','NAAG','PE','Tau'}; % CRLB < 20
    %     options.which_metab = {'MacY','Asc','GPC','Cr','PCr','GABA','Gln','Glu','GSH','Ins','NAA','NAAG','PE','Tau','fgray','fwhite','fcsf'}; % CRLB < 20
%     options.which_metab = {'MacY','Asc','Asp','PCho','GPC','Cr','PCr','GABA','Glc','Gln','Glu','GSH','Ins','Lac','NAA','NAAG','PE','sIns','Tau'}; % all
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
    options.quality_file = []; 
    % e.g., options.quality_file = '/home/shaw-raid1/data/MRS/processed_data/data_quality/data_coil_parameters_performance_OCC_20200420.csv';
    % e.g., options.quality_file = '/home/shaw-raid1/data/MRS/processed_data/data_quality/data_coil_parameters_performance_PFC_20200420.csv';
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
if ~isfield (options, 'subj_group_def')
    options.subj_group_def = 1; % 1 = controls, relatives, probands;
    % 2 = controls, SZ, BP
    % 3 = SZ, schizoaffective (SCA), BP;
    % 4 = healthy (con+rel), SZ+SCA, bipolar,
    % 5 = controls, probands, relatives (flip order of P & R)
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
    
    
    for iM = 1:numel(options.which_metab)
        metab = options.mrs_struct.(options.which_metab{iM})(1:end-2); % toss last 2 values - mean and sd/mean
        metab(toss_idx) = []; % toss excluded subj
%         mrs_file
        avg_metab = -1*ones(size(metab));
        avg_metab(unique_idx) = metab(unique_idx);
        avg_metab(repeat_scan_idx-1) = nanmean([metab(repeat_scan_idx-1) ...
            metab(repeat_scan_idx)],2); % avg repeated subj
        metab = avg_metab(~(avg_metab < 0)); % do it this way, to deal with NaNs for tossed values with CRLB > limit
        
        data.(options.which_metab{iM}) = metab;
    end
    
else
    for iM = 1:numel(options.which_metab)
        data.(options.which_metab{iM}) = options.mrs_struct.(options.which_metab{iM})(1:end-2); % toss last 2 values - mean and sd/mean
    end
end

%% Subject group indices (from refit_COP_dataRK) added by MPS 13 NOV 2019
subj_number = MRS_subj_date(:,1);
date_number = MRS_subj_date(:,2);

demog_opts.target_file = '/home/shaw-raid1/data/7T/demographics/PHCP7TfMRIDemo.csv';
demog_opts.subj_number = subj_number;
demog_opts.date_number = date_number;
demog_data = read_in_demog_data(demog_opts);


dx_list = nan(numel(subj_number),1);
missing_dx_list = [];
for iSubj = 1:numel(subj_number)MRS_subj_date(:,1) < 6000000 
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

%% do decoding...
n_chunks = 5; % n-fold cross validation
n_reps = 1000; % repeat m times, with different random sub-sampling
        
n_test_g1 = floor(numel(g1_idx)/n_chunks);
n_test_g2 = floor(numel(g2_idx)/n_chunks);
n_test_g3 = floor(numel(g3_idx)/n_chunks);

train_n_subj_min = min([numel(g1_idx) - n_test_g1 ...
    numel(g2_idx) - n_test_g2 ...
    numel(g3_idx) - n_test_g3]); % use same number of subjects in each group,
                                 % based on this minimum -- randomly
                                 % sub-sample subj within group for 
                                 % training

accuracy = zeros(n_reps, n_chunks);
save_predicted_labels = zeros(n_reps, n_chunks, [n_test_g1 + n_test_g2 + n_test_g3]);
accuracy_shuffle = zeros(n_reps, n_chunks);
save_predicted_labels_shuffle = zeros(n_reps, n_chunks, [n_test_g1 + n_test_g2 + n_test_g3]);
save_test_labels = zeros(n_reps, n_chunks, [n_test_g1 + n_test_g2 + n_test_g3]);

parfor iRep = 1:n_reps
    shuffle_g1_idx = g1_idx(randperm(numel(g1_idx)));
    shuffle_g2_idx = g2_idx(randperm(numel(g2_idx)));
    shuffle_g3_idx = g3_idx(randperm(numel(g3_idx)));
    
    for iChunk = 1:n_chunks
        
        test_g1 = shuffle_g1_idx((1:n_test_g1) + (iChunk-1)*n_test_g1); % N x cross validation
        test_g2 = shuffle_g2_idx((1:n_test_g2) + (iChunk-1)*n_test_g2);
        test_g3 = shuffle_g3_idx((1:n_test_g3) + (iChunk-1)*n_test_g3);
        
        train_g1 = shuffle_g1_idx(~ismember(shuffle_g1_idx, test_g1)); % find g1_idx that are not in test
        train_g2 = shuffle_g2_idx(~ismember(shuffle_g2_idx, test_g2));
        train_g3 = shuffle_g3_idx(~ismember(shuffle_g3_idx, test_g3));
        
        rand_idx_g1 = randperm(numel(train_g1));
        rand_idx_g1 = rand_idx_g1(1:train_n_subj_min); % use the same number of subjects in each group for training...
        train_g1 = train_g1(rand_idx_g1);
        
        rand_idx_g2 = randperm(numel(train_g2));
        rand_idx_g2 = rand_idx_g2(1:train_n_subj_min); % use the same number of subjects in each group for training...
        train_g2 = train_g2(rand_idx_g2);
        
        rand_idx_g3 = randperm(numel(train_g3));
        rand_idx_g3 = rand_idx_g3(1:train_n_subj_min); % use the same number of subjects in each group for training...
        train_g3 = train_g3(rand_idx_g3);
        
        test_data = [];
        train_data = [];
        
        for iM = 1:numel(options.which_metab)
            
            test_data = [test_data ; data.(options.which_metab{iM})(test_g1)' ...
                data.(options.which_metab{iM})(test_g2)' ...
                data.(options.which_metab{iM})(test_g3)' ];
            
            train_data = [train_data ; data.(options.which_metab{iM})(train_g1)' ...
                data.(options.which_metab{iM})(train_g2)' ...
                data.(options.which_metab{iM})(train_g3)' ];
            
        end
        
        
        test_labels = [ ones([numel(test_g1) 1]) ; ...
            2*ones([numel(test_g2) 1]) ; 3*ones([numel(test_g3) 1]) ];
        
        train_labels = [ones([numel(train_g1) 1]) ; ...
            2*ones([numel(train_g2) 1]) ; 3*ones([numel(train_g3) 1]) ];
        
        train_labels_shuffle = train_labels(randperm(numel(train_labels)));
        
        mdl = fitcecoc(train_data', train_labels, 'Coding', 'onevsall', 'Learners', 'SVM' ); % train support vector machine
        mdl_shuffle = fitcecoc(train_data', train_labels_shuffle, 'Coding', 'onevsall', 'Learners', 'SVM' ); % train support vector machine using shuffled labels
        
        predicted_labels = predict(mdl, test_data'); % predict conditions for new data
        predicted_labels_shuffle = predict(mdl_shuffle, test_data'); % predict conditions for new data w/ shuffled labels
        
        % Save iterations 
        save_predicted_labels(iRep,iChunk,:) = predicted_labels;
        save_test_labels(iRep,iChunk,:) = test_labels;
        save_predicted_labels_shuffle(iRep,iChunk,:) = predicted_labels_shuffle;
        
        accuracy(iRep,iChunk) = mean(predicted_labels == test_labels);
        accuracy_shuffle(iRep,iChunk) = mean(predicted_labels_shuffle == test_labels);
    end
end

output.predicted_labels = save_predicted_labels;
output.accuracy = accuracy;
output.predicted_labels_shuffle = save_predicted_labels_shuffle;
output.accuracy_shuffle = accuracy_shuffle;

output.p_val_boot = sum( sum( accuracy_shuffle > accuracy ) ) / (n_reps * n_chunks);

% Compre the two distributions using K-S test % KWK - 20201228
% https://en.m.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
[output.ks_test.h, output.ks_test.p, output.ks_test.stats] = ...
    kstest2(nanmean(output.accuracy,2), nanmean(output.accuracy_shuffle,2));

if options.displayFigs
    % Plot the two distributions % KWK - 20201228
    figure()
    hist([nanmean(output.accuracy_shuffle,2) nanmean(output.accuracy,2)],20)
    set(gcf,'color','w')
%     ylabel('')
%     xlim([-10 n_reps+10])
%     xlabel(['Repetitions (n=' num2str(n_reps) ')'])
    box off
    
    figure; 
    subplot(1,4,1:3)
    errorbar(mean(accuracy,2),std(accuracy,0,2)/sqrt(n_chunks),'k-')
    hold on;
    plot([-10 n_reps+10],[0.33 0.33],'r--'); 
    set(gcf,'color','w')
    ylabel('Decoding accuracy (%)')
    xlim([-10 n_reps+10])
    xlabel(['Repetitions (n=' num2str(n_reps) ')'])
    box off
    
    subplot(1,4,4)
    bar(mean(mean(accuracy,2)));
    hold on
    plot([.5 1.5],[0.33 0.33],'r--');
    errorbar(mean(mean(accuracy,2)),std(mean(accuracy,2)),'k-');
    set(gcf,'color','w')
    ylabel('Average decoding accuracy (%)')
    ylim([0 .55])
    box off
end

%% confusion
confusion_bins = zeros(3, 3);
for iCond = 1:3 % test, true
    for iCond2 = 1:3 % predicted
        find_confusion = ( save_test_labels == iCond ) & ...
            ( save_predicted_labels == iCond2);
        confusion_bins(iCond,iCond2) = sum(find_confusion(:));
    end
end
output.confusion_matrix = confusion_bins./repmat(sum(confusion_bins,2),...
    [1 size(confusion_bins,1)]);

% Set up the c axis min max - KWK 20210421
cMin = 0;
cMax = .65;   % Hardcoded as the current max value of the 4 decoding conditions - KWK 20210421

if options.displayFigs
    labels = {g1_short, g2_short, g3_short};
    figure
    imagesc(output.confusion_matrix, [0 max(output.confusion_matrix(:))])
    caxis([cMin, cMax]);
    colorbar
    caxis([cMin, cMax]);
    set(gcf,'color','w')
    box off
    xlabel('Predicted')
    ylabel('True')
    set(gca,'XTick',1:3,'XTickLabel',labels,'YTick',1:3,'YTickLabel',labels)
    % ('Accuracy')
end
%% out

output.subj_date = MRS_subj_date;
output.data = data;
output.options = options;
output.date_run = datestr(now);

% Save
clockHolder = clock;
save(sprintf('%s%04.f%02.f%02.f','./mrsDecode_',...
    clockHolder(1),clockHolder(2),clockHolder(3)),'-struct','output');

end