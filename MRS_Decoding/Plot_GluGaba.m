% plot 5 randomly selected patients and plot their GABA and Glu values
% To be run after 'decode_mrs', needs data from that...

function [data,options] = Plot_GluGaba(options)

options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20200414_phcp_OCC_139subj_H2O_scaled_tissue.csv';
options.quality_file = '/home/shaw-raid1/data/MRS/processed_data/data_quality/data_coil_parameters_performance_OCC_20200420.csv';
this = decode_mrs(options);

% Randomly grab 5 patients
partIdxHolder = find(this.subj_date >= 6000000);   % Grab all the patients
data.partIdx = partIdxHolder(randi([1 length(partIdxHolder)],[1 5]));

data.gaba = this.data.GABA(data.partIdx);   % This = output from decode_mrs
data.glu = this.data.Glu(data.partIdx);

data.aveGaba = nanmean(data.gaba);   % This = output from decode_mrs
data.aveGlu = nanmean(data.glu);

% Plot
x_labels = {'GABA','Glutamate'};

addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))

% Reformat for plotting
data.combinedGabaGlu = [data.gaba, data.glu];
data.combinedGrouping = [ones([length(data.gaba) 1]), ones([length(data.gaba) 1])+1];

for i =1:size(data.combinedGabaGlu,2)
    figure()
    figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.switchRate.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.switchRate.figSize = [0 0 ...
        figSize.switchRate.baseSize(4)*...
        (figSize.switchRate.aspectRatio(1)/figSize.switchRate.aspectRatio(2))...
        figSize.switchRate.baseSize(4)];   % Size/postion of fig
    
    % Boxplots
    hb = boxplot(data.combinedGabaGlu(:,i),data.combinedGrouping(:,i));
    set(gca,'XTick',1,'XTickLabel',x_labels{i},'fontsize',15)
    set(hb,'linewidth',2)
    % hb2 = findobj(gca,'type','line');
    % for iHB = 1:size(hb,2)
    %     set(hb2((iHB)+3:3:end),'color',[0 0 1])
    %     set(hb2((iHB)+3:3:end),'MarkerEdgeColor',[0 0 1])
    % end
    hold on
    
    % Beeswarm
    x_val = [1];
    bee_bin_width = .1;
    bee_spread_width = .5;
    beePlot = plotSpread(data.combinedGabaGlu(data.combinedGrouping==i),...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.5 .5 .5]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',30)
    hold on
    
    title('Average GABA and Glutamate In Psychosis','fontsize',18)
    box off
    ylabel('Average Concentration','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.switchRate.figSize)
    
end


end