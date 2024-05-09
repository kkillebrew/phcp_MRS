%% find the dicoms on range5-raid1

dicom_dir = '/home/range5-raid1/dicom';
% /home/range5-raid8/hcp/projects/CCF_PHCP_ITK/

all_dicoms = dir(fullfile(dicom_dir,'*PHCPP*7*MRS'));

%% exclude the following data sets from analysis - known issues
exclude_list = {'20180616-ST001-PHCPP6004687_V3_7ZMRS',...
                '20180625-ST001-PHCPP6001501_V3_7ZMRS',...
                '20190525-ST001-PHCPP6010671_V3_7ZMRS',...
                '20190410-ST001-PHCPP6010926_V2_7BMRS',...
                '20190715-ST001-PHCPP101538_V2_7BMRS',...
                '20181204-ST001-PHCPP6010526_V2_7BMRS'};

%% link all
% N.B. this checks whether there is at least 1 'metab' folder in the
% original dicom folder, and if not, no link is created...

MRS_dicom_dir = '/home/shaw-raid1/data/MRS/dicom_data';

for iD = 1:numel(all_dicoms)
    date_str = all_dicoms(iD).name(1:8);
    subj_id = all_dicoms(iD).name(20:27);
    
    source_dir = fullfile(dicom_dir,all_dicoms(iD).name);
    target_dir = fullfile(MRS_dicom_dir,[subj_id '_' date_str] );
    
    find_metab = dir(fullfile(source_dir,'*metab*'));
    
    if ~exist(target_dir,'dir') && ~sum(strcmp(source_dir,fullfile(...
            dicom_dir,exclude_list))) && ~isempty(find_metab)
        
        eval(['! ln -s ' source_dir ' ' target_dir]);
    end
end
