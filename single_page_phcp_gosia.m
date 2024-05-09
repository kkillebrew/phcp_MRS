% to add romain's code to paths - added by mps 20180816
addpath(genpath('/home/shaw-raid1/matlab_tools/matspec/'))
% addpath(genpath('/home/range4-raid1/gosia/matlab/from_Romain_151126/matspec/spect_processing'))

% selection of data

d=get_subdir_regex('/home/shaw-raid1/data/MRS/dicom_data/','^P','metab_RES');
dw=get_subdir_regex('/home/shaw-raid1/data/MRS/dicom_data/','^P','w1$');


%d=get_subdir_regex('/home/range2-raid2/aging-data/antiox/young','^20','metab');
%dw=get_subdir_regex('/home/range2-raid2/aging-data/antiox/young','^20','w1');

% reading in the data and organizing 
f=explore_spectro_data(char(d));
f=concatenate_fid(f);

fw=explore_spectro_data(char(dw));

% initial processing
% eddy current compensation
f_ecc=eddycorrection(f,fw);
fw_ecc=eddycorrection(fw,fw);

% remove first point and add extra 0 at the end (to have the same number of
% points as initially
fo=remove_first_points_fillup(f_ecc,1);
%fo=f_ecc;

% setting up processing parameters
par=processing_spec;
par.mean_line_broadening=4;
par.method='correlation';
fc=processing_spec(fo,par);

plot_spectrum(fo,par)
plot_spectrum(fc,par)

% par.method='max';
% par.correct_freq_mod='abs';
% fc2=processing_spec(fo,par);

% to obtain water linewidth
fw_ecc=get_water_width(fw_ecc);
fw_ecc.water_width_no_cor 

% preparing for LCModel analysis:

pp.root_dir='/home/rosalind2-raid4/aging-data/MRSofDMN/LCModel/data';
% it uses fc.sujet_name and fc.SerDescr for naming RAW files
pp.gessfile=5;

%write_fid_to_lcRAW(fc,pp,fw_ecc);

