% selection of data

% find metab spectra
d=get_subdir_regex('/home/shaw-raid1/data/MRS'),'^20141104-ST001-gosia_AS_E-017-D','metab');
% find water spectra
dw=get_subdir_regex('/home/shaw-raid1/data/MRS'),'^20141104-ST001-gosia_AS_E-017-D','w1');


%d=get_subdir_regex('/home/range2-raid2/aging-data/antiox/young','^20','metab');
%dw=get_subdir_regex('/home/range2-raid2/aging-data/antiox/young','^20','w1');

% reading in the data and organizing 
f=explore_spectro_data(char(d));
% f=concatenate_fid(f);

fw=explore_spectro_data(char(dw));

% initial processing
% eddy current compensation
f_ecc=eddycorrection(f,fw);
fw_ecc=eddycorrection(fw,fw);

% remove first point and add extra 0 at the end (to have the same number of
% points as initially
fo=remove_first_points_fillup(f_ecc,1);

% setting up processing parameters
par=processing_spec;
par.mean_line_broadening=4;
par.method='correlation';
fc=processing_spec(fo,par);

plot_spectrum(fo,par) % before processing
plot_spectrum(fc,par) % after processing
% to know if frequency and phase correction were successful, look for
% variability in the peak freq of NAA (~2ppm - for freq) and variability in
% the shape (e.g., weird dips) on the flanks of the NAA peak (~1.95ppm -
% for phase).

%par.method='max';
%par.correct_freq_mod='abs';
%fc2=processing_spec(fo,par);

% to obtain water linewidth
fw_ecc=get_water_width(fw_ecc);
fw_ecc.water_width_no_cor % gosia does not trust the correction, so don't use it
% (already have done eddy current compensation) - will be a bit smaller
% than the linewidth of water we measured at the scanner

% preparing for LCModel analysis:

pp.root_dir='/home/shaw-raid1/data/MRS//processed_data/OCC';
% it uses fc.sujet_name and fc.SerDescr for naming RAW files
pp.gessfile=5;

write_fid_to_lcRAW(fc,pp,fw_ecc);

