% selection of data

d=get_subdir_regex('/home/range2-raid2/aging-data/antiox/elderly','^20131119-ST001-gosia_AS_E-007','metab');
dw=get_subdir_regex('/home/range2-raid2/aging-data/antiox/elderly','^20131119-ST001-gosia_AS_E-007','w1');


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

% setting up processing parameters
par=processing_spec;
par.mean_line_broadening=4;
par.method='max';
par.figure=0;

% frequency and phase correction
par.correct_freq_mod='abs';
fc=processing_spec(fo,par);
par.correct_freq_mod='abs';
fc2=processing_spec(fc,par);
par.correct_freq_mod='real';
fc3=processing_spec(fc2,par);

plot_spectrum(fc3,par)

% to obtain water linewidth
fw_ecc=get_water_width(fw_ecc);

% fw_ecc.water_width_no_cor gives the linewidth of water

fw_ecc.water_width_no_cor 

pp.root_dir='/home/range2-raid2/aging-data/antiox/LCModel';
% it uses fc.sujet_name and fc.SerDescr for naming RAW files
pp.gessfile=5;

write_fid_to_lcRAW(fc,pp,fw_ecc);


