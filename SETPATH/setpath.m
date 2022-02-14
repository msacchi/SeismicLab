%SETPATH: program to set path for SeismicLab and Demos 
 
% --------------------------------------------------------
% You need to modify Dir and DirD according to your system  
% --------------------------------------------------------
 

 Dir='/Users/msacchi/Dropbox/SeismicLab/codes/';
 DirD='/Users/msacchi/Dropbox/SeismicLab/';


 path(path, strcat(Dir,'bp_filter'));
 path(path, strcat(Dir,'decon'));
 path(path, strcat(Dir,'dephasing'));
 path(path, strcat(Dir,'fx_decon'));
 path(path, strcat(Dir,'fxy_eigen'));
 path(path, strcat(Dir,'interpolation'));
 path(path, strcat(Dir,'kl_transform'));
 path(path, strcat(Dir,'linear_operators'));
 path(path, strcat(Dir,'med'));
 path(path, strcat(Dir,'mssa'));
 path(path, strcat(Dir,'mssa'));
 path(path, strcat(Dir,'mwni'));
 path(path, strcat(Dir,'pmf'));
 path(path, strcat(Dir,'pocs'));
 path(path, strcat(Dir,'radon'));
 path(path, strcat(Dir,'rank_reduction'));
 path(path, strcat(Dir,'scaling_tapering'));
 path(path, strcat(Dir,'segy'));
 path(path, strcat(Dir,'seismic_plots'));
 path(path, strcat(Dir,'solvers'));
 path(path, strcat(Dir,'spectra'));
 path(path, strcat(Dir,'synthetics'));
 path(path, strcat(Dir,'utilities'));
 path(path, strcat(Dir,'velan_nmo'));

 path(path, DirD);
 path(path, strcat(DirD,'SeismicLab_demos'));
 path(path, strcat(DirD,'SeismicLab_data'));
