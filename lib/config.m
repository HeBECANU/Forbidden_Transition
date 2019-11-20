function [opts,const] = master_transition_config(data_dir)


if ~strcmp(data_dir(end),filesep)
    data_dir = [data_dir,filesep];
end

%set up an output dir %https://gist.github.com/ferryzhou/2269380
opts.dir = data_dir;
if (exist(fullfile(opts.dir,'out'), 'dir') == 0), mkdir(fullfile(opts.dir,'out')); end
%make a subfolder with the ISO timestamp for that date
opts.out_dir=sprintf('%sout\\%s\\', opts.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
if (exist(opts.out_dir, 'dir') == 0), mkdir(opts.out_dir); end
opts.cache_dir=fullfile(opts.dir,'cache');
if (exist(opts.cache_dir, 'dir') == 0), mkdir(opts.cache_dir); end
%% Frequently adjusted quantities

% These variables are set here for quick access to quantities passed to
% functions in fields defined later in this code
opts.ignorefiles = nan; %Index of shot after import (i.e. chronological position, NOT TDC number
opts.freq_bin_size = 1;

% Experimental parameters
opts.probe_set_pt=0.4;
opts.trig_ai_in=20;
opts.trig_dld=20.3;
opts.dld_aquire=10;
opts.aquire_time=10;
opts.aom_freq=189;%190*1e6;%Hz %set to zero for comparison with previous data runs
% Experimental constants
%% Oscillations
opts.osc_fit.binsx=1000;
opts.osc_fit.blur=1;
opts.osc_fit.xlim=[-20,20]*1e-3;
opts.osc_fit.tlim=[0.86,1.08];
opts.osc_fit.dimesion=2; %Sel ect coordinate to bin. 1=X, 2=Y.
opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq
opts.osc_fit.appr_osc_freq_guess=[52,52,52];
opts.osc_fit.freq_fit_tolerance=5;
opts.osc_fit.plot_fits=false;
opts.osc_fit.plot_err_history=false;
opts.osc_fit.plot_fit_corr=false;
%% LabView import

opts.lv.plots = true;
opts.lv.dir = opts.dir;
opts.lv.out_dir = opts.out_dir;
%% Analog import
%opts.ai.num_files = 500;%nan

% opts.ai.force_recalc_import = true;
opts.ai.log_name='log_analog_in_';
opts.ai.verbose = 1;
opts.ai.plots=1;

% Default values, should be overriden by local
opts.ai.num_files = nan; %nan
opts.ai.pd_offset = -0.065;
opts.ai.high_thresh = 8e-3;
opts.ai.srate = 20000.00;
opts.ai.exposure = 0.95*0.15;
opts.ai.pd_setpoint = .043;

opts.ai.post_fun = @transition_post_ai;
% Options for global import
opts.ai.cache_import.verbose=0;
opts.ai.cache_import.force_cache_load=false;
opts.ai.cache_import.force_recalc=false;
opts.ai.cache_import.dir = opts.cache_dir;
% Options for single imports (unnecessary unless they're cached, ill-advised usually)
opts.ai.cache_single_import = false;
opts.ai.cache_single.verbose=0;
opts.ai.cache_single.force_recalc=false;
opts.ai.cache_single.path_directions={1,'dir'};
opts.ai.args_single.cmp_multiplier_disp=50; %multiplier to display the compressed data better

opts.ai.dir=opts.dir;
opts.ai.out_dir = opts.out_dir;



%% Wavemeter log importing

opts.wm.force_reimport=false;
opts.wm.num_logs = nan;
opts.wm.plots = true;
opts.wm.plot_all=1;
opts.wm.plot_failed=false;

opts.wm.wm_log_name='log_wm_';

opts.wm.cache_import.verbose=0;
opts.wm.cache_import.force_recalc=0;
opts.wm.cache_import.dir = opts.cache_dir;
opts.wm.cache_import.save_compressed=true;%needed otherwise save takes a very long time
opts.wm.cache_import.path_directions={1,'dir'};

opts.wm.plot_all=true;
opts.wm.plot_failed=false;
opts.wm.force_reimport=false;




opts.wm.dir=opts.dir;
wm_logs=dir([opts.wm.dir,opts.wm.wm_log_name,'*.txt']);
opts.wm.names={wm_logs.name};
opts.wm.out_dir = opts.out_dir;
%% TDC import


opts.tdc.plots = true;
opts.tdc.file_name='d';
opts.tdc.force_load_save=false;   %takes precidence over force_reimport
opts.tdc.force_reimport=false;
opts.tdc.force_forc=false;
opts.tdc.dld_xy_rot=0.61;
opts.tdc.dir = opts.dir;
opts.tdc.out_dir = opts.out_dir;

opts.tdc.cache_opts.dir = opts.cache_dir;

%Should probably try optimizing these
tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,10];
opts.tdc.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.max_runtime=inf;%inf%cut off the data run after some number of hours, should bin included as its own logic not applied to the atom number ok

%% Atom Laser stuff
opts.atom_laser.pulsedt=8.000e-3;
opts.atom_laser.t0=6.41784; %center i ntime of the first pulse
opts.atom_laser.start_pulse=1; %atom laser pulse to start with
opts.atom_laser.pulses=100;

opts.atom_laser.appr_osc_freq_guess=[52,40,40];
opts.atom_laser.pulse_twindow=opts.atom_laser.pulsedt*0.95;

opts.atom_laser.xylim=opts.tdc.txylim(2:3,:); %set same lims for pulses as import
%% Error handling
opts.check.min_counts = 1e3;
opts.check.wm_tolerance = 200;
opts.check.num_cal_bins = 20;

%% Peak detection
opts.peak.plot = true;
opts.peak.cutoff_thresh = 5e3;
opts.peak.smooth_width = 3;
opts.peak.saturation_threshold = 0.975;
opts.peak.fitwidth = 15;

%% Plotting

opts.num_freq_bins = 30;

%% Physical constants
opts.Bfield = [18.25,11.43]; % Gauss
opts.Bfield_f_unc = [1,1]; %MHz. From memory, this is about 2% unc
% So temp set this manually:
opts.Bfield_unc = opts.Bfield./[51,32]; % Gauss

const.mu = 9.27e-28; %J/G
const.h = 6.63e-34;
const.hbar = const.h/(2*pi);
const.f_mu = const.mu/const.h;
const.w_mu = const.mu/const.hbar;
const.c = 299792458;
const.q = 1.602e-19;
const.g0 = 9.81;
% Notation & lookup
const.terms = {'S','P','D','F','G'};
%% Reference values
const.f_table.g_2_3P_2.e_5_3S_1 = 727.3032446e12;
% Misc transitions - what do the stars mean?
const.f_table.g_2_3P_2.e_5_3P_0 = 1e9*const.c/404.628937550957;
const.f_table.g_2_3P_2.e_5_3P_1 = 1e9*const.c/404.629844755577;
const.f_table.g_2_3P_2.e_5_3P_2 = 1e9*const.c/404.629918705477; 
% Historically controversial transitions
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620968e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622889e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651246e12; 
% Singlet-triplet transitions
const.f_table.g_2_3P_2.e_5_1S_0 = 1e9*const.c/406.8886971706;
const.f_table.g_2_3P_2.e_5_1P_1 = 1e9*const.c/402.322271224483;
const.f_table.g_2_3P_2.e_5_1D_2 = 744.43034335e12; % 402.7nm

%Fitted valuse for the 5^3D's
const.f_table.g_2_3P_2.e_5_3D_3 = 744.39620836e12;
const.f_table.g_2_3P_2.e_5_3D_2 = 744.39622758e12;
const.f_table.g_2_3P_2.e_5_3D_1 = 744.39651114e12;


%% Experimental constants
const.fall_time=0.417;
const.global.qe=0.09;

%% Append to options
opts.const = const;
opts.osc_fit.global=opts.const;
opts.osc_fit.global.velocity=const.g0*const.fall_time;







    
end


