%% Data analysis for Helium forbidden transition spectroscopy

% clear all;
% Remove old data dirs from path

data = [];


%% set up path
% find this .m file's path, this must be in the project root dir
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));%add all subfolders to the path to find genpath_exclude
path_to_genpath=fileparts(which('genpath_exclude'));
path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(this_folder)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(this_folder,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(this_folder,'dev'),'\.'))
%addpath(genpath_exclude(fullfile(this_folder,'bin'),'\.'))

%% variables
% add all subfolders to the path
% % Setting up

anal_opts=[]; %reset the options (would be good to clear all variables except the loop config
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190704_forbidden_long interrogation\';
anal_opts.tdc_import.save_cache_in_data_dir=true;
tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
tlim=[0,inf];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];
anal_opts.dld_aquire=10;
anal_opts.trig_dld=20.3;
anal_opts.probe.t0=20;
anal_opts.probe.duration=6;


cli_format_text('','c',2)
cli_format_text('STARTING ANALYSIS','c',2)
cli_format_text('','c',2)

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import TDC files
%add a file seperator to the end of the import path
if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end

import_opts.shot_num=find_data_files(anal_opts.tdc_import);

%set up an output dir %https://gist.github.com/ferryzhou/2269380
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir(fullfile(anal_opts.tdc_import.dir,'out')); end
%make a subfolder with the ISO timestamp for that date
anal_out.dir=sprintf('%sout\\%s\\',...
    anal_opts.tdc_import.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;
%start up the diary of stdout
diary([anal_out.dir,'anal.txt'])
%import the data

anal_opts.tdc_import.force_reimport=true;
data.mcp_tdc = import_mcp_tdc_data(anal_opts.tdc_import);
cli_format_text('Data import complete!','h',1)

%% Import LabView log
import_opts.labview.dir=anal_opts.tdc_import.dir;
import_opts.labview.out_dir=anal_opts.global.out_dir;
import_opts.labview.plots=true
data.labview = import_labview_log(import_opts);

data=match_labview_log(data,anal_opts)


%% Import the analog log files
%TODO
% modify to return the integrated pd signal


% anal_opts.ai_log.dir=anal_opts.tdc_import.dir;
% anal_opts.ai_log.force_reimport=false;
% anal_opts.ai_log.force_load_save=false;
% anal_opts.ai_log.log_name='log_analog_in_';
% %because im only passing the ai_log feild to aviod conflicts forcing a reimport i need to coppy these feilds
% anal_opts.ai_log.calibration=data.mcp_tdc.probe.calibration;
% anal_opts.ai_log.pd.set_probe=anal_opts.probe_set_pt;
% anal_opts.ai_log.trig_dld=anal_opts.trig_dld;
% anal_opts.ai_log.dld_aquire=anal_opts.dld_aquire;
% anal_opts.ai_log.aquire_time=anal_opts.dld_aquire;
% anal_opts.ai_log.trig_ai_in=anal_opts.trig_ai_in;
% % set time matching conditions
% anal_opts.ai_log.aquire_time=4;
% anal_opts.ai_log.pd.diff_thresh=0.05;
% anal_opts.ai_log.pd.std_thresh=0.05;
% anal_opts.ai_log.pd.time_start=0.2;
% anal_opts.ai_log.pd.time_stop=2;
% anal_opts.ai_log.time_match_valid=8; %how close the predicted start of the shot is to the actual
% %sfp options
% anal_opts.ai_log.scan_time=1/20; %fast setting 1/100hz %estimate of the sfp scan time,used to set the window and the smoothing
% anal_opts.ai_log.sfp.num_checks=inf; %how many places to check that the laser is single mode, inf=all scans
% anal_opts.ai_log.sfp.peak_thresh=[-0.005,-0.005];%[0,-0.008]*1e-3; %theshold on the compressed signal to be considered a peak
% anal_opts.ai_log.sfp.pzt_dist_sm=4.5;%minimum (min peak difference)between peaks for the laser to be considered single mode
% anal_opts.ai_log.sfp.pzt_peak_width=0.15; %peak with in pzt voltage used to check that peaks are acually different and not just noise
% anal_opts.ai_log.plot.all=false;
% anal_opts.ai_log.plot.failed=true;
% 
% %do the ac waveform fit
% anal_opts.ai_log.do_ac_mains_fit=false;
% 
% data.ai_log=ai_log_import(anal_opts.ai_log,data);
%%

anal_opts.wm_log.dir=anal_opts.tdc_import.dir;
anal_opts.wm_log.force_reimport=false;
wm_log_name='log_wm_';
wm_logs=dir(fullfile(anal_opts.wm_log.dir,cat(2,wm_log_name,'*.txt')));
anal_opts.wm_log.names={wm_logs.name};
data.wm_log.raw=wm_log_import(anal_opts.wm_log);

%% process
anal_opts.wm_log.plot.all=false;
anal_opts.wm_log.plot.failed=false;
anal_opts.wm_log.force_reimport=false;

anal_opts.wm_log.time_pd_padding=4; %check this many s each side of probe
anal_opts.wm_log.time_blue_padding=1; %check this many seconde each side of probe
anal_opts.wm_log.time_probe=6;
anal_opts.wm_log.ecd_volt_thresh=0.5;

anal_opts.wm_log.red_sd_thresh=5; %allowable standard deviation in MHz
anal_opts.wm_log.red_range_thresh=10; %allowable range deviation in MHz
anal_opts.wm_log.rvb_thresh=20; %allowable value of abs(2*red-blue)

anal_opts.wm_log.global=anal_opts.global;

data.wm_log.proc=wm_log_process(anal_opts,data);
clear('sub_data')



%% Match the timestamps    
%% data.sync = match_timestamps(data,import_opts);
%% bryce change, this is now in the wavemeter proceesing code


%%
import_opts.signal=[];
import_opts.signal.plot.lim.x=[-45,45]*1e-3;
import_opts.signal.plot.lim.y=[-45,45]*1e-3;
import_opts.signal.plot.nbins=1e3;
import_opts.signal.plot.blur=3;
import_opts.signal.plot.cmp_dyn_range=true;
tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
%tlim=[0.5,6.2];
tlim=[5,22.5];
import_opts.signal.square_mask=[tlim;tmp_xlim;tmp_ylim];
import_opts.signal.circ_mask=[[0,0,35e-3,1];
                              [35e-3,5e-3,7e-3,0];
                              [25.05e-3,-19e-3,8e-3,0];
                              [26.94e-3,21.98e-3,5e-3,0];
                              [19.1e-3,-28.0e-3,4e-3,0];
                              [2.973e-3,-33.15e-3,4e-3,0];
                              [6.216e-3,34.41e-3,3e-3,0];
                              [13.78e-3,31.62e-3,3e-3,0];
                              [21.26e-3,25.95e-3,3e-3,0];
                              [30.36e-3,-12.52e-3,5e-3,0];
                              ];

data.signal.masked_number = forbidden_signal_masked_num(data,import_opts.signal);





%% Generate signal
data.signal.total_num = signal_process(data,import_opts);



% Create a calibration model
data.cal = make_calibration_model(data,import_opts);

%% Mask out shots which failed
% data.check = check_for_errors(data,opts);

% %% Break data into categories
% data.cat = categorize_shots(data,opts);
 
% %% Peak detection
%data = auto_peak_detect(data,opts);

% %% Fit the detected peaks
%data = fit_detected_peaks(data,opts);
% data = fancy_fits(data,opts);
% 
% %% Zeeman shift correction
%data = zeeman_correction(data,opts);

%% Presentation plots
%data = present_plots(data,opts);

%     
% 
% end
% % % RESULTS
% Transition name
% Directory
% Various stat errors
%%
% if iscell(opts.e_state)
%     num_pks=numel(opts.e_state);
% else
%     num_pks = 1;
% end
% for pidx=1:num_pks
%         if ~iscell(opts.e_state)
%             e_state = opts.e_state;
%         else
%             e_state = opts.e_state{pidx};
%         end
%         e_level = e_state(1:6);
%         fmt_name = strrep(e_level,'^','_');
%         f_pred = opts.const.f_table.g_2_3P_2.(sprintf('e_%s',fmt_name))/1e6; %MHz
%     fitted_freqs = cellfun(@(x) x.zeeman.corrected(pidx), data.cat);
%     stat_err = cellfun(@(x) x.zeeman.stat_unc(pidx), data.cat);
%     fprintf('===COMBINED MEASUREMENT===\n')
%     fprintf('%s measured value:  %.3f(%.3f) MHz\n',e_level,mean(fitted_freqs),sum(stat_err))
%     fprintf('Theory difference:  %.3f MHz\n',mean(fitted_freqs)-f_pred)
% end
% 
% cli_header({1,'Done.'})
% 
% 
% 
% % %% Save to output
% cli_header({0,'Saving output...'})
% out_data.data = data.cat;
% out_data.options = opts;
% filename = fullfile(opts.out_dir,'output_and_options.mat');
% save(filename,'out_data','-v7.3')
% fwtext('All Done!')