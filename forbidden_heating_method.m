%% Data analysis for Helium forbidden transition spectroscopy using heating method

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
anal_opts.tdc_import.dir='X:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden427_overnight_heating_method';
anal_opts.tdc_import.save_cache_in_data_dir=true;
tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
tlim=[0,inf];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];
anal_opts.dld_aquire=25;
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
hebec_constants
anal_opts.global.fall_velocity=const.g0*anal_opts.global.fall_time;

%% Import TDC files
%add a file seperator to the end of the import path
if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end

anal_opts.shot_num=find_data_files(anal_opts.tdc_import);

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
anal_opts.labview.dir=anal_opts.tdc_import.dir;
anal_opts.labview.out_dir=anal_opts.global.out_dir;
anal_opts.labview.plots=true
data.labview = import_labview_log(anal_opts);

data=match_labview_log(data,anal_opts)


%% Import the analog log files

use_ai = 1;
if use_ai
anal_opts.ai_log.dir=anal_opts.tdc_import.dir;
anal_opts.ai_log.force_reimport=false;
anal_opts.ai_log.force_load_save=false;
anal_opts.ai_log.log_name='log_analog_in_';
%because im only passing the ai_log feild to aviod conflicts forcing a reimport i need to coppy these feilds
anal_opts.ai_log.calibration=data.mcp_tdc.probe.calibration;
anal_opts.ai_log.pd.set_probe=1.1;%anal_opts.probe_set_pt;
anal_opts.ai_log.trig_dld=anal_opts.trig_dld;
anal_opts.ai_log.dld_aquire=anal_opts.dld_aquire;
anal_opts.ai_log.aquire_time=anal_opts.dld_aquire;
anal_opts.ai_log.trig_ai_in=3.0;%anal_opts.trig_ai_in;
% set time matching conditions
anal_opts.ai_log.aquire_time=21;
anal_opts.ai_log.pd.diff_thresh=0.5;
anal_opts.ai_log.pd.std_thresh=0.5;
anal_opts.ai_log.pd.time_start=0.0;
anal_opts.ai_log.pd.time_stop=21;
anal_opts.ai_log.time_match_valid=16; %how close the predicted start of the shot is to the actual
%sfp options
anal_opts.ai_log.scan_time=1.1; %fast setting 1/100hz %estimate of the sfp scan time,used to set the window and the smoothing
anal_opts.ai_log.sfp.num_checks=inf; %how many places to check that the laser is single mode, inf=all scans
anal_opts.ai_log.sfp.peak_thresh=[0,0];%[0,-0.008]*1e-3; %theshold on the compressed signal to be considered a peak
anal_opts.ai_log.sfp.pzt_dist_sm=4.5;%minimum (min peak difference)between peaks for the laser to be considered single mode
anal_opts.ai_log.sfp.pzt_peak_width=0.5; %peak with in pzt voltage used to check that peaks are acually different and not just noise
anal_opts.ai_log.plot.all=0;
anal_opts.ai_log.plot.failed=0;

%do the ac waveform fit
anal_opts.ai_log.do_ac_mains_fit=false;

data.ai_log=ai_log_import(anal_opts.ai_log,data);
%data.ai_log.pd.mean gives the average pd value
end
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

anal_opts.wm_log.red_sd_thresh=2; %allowable standard deviation in MHz
anal_opts.wm_log.red_range_thresh=4; %allowable range deviation in MHz
anal_opts.wm_log.rvb_thresh=20; %allowable value of abs(2*red-blue)

anal_opts.wm_log.global=anal_opts.global;

data.wm_log.proc=wm_log_process(anal_opts,data);
clear('sub_data')

%% Match the timestamps    
%% data.sync = match_timestamps(data,import_opts);
%% bryce change, this is now in the wavemeter proceesing code

tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
%tlim=[0.5,6.2];
tlim=[0,inf];
anal_opts.hotspot_mask.square_mask=[tlim;tmp_xlim;tmp_ylim];
anal_opts.hotspot_mask.circ_mask=[[0,0,35e-3,1];
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
skip_mask=false;

if ~skip_mask
    num_shots=numel(data.mcp_tdc.counts_txy);
    empty_shots=cellfun(@isempty,data.mcp_tdc.counts_txy);
    data.mcp_tdc.masked.counts_txy={};
    data.mcp_tdc.masked.num_counts=data.mcp_tdc.num_counts*nan;
    fprintf('masking shots %04u:%04u',num_shots,0)
    for ii=1:num_shots
        txy_shot=data.mcp_tdc.counts_txy{ii};
        if ~isempty(txy_shot)
        txy_shot=masktxy_square(txy_shot,anal_opts.hotspot_mask.square_mask);
        txy_shot=masktxy_2d_circle(txy_shot,anal_opts.hotspot_mask.circ_mask);
        data.mcp_tdc.masked.num_counts(ii)=numel(txy_shot);
        data.mcp_tdc.masked.counts_txy{ii}=txy_shot;
        else
            warning('empty shot')
        end
        if mod(ii,10)==0,fprintf('\b\b\b\b%04u',ii),end 
    end
    fprintf('...Done\n') 
else
    data.mcp_tdc.masked.counts_txy=data.mcp_tdc.counts_txy;
    data.mcp_tdc.masked.num_counts=data.mcp_tdc.num_counts;
end



%%
%% warning if you do this with all the counts will likely run out of mem
anal_opts.plot2d.do=false;
if anal_opts.plot2d.do
anal_opts.signal=[];
anal_opts.plot2d.lim.x=[-45,45]*1e-3;
anal_opts.plot2d.lim.y=[-45,45]*1e-3;
anal_opts.plot2d.nbins=1e3;
anal_opts.plot2d.blur=3;
anal_opts.plot2d.cmp_dyn_range=true;


not_empty_shots=~cellfun(@isempty,data.mcp_tdc.masked.counts_txy);
all_txy=cat(1,data.mcp_tdc.masked.counts_txy{not_empty_shots});
dyn_range_pow=0.2;
num_bins=anal_opts.plot2d.nbins;
edge_x=linspace(min(anal_opts.plot2d.lim.x),max(anal_opts.plot2d.lim.x),num_bins);
edge_y=linspace(min(anal_opts.plot2d.lim.y),max(anal_opts.plot2d.lim.y),num_bins);

spatial_blur=anal_opts.plot2d.blur;
bin_area=(range(anal_opts.plot.lim.x)/num_bins)*(range(anal_opts.plot.lim.x)/num_bins);
[counts,centers]=hist3(all_txy(:,2:3),'edges',{edge_x,edge_y});
counts=counts/bin_area;
counts=counts/num_shots;

%imagesc seems to plot the wrong way round so we transpose here

if anal_opts.plot2d.cmp_dyn_range
    counts=counts.^dyn_range_pow;
end
if  ~spatial_blur==0
    counts=imgaussfilt(counts,spatial_blur);
end
stfig('counts during probe')
imagesc(10^3*centers{1},10^3*centers{2},transpose(counts))
colormap(viridis())
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
title('Spatial Dist. TOP')
xlabel('X(mm)')
ylabel('Y(mm)')
h=colorbar;
if anal_opts.plot.cmp_dyn_range
xlabel(h,sprintf('Count Density^{%.2f} (m^{-2})',dyn_range_pow))
else
xlabel(h,'Count Density (m^{-2})')
end
end
    

%% fitting the atom laser pulses
%now find the mean position of each pulse of the atom laser in each shot

anal_opts.atom_laser=[];
% settings for older heating runs
% anal_opts.atom_laser.pulsedt=120e-3;
% anal_opts.atom_laser.t0=1.637327; %center i ntime of the first pulse
% anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
% anal_opts.atom_laser.pulses=189;

anal_opts.atom_laser.pulsedt=240e-3;%120e-3;%
anal_opts.atom_laser.t0=1.71056; %1.6373;%center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=35;

anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
tmp_xlim=[-45,45];    
tmp_ylim=[-45, 45];
anal_opts.atom_laser.xylim=[tmp_xlim;tmp_ylim];
%anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import
anal_opts.atom_laser.plot.all=false;
anal_opts.atom_laser.fit.hist_sigma=3e-4;
anal_opts.atom_laser.global=anal_opts.global;


al_mcp_data=data.mcp_tdc.masked;
%al_mcp_data.counts_txy=al_mcp_data.counts_txy(1:20); % good for debuging
al_mcp_data.all_ok=al_mcp_data.num_counts>1e3;

tic
data.al_pulses=fit_gauss_to_al_pulses(anal_opts.atom_laser,al_mcp_data);
toc


% save('20190708_al_fits_done.mat','-v7.3')




%%
anal_opts.heating_fit=[];
data.mcp_tdc.all_ok=al_mcp_data.all_ok;

data.signal.heating = forbidden_signal_heating(data,anal_opts.heating_fit);
%% Generate signal
% data.signal.total_num = signal_process(data,anal_opts);


anal_opts.cal_model=[];
anal_opts.cal_model.smooth_time=60*20;
anal_opts.cal_model.plot=true;
% Create a calibration model
data.cal = make_cal_model_heating(anal_opts.cal_model,data);

%
signal_bined=[];
signal_unbinned.msr.val = data.cal.calibrated_signal.val;
signal_unbinned.msr.freq = data.signal.heating.msr.freq;
%probe_freq_bins=col_vec(linspace(-80,80,60));
probe_freq_bins=col_vec(linspace(-5,15,15));
iimax=numel(probe_freq_bins)-1;
for ii=1:iimax
    signal_bined.freq_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=signal_unbinned.msr.freq<probe_freq_bins(ii+1) & signal_unbinned.msr.freq>probe_freq_bins(ii);
    signal_bined.val(ii,:)=nanmean(signal_unbinned.msr.val(bin_mask,:),1);
    %sum(bin_mask)
    %std(signal_unbinned.val(bin_mask))
    signal_bined.unc_val(ii,:)=nanstd(signal_unbinned.msr.val(bin_mask,:),[],1)./sqrt(sum(bin_mask));
    
    signal_bined.freq(ii)=nanmean(probe_freq_bins(ii:ii+1));
    signal_bined.freq_lims_diff(ii,:)=abs(signal_bined.freq_lims(ii,:)-signal_bined.freq(ii));
end

stfig('Heating method calibrated')
clf

tot_plots=size(signal_unbinned.msr.val,2);
signal_idx=1;
for signal_idx=1:tot_plots
%     ymultipler=signal_unbinned.ymult(signal_idx);
%     ylabel_str=signal_unbinned.ystr(signal_idx);
    subplot(2,tot_plots,0+signal_idx)
    plot(signal_unbinned.msr.freq,signal_unbinned.msr.val(:,signal_idx).*1e9,'x')
    xlabel('freq-theory (MHz)')
    %ylabel('heating rate (nk/s)')
    ylabel('dif in heating rate (nK/s)')
%     title(signal_unbinned.names{signal_idx})
    xlim([min(probe_freq_bins),max(probe_freq_bins)])
    xl=xlim;
    subplot(2,tot_plots,tot_plots+signal_idx)
    errorbarxy(col_vec(signal_bined.freq),signal_bined.val(:,signal_idx).*1e9,...
        signal_bined.freq_lims_diff,signal_bined.unc_val(:,signal_idx).*1e9);
    xlim(xl)
    xlabel('freq-theory (MHz)')
    ylabel('dif in heating rate (nK/s)')
%     title(signal_unbinned.names{signal_idx})
end
saveas(gca,fullfile(anal_opts.global.out_dir,'signal_fits.png'))
%% Save to output
cli_header({0,'Saving output...'})
out_data.data.cal = data.cal;
out_data.data.signal = data.signal.heating;
out_data.data.al_pulses = data.al_pulses;
out_data.data.num = data.mcp_tdc.num_counts;
out_data.data.all_ok = data.mcp_tdc.all_ok;
out_data.data.masked_num = data.mcp_tdc.masked.num_counts;
out_data.data.ai_log = data.ai_log;
out_data.options = anal_opts;
save(fullfile(anal_opts.global.out_dir,'data_results.mat'),'out_data')
fwtext('All Done!')

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
