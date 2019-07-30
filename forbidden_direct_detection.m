%% Data analysis for Helium forbidden transition spectroscopy
%% To do
% - Remove correlations to atom number []
% - reintroduce ai import [x]
%   - add total integrated power calibration [x]
%   - add ai checks (mainly when laser is multi mode) []
% - Create functionanilty to combine multiple runs together (wm drift model
% mainly) []
% - ensure masking is optimal
% - general organisation and commenting []
%%
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
% anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190713_forbidden427_direct_det_narrow\';
anal_opts.tdc_import.dir='Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190711_forbidden427_direct_det_bad_atom_num\'
anal_opts.tdc_import.save_cache_in_data_dir=true;
tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
tlim=[0,inf];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];
anal_opts.dld_aquire=22;
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
% fix multimode detection
% 
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

anal_opts.wm_log.red_sd_thresh=5; %allowable standard deviation in MHz
anal_opts.wm_log.red_range_thresh=10; %allowable range deviation in MHz
anal_opts.wm_log.rvb_thresh=20; %allowable value of abs(2*red-blue)

anal_opts.wm_log.global=anal_opts.global;

data.wm_log.proc=wm_log_process(anal_opts,data);
clear('sub_data')



%% Match the timestamps    
%% data.sync = match_timestamps(data,import_opts);
%% bryce change, this is now in the wavemeter proceesing code
fprintf('saving status...')
%save('20190704_data_imported.mat','-v7.3')
fprintf('done\n')

%% masking out hot spots

tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
%tlim=[0.5,6.2];
tlim=[1.5,22.5];
tlim_tot=[22.6,25];%0 to 25 or 22.6 to 25?

anal_opts.hotspot_mask.square_mask=[tlim;tmp_xlim;tmp_ylim];
anal_opts.hotspot_mask.square_mask_tot=[tlim_tot;tmp_xlim;tmp_ylim];
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
                            [-10.6e-3,15.56e-3,1.5e-3,0];
                              ];
do_mask=true;

if do_mask
    num_shots=numel(data.mcp_tdc.counts_txy);
    empty_shots=cellfun(@isempty,data.mcp_tdc.counts_txy);
    data.mcp_tdc.masked.counts_txy={};
    data.mcp_tdc.masked.num_counts=data.mcp_tdc.num_counts*nan;
    
    data.mcp_tdc.masked.tot_counts_txy={};
    data.mcp_tdc.masked.tot_num_counts=data.mcp_tdc.num_counts*nan;
    
    fprintf('masking shots %04u:%04u',num_shots,0)
    for ii=1:num_shots
        txy_shot=data.mcp_tdc.counts_txy{ii};
        if ~isempty(txy_shot)
            txy_shot_sig=masktxy_square(txy_shot,anal_opts.hotspot_mask.square_mask);
            txy_shot_sig=masktxy_2d_circle(txy_shot_sig,anal_opts.hotspot_mask.circ_mask);
            data.mcp_tdc.masked.num_counts(ii)=numel(txy_shot_sig);
            data.mcp_tdc.masked.counts_txy{ii}=txy_shot_sig;
            
            txy_shot_tot=masktxy_square(txy_shot,anal_opts.hotspot_mask.square_mask_tot);
            txy_shot_tot=masktxy_2d_circle(txy_shot_tot,anal_opts.hotspot_mask.circ_mask);
            data.mcp_tdc.masked.tot_num_counts(ii)=numel(txy_shot_tot);
            data.mcp_tdc.masked.tot_counts_txy{ii}=txy_shot_tot;
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
%% warning if you do this with all the counts will likely run out of memory
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

%% creation of the signal

data.mcp_tdc.ok.all=col_vec(data.mcp_tdc.num_counts)>2.5e3;
data.signal.raw=[];
offset = 0;
data.signal.raw.val=col_vec((data.mcp_tdc.masked.num_counts-offset)./data.mcp_tdc.num_counts);% %./data.mcp_tdc.num_counts
% data.signal.raw.val=col_vec(-data.mcp_tdc.masked.num_counts+data.mcp_tdc.num_counts.*0.1602-6.786e3);
% data.signal.raw.val=col_vec(data.mcp_tdc.num_counts.*0.1602-6.786e3);
%data.signal.raw.val=col_vec(data.mcp_tdc.masked.num_counts./data.mcp_tdc.masked.tot_num_counts);% %./data.mcp_tdc.num_counts
%-1/6.3508.*(data.mcp_tdc.masked.tot_num_counts-1.190409565626061e+05)
% data.signal.raw.unc=sqrt(data.signal.raw.val);
data.signal.raw.unc=ones(size(data.signal.raw.val)).*1e-10;

%estimate total incident energy from integrated pd voltage
data.incident_power = data.ai_log.pd.integrated.*3.2616e-5;

anal_opts.cal_model=[];
anal_opts.cal_model.smooth_time=300;
anal_opts.cal_model.plot=true;

out_calibration=make_cal_model_signal(anal_opts.cal_model,data);
data.signal.cal=out_calibration;


% bin data points and fit

predicted_freq=700939267; %MHz


%set up the colors to use
colors_main=[[233,87,0];[33,188,44];[0,165,166]];
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=14;

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);




is_cal=data.mcp_tdc.probe.calibration;
is_cal(isnan(is_cal))=1;

if ~use_ai
    is_shot_good=data.wm_log.proc.probe.freq.act.std<5 &...
        (data.wm_log.proc.probe.freq.set-data.wm_log.proc.probe.freq.act.mean)<10 &...
        ~is_cal;
else
    is_shot_good=data.wm_log.proc.probe.freq.act.std<5 &...
        (data.wm_log.proc.probe.freq.set-data.wm_log.proc.probe.freq.act.mean)<10 &...
        ~is_cal & ~data.ai_log.ok.sfp;
end

probe_freq=data.wm_log.proc.probe.freq.act.mean*2;%freq in blue
probe_freq=probe_freq-predicted_freq;

signal_unbinned.val=cat(2,data.signal.cal.calibrated_signal.val(is_shot_good));
signal_unbinned.freq=probe_freq(is_shot_good);

sigma_from_mean=(signal_unbinned.val-nanmean(signal_unbinned.val,1))/nanstd(signal_unbinned.val,1);
is_not_oulier=sigma_from_mean<7;
signal_unbinned.val=signal_unbinned.val(is_not_oulier);
signal_unbinned.freq=signal_unbinned.freq(is_not_oulier);

signal_unbinned.names={'RF knife method'};
 signal_unbinned.ystr={'signal (J^{-1})'};
signal_unbinned.ymult=[1];


fprintf('unbinned mean of counts               %f \n',nanmean(signal_unbinned.val(:,1)))
fprintf('unbinned standard deviation of counts %f \n',nanstd(signal_unbinned.val(:,1)))


signal_bined=[];
[grouped_values,group_pop,group_idx]=uniquetol(signal_unbinned.freq,0.03);
min_group_pop=2;
group_mask=group_pop>min_group_pop;
group_mask_index=col_vec(1:numel(grouped_values));
group_mask_index=group_mask_index(group_mask);
operation_on_group=@(y) arrayfun(@(x) y(signal_unbinned.freq(group_idx==x)),group_mask_index);


suggested_limits=[min(grouped_values),max(grouped_values)]
suggested_bins=numel(grouped_values)+1;
suggested_bin_width=mean(diff(grouped_values));
bin_width=2;

%probe_freq_bins=col_vec(linspace(-70,70,42));
%probe_freq_bins=col_vec(linspace(min(suggested_limits)-bin_width/2,max(suggested_limits)+bin_width/2,suggested_bins));
%probe_freq_bins=col_vec(linspace(-33.4-bin_width/2,20.7+bin_width/2,42));
%probe_freq_bins=col_vec(linspace(-40-bin_width/2,40+bin_width/2,81));
% probe_freq_bins=[-10,10,100]
%probe_freq_bins=col_vec(linspace(-25,20,9));

probe_freq_bins=cat(1,operation_on_group(@min)-100*eps,inf);

iimax=numel(probe_freq_bins)-1;
signal_bined.freq_std=nan(iimax,1);
signal_bined.val=nan(iimax,1);
signal_bined.unc_val=nan(iimax,1);
signal_bined.freq_mean=nan(iimax,1);
signal_bined.freq_obs_min_max_mean_diff=nan(iimax,2);
for ii=1:iimax
    signal_bined.freq_bin_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=signal_unbinned.freq<=probe_freq_bins(ii+1) & signal_unbinned.freq>probe_freq_bins(ii);
     signal_bined.freq_bin_cen(ii)=nanmean(probe_freq_bins(ii:ii+1));
    if sum(bin_mask)==0
        warning('no elements')
        signal_bined.num_bin(ii)=0;
    else
        signal_bined.num_bin(ii)=sum(bin_mask);
        signal_bined.val(ii,:)=nanmean(signal_unbinned.val(bin_mask,:),1);
        %sum(bin_mask)
        %std(signal_unbinned.val(bin_mask))
        signal_bined.unc_val(ii,:)=nanstd(signal_unbinned.val(bin_mask,:),[],1)./sqrt(sum(bin_mask));
        signal_bined.freq_mean(ii)=nanmean(signal_unbinned.freq(bin_mask));
        signal_bined.freq_std(ii)=nanstd(signal_unbinned.freq(bin_mask));
        signal_bined.freq_obs_min_max(ii,:)=[min(signal_unbinned.freq(bin_mask)),max(signal_unbinned.freq(bin_mask))];
        signal_bined.freq_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_bin_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_obs_min_max_mean_diff(ii,:)=abs(signal_bined.freq_obs_min_max(ii,:)-signal_bined.freq_mean(ii));
     end
end

stfig('counts vs probe freq normalised');
clf

gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2)); 
coeff_names={'amp','mu','sig'};

tot_plots=size(signal_unbinned.val,2);
signal_idx=1;
for signal_idx=1:tot_plots
   
    ymultipler=signal_unbinned.ymult(signal_idx);
    ylabel_str=signal_unbinned.ystr(signal_idx);
    
    subplot(2,tot_plots,0+signal_idx)
    plot(signal_unbinned.freq,signal_unbinned.val(:,signal_idx)*ymultipler,'x')
    hold on
    plot(probe_freq_bins,probe_freq_bins*0,'s')
    hold off
    xlabel('freq-theory (MHz)')
    %ylabel('heating rate (nk/s)')
    ylabel(ylabel_str)
    title(signal_unbinned.names{signal_idx})
    xlim([min(probe_freq_bins),max(probe_freq_bins)])
    xl=xlim;
    
    xdata=col_vec(signal_bined.freq_mean);
    ydata=signal_bined.val(:,signal_idx);
    yunc=signal_bined.unc_val(:,signal_idx);
    %mean(ydata-predict(fitobject,xdata))
    subplot(2,tot_plots,tot_plots+signal_idx)
    %plot(col_vec(signal_bined.freq_bin_cen),ydata*ymultipler,'x')
    plot(xdata,ydata*ymultipler,'o','MarkerSize',5,'MarkerFaceColor',colors_detail(1,:))
    errorbar(xdata,ydata*ymultipler,...
        signal_bined.unc_val(:,signal_idx)*ymultipler,signal_bined.unc_val(:,signal_idx)*ymultipler,...
         signal_bined.freq_obs_min_max_mean_diff(:,1), signal_bined.freq_obs_min_max_mean_diff(:,2),...
        'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
         'MarkerFaceColor',colors_detail(3,:),'LineWidth',2.5);
    
    %errorbarxy(xdata,ydata*ymultipler,...
    %    signal_bined.freq_obs_min_max_mean_diff,signal_bined.unc_val(:,signal_idx)*ymultipler);
    xlim(xl)
    xlabel('freq-theory (MHz)')
    ylabel(ylabel_str)
    title(signal_unbinned.names{signal_idx})
    
    % fit a gaussian
   
    is_data_not_nan=~isnan(xdata) & ~isnan(ydata) & ~isnan(yunc);
    xdata=xdata(is_data_not_nan);
    ydata=ydata(is_data_not_nan);
    yunc=yunc(is_data_not_nan);
    
    amp_guess=max(ydata);
    ydata_shifted =ydata-min(ydata);
    mu_guess=-7.8;%wmean(xdata,ydata_shifted); %compute the weighted mean
    %sig_guess=sqrt(nansum((xdata-mu_guess).^2.*ydata_shifted)/nansum(ydata_shifted)); %compute the mean square weighted deviation
    sig_guess=5;
    fo = statset('TolFun',10^-6,...
        'TolX',1e-4,...
        'MaxIter',1e4,...
        'UseParallel',1);
    % 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
    inital_guess=[amp_guess,mu_guess,sig_guess];
    fitobject=fitnlm(xdata,ydata,...
        gauss_fun1d,...
         inital_guess,...
        'CoefficientNames',coeff_names,'Options',fo);
    fit_coeff=fitobject.Coefficients.Estimate;
    fit_se=fitobject.Coefficients.SE;
    x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
    [ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    hold on
    plot(x_sample_fit,ysamp_val*ymultipler,'r')
    drawnow
    yl=ylim;
    plot(x_sample_fit,ysamp_ci*ymultipler,'color',[1,1,1].*0.5)
    ylim(yl)
    % show the inital guess
    %plot(x_sample_fit,gauss_fun1d(inital_guess,x_sample_fit)*ymultipler)
    fitobject
    
    amp_str=string_value_with_unc(fitobject.Coefficients.Estimate(1),fitobject.Coefficients.SE(1),'b');
    cen_str=string_value_with_unc(fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2),'b');
    width_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject.Coefficients.SE(3),'b');
    offset_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject.Coefficients.SE(3),'b');
    width_units='MHz';
    offset_units='counts';
    amp_units='counts';
    str=sprintf('Gauss fit \n   Cen    %s %s \n   Width %s %s \n   Amp   %s %s \n   Offset %s %s',...
        cen_str,width_units,width_str,width_units,amp_str,amp_units,offset_str,offset_units);
    text(0.01,0.9,str,'Units','normalized'); 
    
    fprintf('transition frequnency %s\n',string_value_with_unc(predicted_freq+fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2)))
                    
    
end

saveas(gca,fullfile(anal_opts.global.out_dir,'signal_fits.png'))
clear out_data
out_data.signal = data.signal;
out_data.freq = probe_freq;
out_data.atom_num = data.mcp_tdc.masked.tot_num_counts;
out_data.integrated_pd = data.ai_log.pd.integrated;
out_data.probe_num = data.mcp_tdc.masked.num_counts;
out_data.time = data.mcp_tdc.time_create_write(:,2);
out_data.is_shot_good = is_shot_good;

out_data.opts = anal_opts;

save(fullfile(anal_opts.global.out_dir,'data_results.mat'),'out_data')

%% residuals
freq_window = [3.5,7];
freq_window = [-60,60];
        
not_empty_shots=~cellfun(@isempty,data.mcp_tdc.masked.counts_txy);
freq_mask=signal_unbinned.freq<=freq_window(2) & signal_unbinned.freq>freq_window(1);

int_pd = data.ai_log.pd.integrated(is_shot_good);
int_pd = int_pd(is_not_oulier);
int_num = data.mcp_tdc.masked.num_counts(is_shot_good);
int_num = int_num(is_not_oulier);
atom_num = data.mcp_tdc.masked.tot_num_counts(is_shot_good);
atom_num = atom_num(is_not_oulier);
ypred_val=predict(fitobject,signal_unbinned.freq,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
res = signal_unbinned.val - ypred_val;
sfigure(345);
clf
plot(res)
xlabel('residuals')
ylabel('shot indx')
sfigure(456);
clf
corr_plot(int_num,res,ones(size(res)))
ylabel('res')
xlabel('background counts')
sfigure(567);
clf
corr_plot(atom_num,res,ones(size(res)))
ylabel('res')
xlabel('total atom counts')
sfigure(678);
clf
corr_plot(signal_unbinned.freq,res,ones(size(res)))
ylabel('res')
xlabel('freq - theory (MHz)')
sfigure(123)
clf
corr_plot(int_pd(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('average pd voltage')
sfigure(4321)
clf
res_mdl =  corr_plot(int_num(freq_mask)./atom_num(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('atom num ratio')
%%
if true
stfig('residual correlations removed');
clf
xdata = signal_unbinned.freq;
ydata = signal_unbinned.val - predict(res_mdl,int_num'./atom_num');
% amp_guess=max(ydata);
% ydata_shifted =ydata-min(ydata);
% mu_guess=wmean(xdata,ydata_shifted); %compute the weighted mean
% %sig_guess=sqrt(nansum((xdata-mu_guess).^2.*ydata_shifted)/nansum(ydata_shifted)); %compute the mean square weighted deviation
% sig_guess=10;
% fo = statset('TolFun',10^-6,...
%     'TolX',1e-4,...
%     'MaxIter',1e4,...
%     'UseParallel',1);
% % 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
% inital_guess=[amp_guess,mu_guess,sig_guess,0];
% fitobject_adj=fitnlm(xdata,ydata,...
%     gauss_fun1d,...
%      inital_guess,...
%     'CoefficientNames',coeff_names,'Options',fo);
% fit_coeff=fitobject_adj.Coefficients.Estimate;
% fit_se=fitobject_adj.Coefficients.SE;
x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
[ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on
plot(x_sample_fit,ysamp_val*ymultipler,'r')
drawnow
yl=ylim;
plot(x_sample_fit,ysamp_ci*ymultipler,'color',[1,1,1].*0.5)
ylim(yl)
amp_str=string_value_with_unc(fitobject.Coefficients.Estimate(1),fitobject_adj.Coefficients.SE(1),'b');
cen_str=string_value_with_unc(fitobject.Coefficients.Estimate(2),fitobject_adj.Coefficients.SE(2),'b');
width_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject_adj.Coefficients.SE(3),'b');
offset_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject_adj.Coefficients.SE(3),'b');
width_units='MHz';
offset_units='counts';
amp_units='counts';
str=sprintf('Gauss fit \n   Cen    %s %s \n   Width %s %s \n   Amp   %s %s \n   Offset %s %s',...
    cen_str,width_units,width_str,width_units,amp_str,amp_units,offset_str,offset_units);
text(0.01,0.9,str,'Units','normalized');

% xtemp=col_vec(signal_bined.freq_mean);
% ytemp=signal_bined.val(:,signal_idx);
% plot(xtemp,ytemp*ymultipler,'x','MarkerSize',5)
hold on
for ii=1:iimax
    signal_bined.freq_bin_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=signal_unbinned.freq<=probe_freq_bins(ii+1) & signal_unbinned.freq>probe_freq_bins(ii);
     signal_bined.freq_bin_cen(ii)=nanmean(probe_freq_bins(ii:ii+1));
    if sum(bin_mask)==0
        warning('no elements')
        signal_bined.num_bin(ii)=0;
    else
        signal_bined.num_bin(ii)=sum(bin_mask);
        signal_bined.val(ii,:)=nanmean(ydata(bin_mask,:),1);
        %sum(bin_mask)
        %std(signal_unbinned.val(bin_mask))
        signal_bined.unc_val(ii,:)=nanstd(ydata(bin_mask,:),[],1)./sqrt(sum(bin_mask));
        signal_bined.freq_mean(ii)=nanmean(signal_unbinned.freq(bin_mask));
        signal_bined.freq_std(ii)=nanstd(signal_unbinned.freq(bin_mask));
        signal_bined.freq_obs_min_max(ii,:)=[min(signal_unbinned.freq(bin_mask)),max(signal_unbinned.freq(bin_mask))];
        signal_bined.freq_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_bin_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_obs_min_max_mean_diff(ii,:)=abs(signal_bined.freq_obs_min_max(ii,:)-signal_bined.freq_mean(ii));
     end
end

xdata=col_vec(signal_bined.freq_mean);
ydata=signal_bined.val(:,signal_idx);
yunc=signal_bined.unc_val(:,signal_idx);


%plot(col_vec(signal_bined.freq_bin_cen),ydata*ymultipler,'x')
plot(xdata,ydata*ymultipler,'o','MarkerSize',5,'MarkerFaceColor',colors_detail(1,:))
% errorbar(xdata,ydata*ymultipler,...
%     signal_bined.unc_val(:,signal_idx)*ymultipler,signal_bined.unc_val(:,signal_idx)*ymultipler,...
%     signal_bined.freq_obs_min_max_mean_diff(:,1), signal_bined.freq_obs_min_max_mean_diff(:,2),...
%     'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
%     'MarkerFaceColor',colors_detail(3,:),'LineWidth',2.5);
%ylim([min(ydata),max(ydata)])
%xlim([min(xdata),max(xdata)])
end
%% structure in the noise
anal_opts.plot2d.do=false;
if anal_opts.plot2d.do
    anal_opts.signal=[];
    anal_opts.plot2d.lim.x=[-45,45]*1e-3;
    anal_opts.plot2d.lim.y=[-45,45]*1e-3;
    anal_opts.plot2d.nbins=4e2;
    anal_opts.plot2d.blur=3;
    anal_opts.plot2d.cmp_dyn_range=false;
    freq_window = [4,6];
    
    
    
    not_empty_shots=~cellfun(@isempty,data.mcp_tdc.masked.counts_txy);
    bin_mask=signal_unbinned.freq<=freq_window(2) & signal_unbinned.freq>freq_window(1);
    
    temp_txy=data.mcp_tdc.masked.counts_txy(is_shot_good);
    temp_txy=temp_txy(is_not_oulier);
    temp_txy=temp_txy(bin_mask);
    num_shots_probe = size(temp_txy,1);
    probe_txy=cat(1,temp_txy{:});
    
    cal_txy=cat(1,data.mcp_tdc.masked.counts_txy{logical(is_cal)});
    num_shots_cal = size(cal_txy,1);
    
    all_txy = probe_txy;
    dyn_range_pow=0.2;
    num_bins=anal_opts.plot2d.nbins;
    edge_x=linspace(min(anal_opts.plot2d.lim.x),max(anal_opts.plot2d.lim.x),num_bins);
    edge_y=linspace(min(anal_opts.plot2d.lim.y),max(anal_opts.plot2d.lim.y),num_bins);
    
    spatial_blur=anal_opts.plot2d.blur;
    bin_area=(range(anal_opts.plot2d.lim.x)/num_bins)*(range(anal_opts.plot2d.lim.x)/num_bins);
    [counts,centers]=hist3(all_txy(:,2:3),'edges',{edge_x,edge_y});
    counts=counts/bin_area;
    counts=counts/num_shots_probe;
    
    [cal_counts,cal_centers]=hist3(cal_txy(:,2:3),'edges',{edge_x,edge_y});
    cal_counts=cal_counts/bin_area;
    cal_counts=cal_counts/num_shots_cal;
    
    %imagesc seems to plot the wrong way round so we transpose here
    
    if anal_opts.plot2d.cmp_dyn_range
        counts=counts.^dyn_range_pow;
    end
    if  ~spatial_blur==0
        counts=imgaussfilt(counts,spatial_blur);
        cal_counts=imgaussfilt(counts,spatial_blur);
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
    if anal_opts.plot2d.cmp_dyn_range
        xlabel(h,sprintf('Count Density^{%.2f} (m^{-2})',dyn_range_pow))
    else
        xlabel(h,'Count Density (m^{-2})')
    end
    stfig('counts during cal')
    imagesc(10^3*centers{1},10^3*centers{2},transpose(cal_counts))
    colormap(viridis())
    set(gca,'Ydir','normal')
    set(gcf,'Color',[1 1 1]);
    title('Spatial Dist. TOP')
    xlabel('X(mm)')
    ylabel('Y(mm)')
    h=colorbar;
    if anal_opts.plot2d.cmp_dyn_range
        xlabel(h,sprintf('Count Density^{%.2f} (m^{-2})',dyn_range_pow))
    else
        xlabel(h,'Count Density (m^{-2})')
    end
    stfig('dif in counts')
    imagesc(10^3*centers{1},10^3*centers{2},transpose(counts-cal_counts))
    colormap(viridis())
    set(gca,'Ydir','normal')
    set(gcf,'Color',[1 1 1]);
    title('Spatial Dist. TOP')
    xlabel('X(mm)')
    ylabel('Y(mm)')
    h=colorbar;
    if anal_opts.plot2d.cmp_dyn_range
        xlabel(h,sprintf('Count Density^{%.2f} (m^{-2})',dyn_range_pow))
    else
        xlabel(h,'Count Density (m^{-2})')
    end
end
%%
% import_opts.signal=[];
% import_opts.signal.plot.lim.x=[-45,45]*1e-3;
% import_opts.signal.plot.lim.y=[-45,45]*1e-3;
% import_opts.signal.plot.nbins=1e3;
% import_opts.signal.plot.blur=3;
% import_opts.signal.plot.cmp_dyn_range=true;
% tmp_xlim=[-50e-3, 50e-3];    
% tmp_ylim=[-50e-3, 50e-3];
% %tlim=[0.5,6.2];
% tlim=[5,22.5];
% import_opts.signal.square_mask=[tlim;tmp_xlim;tmp_ylim];
% import_opts.signal.circ_mask=[[0,0,35e-3,1];
%                               [35e-3,5e-3,7e-3,0];
%                               [25.05e-3,-19e-3,8e-3,0];
%                               [26.94e-3,21.98e-3,5e-3,0];
%                               [19.1e-3,-28.0e-3,4e-3,0];
%                               [2.973e-3,-33.15e-3,4e-3,0];
%                               [6.216e-3,34.41e-3,3e-3,0];
%                               [13.78e-3,31.62e-3,3e-3,0];
%                               [21.26e-3,25.95e-3,3e-3,0];
%                               [30.36e-3,-12.52e-3,5e-3,0];
%                               ];
% 
% data.signal.masked_number = forbidden_signal_masked_num(data,import_opts.signal);
% 
% 
% 
% 
% 
% %% Generate signal
% data.signal.total_num = signal_process(data,import_opts);
% 
% 
% 
% % Create a calibration model
% data.cal = make_calibration_model(data,import_opts);
% 
% %% Mask out shots which failed
% % data.check = check_for_errors(data,opts);
% 
% % %% Break data into categories
% % data.cat = categorize_shots(data,opts);
%  
% % %% Peak detection
% %data = auto_peak_detect(data,opts);
% 
% % %% Fit the detected peaks
% %data = fit_detected_peaks(data,opts);
% % data = fancy_fits(data,opts);
% % 
% % %% Zeeman shift correction
% %data = zeeman_correction(data,opts);
% 
% %% Presentation plots
% %data = present_plots(data,opts);

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