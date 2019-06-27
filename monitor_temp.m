% monitor an experiment that uses measrument of the trap freq
% It would be usefull to get a decent idea of what the trap frequency is doing during a run without
% the need for a full processing of the data as in main
%  - processing each shot as it is made
%  - plot a history of the trap freq

% Known BUGS/ Possible Improvements
%
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-01
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190521_search_forbidden_6_trap_freq\';
data_dir = anal_opts.tdc_import.dir;
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=true;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=100;

anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.95;

anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


anal_opts.osc_fit.binsx=1000;
anal_opts.osc_fit.blur=1;
anal_opts.osc_fit.xlim=[-20,20]*1e-3;
anal_opts.osc_fit.tlim=[0.86,1.08];
anal_opts.osc_fit.dimesion=2; %Sel ect coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=50;

xmin = -25e-3;
xmax = 25e-3;
ymin = -25e-3;
ymax = 25e-3;
xybins = 500;
t_window_min = 0;
t_window_max = 2;

% END USER VAR-----------------------------------------------------------
fclose('all')
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

hebec_constants
anal_opts.tdc_import.mat_save=false;
anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;

if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end

%anal_out.dir=sprintf('%sout\\monitor\\',...
%    anal_opts.tdc_import.dir);

anal_out.dir=[fullfile(anal_opts.tdc_import.dir,'out','monitor'),filesep];
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;



%%
temp_history=[];
temp_history.temp_val=[];
temp_history.temp_unc=[];
temp_history.shot_num=[];
sfigure(1);
set(gcf,'color','w')
clf;
%%
handles.masshe=6.646400000000001e-27;
handles.boltzconst=1.380648520000000e-23;
handles.falltime=4.160000000000000e-01;
handles.trapfreqrad=500;
handles.trapfreqaxial=50;
handles.hbar=1.054571800000000e-34;
%%
% batch_data=[];
% batch_data.shot_num=[];
% anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
% batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
%% LabView import

opts.lv.plots = true;
%% Wavemeter log importing

opts.wm.force_reimport=false;
opts.wm.num_logs = nan;
opts.wm.plots = true;
opts.wm.plot_all=1;
opts.wm.plot_failed=false;


opts.wm.wm_log_name='log_wm_';

opts.wm.cache_import.verbose=0;
opts.wm.cache_import.force_recalc=0;

opts.wm.cache_import.save_compressed=true;%needed otherwise save takes a very long time
opts.wm.cache_import.path_directions={1,'dir'};

opts.wm.plot_all=true;
opts.wm.plot_failed=false;
opts.wm.force_reimport=false;
opts.dir = data_dir;
%% TDC import


opts.tdc.plots = true;
opts.tdc.file_name='d';
opts.tdc.force_load_save=false;   %takes precidence over force_reimport
opts.tdc.force_reimport=false;
opts.tdc.force_forc=false;
opts.tdc.dld_xy_rot=0.61;
%Should probably try optimizing these
tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
opts.tdc.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.max_runtime=inf;%inf%cut off the data run after some number of hours, should bin included as its own logic not applied to the atom number ok
%%
    if (exist([opts.dir,'out'], 'dir') == 0), mkdir([opts.dir,'out']); end
    %make a subfolder with the ISO timestamp for that date
    opts.out_dir=sprintf('%sout\\%s\\',...
        opts.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
    if (exist(opts.out_dir, 'dir') == 0), mkdir(opts.out_dir); end
    opts.lv.dir = opts.dir;
    opts.lv.out_dir = opts.out_dir;
    opts.ai.dir=opts.dir;
    opts.ai.out_dir = opts.out_dir;
    opts.wm.dir=opts.dir;
    opts.wm.out_dir = opts.out_dir;
    opts.tdc.dir = opts.dir;
    opts.tdc.out_dir = opts.out_dir;
    opts.ai.cache_single.mock_working_dir=opts.dir;
    opts.wm.cache_import.mock_working_dir=opts.dir;
    wm_logs=dir([opts.wm.dir,opts.wm.wm_log_name,'*.txt']);
    opts.wm.names={wm_logs.name};
opts.check.ai_override = 1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Import the analog log files
% data.ai = ai_log_import(opts);
% %% Import LabView log
data.lv = import_lv_log(opts);
% % %% Import wavemeter logs
data.wm = wm_log_import(opts);
% % %% Import TDC files
data.tdc = import_mcp_tdc(opts);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Match the timestamps    
data.sync = match_timestamps(data,opts);

%% Create a calibration model
data.cal = make_calibration_model(data,opts);

%% Mask out shots which failed
data.check = check_for_errors(data,opts);

%% Break data into categories
data.cat = categorize_shots(data,opts);

%% Grouping by wavelength 
%data = bin_by_wavelength(data,opts);
%%
temp_history.temp_val = [];
for c_data = data.tdc.counts_txy
    tdc_data = cell2mat(c_data);
mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
    & tdc_data(:,2)>xmin & tdc_data(:,2)<xmax ...
    & tdc_data(:,3)>ymin & tdc_data(:,3)<ymax ;
data_windowed = tdc_data(mask,:);
X_bin=linspace(xmin,xmax,xybins);
X_bin_width=(xmax-xmin)/xybins;
[X1d_counts,X1d_edges]=histcounts(data_windowed(:,2),X_bin);
X1d_centers=mean([X1d_edges(1:end-1);X1d_edges(2:end)]);
X1d_counts=1e-3*X1d_counts/(X_bin_width);%why 1e-3 factor?
sfigure(5);
clf
set(gcf,'color','w')
plot(X1d_centers*1000,X1d_counts/X_bin_width,'k')
xlabel('x (mm)')
ylabel('count rate')
fit_params = fit_cond(handles,X1d_centers*1000,X1d_counts/X_bin_width);
val = fit_params;
temp_history.temp_val = [temp_history.temp_val;val];
%pause(0.5)
end
temp_history.shot_num = data.tdc.shot_num;

%%
sfigure(1);
% errorbar(trap_freq_history.shot_num,...
%     trap_freq_history.trap_freq_val,trap_freq_history.trap_freq_unc,...
%     'kx-','MarkerSize',7,'CapSize',0,'LineWidth',1.5)
scatter(temp_history.shot_num,temp_history.temp_val(:,3),7,'kx')
grid on
h=gca;
grid on    % turn on major grid lines
grid minor % turn on minor grid lines
% Set limits and grid spacing separately for the two directions:
% Must set major grid line properties for both directions simultaneously:
h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
h.GridAlpha=1;  % the default is partially transparent
h.GridColor=[0,0,0]; % here's the color for the major grid lines
% Idem for minor grid line properties:
h.MinorGridLineStyle='-';
h.MinorGridAlpha=0.1;
h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
xlabel('Shot Number')
ylabel('Signal')
pause(1e-6)
%saveas(gcf,fullfile(anal_out.dir,'freq_history.png')
%%
function fit_params=fit_therm(xdata,ydata)
amp_guess=max(ydata);
mu_guess=sum(xdata.*ydata)/sum(ydata); %compute the weighted mean
sig_guess=sum((xdata-mu_guess).^2.*ydata)/sum(ydata); %compute the mean square weighted deviation
fo = statset('TolFun',10^-6,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);
fitobject=fitnlm(xdata,ydata,...
    'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))',...
    [amp_guess,mu_guess,sig_guess],...
    'CoefficientNames',{'amp','mu','sig'},'Options',fo);
fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
end


function out=fit_cond(handles,xdata,ydata)
%idealy this would be done with a constrained fit but ticky to implement in
%matlab
%roughly following https://arxiv.org/pdf/0811.4034.pdf

%dim = [0.2 0.6 0.3 0.3];
%str = {'Straight Line Plot','from 1 to 10'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');

% thermfitparms=fit_therm(xdata,ydata); %{'amp','mu','sig'}
% 
% amp_guess=thermfitparms(1,1);
% mu_guess=thermfitparms(2,1);
% sig_guess=thermfitparms(3,1);

amp_guess=max(ydata);
mu_guess=sum(xdata.*ydata)/sum(ydata); %compute the weighted mean
sig_guess=sum((xdata-mu_guess).^2.*ydata)/sum(ydata); %compute the mean square weighted deviation

fo = statset('TolFun',10^-8,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);
%params in order
%1 parabola height (unnorm)
%2 center
%3 TF rad
%4 gauss peak height
%5 gauss width
modelfun=@bimod;
fitobject=fitnlm(xdata,ydata,modelfun,...
    [amp_guess*2/3,mu_guess,sig_guess/10,amp_guess*1/3,sig_guess*1.1],'Options',fo);


xvalues=linspace(min(xdata),max(xdata),300);
fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
hold on
plot(xvalues,justcond(fit_params(:,1),xvalues),'r--')
plot(xvalues,justtherm(fit_params(:,1),xvalues),'r')
plot(xvalues,bimod(fit_params(:,1),xvalues),'b')
hold off


%to convert this to temp we use the expression 2.41 from pethick
%momentum width=(mkT)^{1/2}
%T=sqrt(pwidth)/mk
%to convert to what we have
%pwidth=xwidth*m/tfall
%for the time axis we must convert to spatail using velocity
%xwidth=twidth*velocity=vwidth*1/2 g t^2
%total expression then T=(twidth^2*velocity^2)/(k*tfall^2*tfall^2)
%handles.falldist = .848; %fall distance of 848mm
%handles.falltime = .416; %fall time of 416ms

    vdet=1/1000; %to cancel the plot being in mm
    units='mm';

%1 parabola height (unnorm)
%2 center
%3 TF rad
%4 gauss peak height
%5 gauss width

fit_params(3,1)=fit_params(3,1)*vdet;
fit_params(3,2)=fit_params(3,2)*vdet;
fit_params(5,1)=fit_params(5,1)*vdet;
fit_params(5,2)=fit_params(5,2)*vdet;


%here i basicaly integrate under the curve and find the ratio
%to give Nc/Ntot beause i have changed the scale of the plots to put things
%into nice units the absolout value of the counts is a bit odd
ampTF=fit_params(1,1);
ampgauss=fit_params(4,1);
countsTF=abs((1/(4*fit_params(3,1)))*ampTF);
countsGauss=abs((1/(sqrt(2*pi)*fit_params(5,1)))*ampgauss);
Ncondfrac=countsTF/(countsTF+countsGauss);
TonTc=(1-Ncondfrac)^(1/3);

T=(abs(fit_params(5,1)))^2 *handles.masshe /(handles.boltzconst*handles.falltime^2);

%should use 2.90 from pethick to correct for our cigar trap
Tc=T/TonTc;

omegabar=(2*pi*2*pi*2*pi*handles.trapfreqrad*handles.trapfreqrad*handles.trapfreqaxial)^(1/3);

Nest=(handles.boltzconst*Tc/(0.94*handles.hbar*omegabar))^3;

str=sprintf('GaussFit Radius %0.2e ± %0.1e %s \n Temp.(no interactions)%0.2e k \n TF radius %0.2e ± %0.1e %s \n Condensate fraction %0.1f%%\nT/Tc %0.1f%%\n Tc %0.2ek\n Est. N %0.2e',...
    abs(fit_params(3,1)),fit_params(3,2),units,T,abs(fit_params(5,1)),abs(fit_params(5,2)),units,Ncondfrac*100,TonTc*100,Tc,Nest);
text(0.02,0.9,str,'Units','normalized','VerticalAlignment','top');

out = [abs(fit_params(3,1)),fit_params(3,2),T,abs(fit_params(5,1)),abs(fit_params(5,2)),Ncondfrac,TonTc,Tc,Nest,sig_guess];

end


function out=bimod(b,x)
%1 parabola height (unnorm)
%2 center
%3 TF rad
%4 gauss peak height
%5 gauss width
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos
parabola=(1-((x(:)-b(2))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
therm=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
out=real(b(1).*max(zerosformax,parabola)+therm);
end

function out=justcond(b,x)
  %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width   
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos

% if b(5)<b(3)*1.5
%     b(1)=0; %this is a very hacky way to set a limit
% end

parabola=(1-((x(:)-b(2))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
out=real(b(1).*max(zerosformax,parabola));
end

function out=justtherm(b,x)
  %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width   
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos

out=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
end