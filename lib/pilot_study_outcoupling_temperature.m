%here i want to see what kin od heating rate we can be sensitive too

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

hebec_constants 

%% variables
% add all subfolders to the path
% % Setting up
% BEGIN USER VAR-------------------------------------------------
anal_opts=[]; %reset the options (would be good to clear all variables except the loop config
anal_opts.tdc_import.dir='C:\Users\Bryce\Dropbox\UNI\project\programs\Forbidden_Transition\data\20190704_forbidden_heating_sample_data\700939267.0-1.7Mhz';
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

% END USER VAR-----------------------------------------------------------
%%


anal_opts.global.fall_velocity=const.g0*anal_opts.global.fall_time;

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


%%

%% BINNING UP THE ATOM LASER PULSES
%now find the mean position of each pulse of the atom laser in each shot

anal_opts.atom_laser.pulsedt=120e-3;
anal_opts.atom_laser.t0=1.637327; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=189;
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;

tmp_xlim=[-30e-3, 30e-3];    
tmp_ylim=[-30e-3, 30e-3];
anal_opts.atom_laser.xylim=[tmp_xlim;tmp_ylim];
%anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import
anal_opts.atom_laser.plot.all=false;

data.mcp_tdc.all_ok=data.mcp_tdc.num_counts> 0;

data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);

%%
thermal_fit=[];
plot_therm_fit=1; 
for shot_idx=1:numel(data.mcp_tdc.shot_num)

    temp_multipler_plot=1e-6;
    pulse_time=data.mcp_tdc.al_pulses.time_cen;
    pulse_time=pulse_time-pulse_time(1);
    temperature_from_time=(anal_opts.global.fall_velocity*data.mcp_tdc.al_pulses.pos.std(shot_idx,:,1)/anal_opts.global.fall_time).^2 *const.mhe/const.kb;
    % if we assume a normal distribution we can estimate the error in the estimated std
    % std(sigma)= \[Sigma] Sqrt[1 - c4^2]
    %std_in_std=data.mcp_tdc.al_pulses.pos.std(5,:,1).*sqrt(1-normal_correction_c4(data.mcp_tdc.al_pulses.num_counts(shot_number,:)).^2);
    %approx way
    frac_std_in_std=1./sqrt(2*(data.mcp_tdc.al_pulses.num_counts(shot_idx,:)-1));
    % the power makes the fractional uncert a factor of 2 more
    temperature_uncert=frac_std_in_std.*2.*temperature_from_time;


    xdata=pulse_time;
    ydata=temperature_from_time;
    yunc=temperature_uncert;

%     modelfun= @(b,x) b(1)+x.*b(2)+x.^2.*b(3);
%     beta0=[0.2,0.2,0];
%     cof_names={'offset','lin','square'};

    modelfun= @(b,x) b(1)+x.*b(2);
    beta0=[0.2,0];
    cof_names={'offset','lin'};

    %exp(-x(:,1).*max(0,b(7))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(8)*x(:,1)+b(5)*x(:,2)+b(6)*x(:,3);
    weights=1./(yunc.^2);
    weights(isnan(weights))=1e-20; %if nan then set to smallest value you can
    weights=weights/sum(weights);
     opt = statset('TolFun',1e-10,'TolX',1e-10,'MaxIter',1e4,...
                'UseParallel',1);

    fitobject=fitnlm(xdata,ydata,modelfun,beta0,...
        'Weights',weights,'CoefficientNames',cof_names,...
        'options',opt);
    thermal_fit.fit_obj{shot_idx}=fitobject;
    thermal_fit.RMSE(shot_idx)=fitobject.RMSE;
    thermal_fit.fit_coeff.SE(shot_idx,:)=fitobject.Coefficients.SE;
    thermal_fit.fit_coeff.Estimate(shot_idx,:)=fitobject.Coefficients.Estimate;
    
    if plot_therm_fit
        stfig('single shot temp. curve');
        clf
        errorbar(pulse_time,temperature_from_time/temp_multipler_plot,temperature_uncert/temp_multipler_plot,'ko','CapSize',0)
        xlabel('pulse time')
        ylabel('Temperature (uk)')

        mdl_samp_x=col_vec(linspace(min(xdata),max(xdata),1e3));
        [mdl_samp_y,mdl_samp_ci]=predict(fitobject,mdl_samp_x);

        hold on
        plot(mdl_samp_x,mdl_samp_y/temp_multipler_plot,'k-')
        plot(mdl_samp_x,mdl_samp_ci/temp_multipler_plot,'-','color',[1,1,1]*0.5)
        hold off
        pause(0.01)
    end
    
end

%%
[heating_mean,heating_unc]=unc_wmean(thermal_fit.fit_coeff.Estimate(:,2),thermal_fit.fit_coeff.SE(:,2));

heating_rate_str=string_value_with_unc(1e9*heating_mean,1e9*heating_unc,'b');
fprintf('%s\n',anal_opts.tdc_import.dir)
fprintf('mean heating rate %s nk/s \n',heating_rate_str)
%%

mean_widths=squeeze(mean(data.mcp_tdc.al_pulses.pos.std,1));

plot(mean_widths(:,1))
