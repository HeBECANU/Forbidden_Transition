% anal_opts.tdc_import.dir=;
clear all




%might need an extra wavemeter point for this guy




%

%

%
%
%Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190710_forbidden427_direct_det\
%Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190710_forbidden427_direct_det_narrow_dither_on\
%Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190713_forbidden427_direct_det_narrow\
%Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190716_forbidden_Rf_2.00\
%Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190716_forbidden_Rf_2.05\
%Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190721_weekend_run\

%Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190715_forbidden427_narrow_changed_pol\
%Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190704_forbidden_long interrogation\



%combine the data from each director (all we need is the signal corrected
%for wvm drift)

%also need to acount for background correlation and pd variations as well

% wm drift model
%format unix time, offset, unc
wm_offset = [1.563278643754e9,-0.71111273765563965,0.001339
1.563183199886e9,-1.1711589097976685,  nan
1.563115846138e9,-1.7073591947555542, nan
1.563016150706e9,-0.57892698049545288, nan
%1.562933976169e9,-2.6968891024589539, nan
1.562861326011e9,-3.2293498516082764, nan
1.562724583967e9,-1.4483662247657776, nan
%1.562724143194e9,-1.1263280510902405, nan
%1.562212485100000e9,-1.1263280510902405, nan %this a bit dodge
%1.563749169809e9,1.8316859006881714,  nan 
1.5637490316e9,-1.85798978805542, nan
%1.563748921638e9,0.541365385055542, nan
%1.561419360929e9,-0.6792532205581665,nan
];

data_dirs = {'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190713_forbidden427_direct_det_narrow\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190710_forbidden427_direct_det_narrow_dither_on\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190710_forbidden427_direct_det\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden_Rf_2.00\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_weekend_run\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden_Rf_2.05\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190715_forbidden_427_narrow_scan_missed\'
    %%
    %'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190719_forbidden_Rf_none\'
    %'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190715_forbidden427_narrow_changed_pol_half_pow\'
    %%
    %'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190711_forbidden427_direct_det_bad_atom_num\'
    %'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190704_forbidden_long interrogation\'
    };
%%
data.signal = [];
data.freq = [];
data.time = [];
data.atom_num = [];
data.probe_num = [];
data.integrated_pd = [];
data.is_shot_good = [];

comp = [1,1,1,1,1,1,1,1,1];
% comp = [3.4,3.2,7.3,6.3,7.8,6.4,1];
% comp = [3.4,4.1,3.2,7.3,6.3,7.8,6.4];
run_time = [];
cen = [5.4,7.2,2.3,6.56,4.2,6.0,4.5,0,-7.8];
cen = [5.4,7.2,2.3,6.56,4.2,6.0,4.5,0];
cen = [5.4,7.2,2.3,6.56,4.2,6.0,4.5];
% cen = [5.4,2.3,6.56,4.2,6.0,4.5]; %dither removed
cen_unc = [0.2,0.2,0.2,0.12,0.2,0.3,0.13];
% cen = [1];
% cen_unc = [0.1];
for loop_idx=1:length(data_dirs)
    current_dir = data_dirs{loop_idx};
    fprintf('importing data from \n %s \n',current_dir)
    if ~strcmp(current_dir(end),'\')
        current_dir = [current_dir,'\'];
    end
    out_dirs=dir(fullfile(current_dir,'out'));
    out_dirs=out_dirs(3:end);
    out_dirs=out_dirs(cat(1,out_dirs.isdir));
    if size(out_dirs,1)==0
        warning(sprintf('dir \n %s \n does not contain any out dirs',current_dir)), 
    else 
        % convert the folder name (iso time) to posix time
        time_posix=cellfun(@(x) posixtime(datetime(datenum(x,'yyyymmddTHHMMSS'),'ConvertFrom','datenum')),{out_dirs.name});
        [~,sort_idx]=sort(time_posix,'descend');
        out_dirs=out_dirs(sort_idx);
        looking_for_data_dir=1;
        folder_index=1;
        %runs through all the out put dirs for a given run and looks for saved data, if none is there
        %skips that data dir
        while looking_for_data_dir
            try
                out_instance_folder_path=fullfile(current_dir,'out',out_dirs(folder_index).name,'data_results.mat');
                if (exist(out_instance_folder_path,'file') || ...
                    exist(out_instance_folder_path,'file')) && ...
                    exist(out_instance_folder_path,'file')
                    looking_for_data_dir=0;
                else
                    folder_index=folder_index+1;
                    if folder_index>numel(out_dirs) %if beyon the end of the folder list return nan;
                         looking_for_data_dir=0;
                         folder_index=nan;
                         warning('did not find a valid output direcory for folder %s',current_dir)
                    end
                end
                %~and(isfile([current_dir,'out\',most_recent_dir.name,'\main_data.mat']),isfile([current_dir,'out\',most_recent_dir.name,'\drift_data.mat']))
                %offset = offset + 1;
                %most_recent_dir=out_dirs(end-offset,1);
                %check = drift_data.avg_coef; %check if it has the avg coefs update

            catch e
                fprintf('\n dir: %s didnt work \n',current_dir)
                msgText = getReport(e)
                continue
            end
        end
        if ~isnan(folder_index)

            load(fullfile(current_dir,'out',out_dirs(folder_index).name,'data_results.mat'))
            % now do some serious data plumbing
            %append to main structure
            
            data.signal = cat(1,data.signal,out_data.signal.cal.calibrated_signal.val./comp(loop_idx));
            data.freq = cat(1,data.freq,out_data.freq);
            data.time = cat(1,data.time,out_data.time);
            data.atom_num = cat(1,data.atom_num,out_data.atom_num');
            data.probe_num = cat(1,data.probe_num,out_data.probe_num');
            data.integrated_pd = cat(1,data.integrated_pd,out_data.integrated_pd);
            data.is_shot_good = cat(1,data.is_shot_good,out_data.is_shot_good);
            run_time = [run_time,mean(out_data.time)];
%             if ~isequal(size(drift_data_compiled.to.val),size(drift_data_compiled.wp.qwp))
%                 error('things are not the same size') 
%             end

        end
    end
end
%% create drift model
data.is_shot_good = logical(data.is_shot_good);
% offset = interp1(wm_offset(:,1),wm_offset(:,2),data.time);%,'spline');
offset = 2.*interp1(wm_offset(:,1),wm_offset(:,2),data.time,'pchip');%spline
offset_run = 2.*interp1(wm_offset(:,1),wm_offset(:,2),run_time,'pchip');% 'makima'pchip 'spline' 'nearest'
% offset_run = smoothdata(interp1(wm_offset(:,1),wm_offset(:,2),run_time,'spline'),'sgolay',1);%
t_temp = linspace(min(wm_offset(:,1)),max(wm_offset(:,1)),4000);
offset_mdl =  2.*interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'pchip');% 'makima' pchip 'spline' 'nearest'
% offset_mdl =  smoothdata(interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'spline'),'sgolay',300);%
%plot drift model
stfig('wm drift model')
clf
%scatter(run_time,cen)
hold on
%plot(data.time,offset)
plot(t_temp,offset_mdl)
scatter(wm_offset(:,1),2.*wm_offset(:,2),'kx')
ylabel('wavemetre offset (MHz)')
xlabel('Posix time (s)')
box on
stfig('centers of distributions')
errorbar(run_time,cen-offset_run,cen_unc,'x')
xlabel('time')
ylabel('corrected cen')
cen_spread_offset = std(cen-offset_run)
cen_spread = std(cen)
%% set up ac stark shifts
Zeeman_shift = 0;%-1.7154;%3.16329586028859e-01,
ac_shift = 0;%2.95062819924799e-01.*data.integrated_pd(data.is_shot_good);
shifts = ac_shift;
%%
%set up the colors to use
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]; %[88,113,219]%[96,144,201]
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=20;

% bin data points and fit

predicted_freq=700939267; %MHz

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);

gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4);
loren_fun1d = @(b,x) b(1)./((x-b(2)).^2+b(3).^2) + b(4);
fun1d = loren_fun1d;
coeff_names={'amp','mu','sig','offset'};

xdata=data.freq(data.is_shot_good) - offset(data.is_shot_good) - shifts;
ydata=data.signal(data.is_shot_good);%./data.atom_num(data.is_shot_good);

pd_tol = 5;
sigma_from_mean=(ydata-nanmean(ydata,1))/nanstd(ydata,1);
is_not_oulier=sigma_from_mean<7 & data.integrated_pd(data.is_shot_good)>pd_tol;
ydata=ydata(is_not_oulier);% - predict(res_mdl_ratio,int_num./atom_num);
xdata=xdata(is_not_oulier);

%many bins wide
num_bins = 200;
probe_freq_bins = linspace(min(xdata),max(xdata),num_bins);
%few bins wide
num_bins = 26;
probe_freq_bins = linspace(min(xdata),max(xdata),num_bins);
%narrow many
% num_bins = 40;
% probe_freq_bins = linspace(-3,17,num_bins);



stfig('combined data')
clf
ylabel_str='Scattered fraction (arb. units)';



amp_guess=max(ydata);
ydata_shifted =ydata-min(ydata);
mu_guess=wmean(xdata,ydata_shifted); %compute the weighted mean
%sig_guess=sqrt(nansum((xdata-mu_guess).^2.*ydata_shifted)/nansum(ydata_shifted)); %compute the mean square weighted deviation
sig_guess=10;
offset_guess = 0;
amp_guess=4.556661697251370e+01;
mu_guess=1;
sig_guess=2.866375907039084e+00;
offset_guess=-3.340314780403577e-01;
fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
% 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
inital_guess=[amp_guess,mu_guess,sig_guess,offset_guess];
fitobject=fitnlm(xdata,ydata,...
    fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo);
fit_coeff=fitobject.Coefficients.Estimate;
fit_se=fitobject.Coefficients.SE;
cen_val =fit_coeff(2);
fitobject

%     amp_str=string_value_with_unc(fitobject.Coefficients.Estimate(1),fitobject.Coefficients.SE(1),'b');
%     cen_str=string_value_with_unc(fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2),'b');
%     width_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject.Coefficients.SE(3),'b');

amp_str=string_value_with_unc(fitobject.Coefficients.Estimate(1)/fitobject.Coefficients.Estimate(3)^2+fitobject.Coefficients.Estimate(4),fitobject.Coefficients.SE(1)/fitobject.Coefficients.Estimate(3)^2,'b');
cen_str=string_value_with_unc(fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2),'b');
width_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject.Coefficients.SE(3),'b');
offset_str=string_value_with_unc(fitobject.Coefficients.Estimate(4),fitobject.Coefficients.SE(4),'b');

width_units='MHz';
offset_units='counts';
amp_units='counts';
str=sprintf('Gauss fit \n   Cen    %s %s \n   Width %s %s \n   Amp   %s %s \n   Offset %s %s',...
    cen_str,width_units,width_str,width_units,amp_str,amp_units,offset_str,offset_units);
%     text(0.01,0.9,str,'Units','normalized');
%wide on edges many on peak
probe_freq_bins =[linspace(min(xdata),fitobject.Coefficients.Estimate(2)-8,8),...
    linspace(fitobject.Coefficients.Estimate(2)-6,fitobject.Coefficients.Estimate(2)+6,8),...
    linspace(fitobject.Coefficients.Estimate(2)+8,max(xdata),8)];
iimax=numel(probe_freq_bins)-1;
signal_bined.freq_std=nan(iimax,1);
signal_bined.val=nan(iimax,1);
signal_bined.unc_val=nan(iimax,1);
signal_bined.freq_mean=nan(iimax,1);
signal_bined.freq_obs_min_max_mean_diff=nan(iimax,2);
for ii=1:iimax
    signal_bined.freq_bin_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=xdata<=probe_freq_bins(ii+1) & xdata>probe_freq_bins(ii);
    signal_bined.freq_bin_cen(ii)=nanmean(probe_freq_bins(ii:ii+1));
    if sum(bin_mask)==0
        warning('no elements')
        signal_bined.num_bin(ii)=0;
    else
        signal_bined.num_bin(ii)=sum(bin_mask);
        signal_bined.val(ii,:)=nanmean(ydata(bin_mask,:),1);
        signal_bined.unc_val(ii,:)=nanstd(ydata(bin_mask,:),[],1)./sqrt(sum(bin_mask));
        signal_bined.freq_mean(ii)=nanmean(xdata(bin_mask));
        signal_bined.freq_std(ii)=nanstd(xdata(bin_mask));
        signal_bined.freq_obs_min_max(ii,:)=[min(xdata(bin_mask)),max(xdata(bin_mask))];
        signal_bined.freq_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_bin_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_obs_min_max_mean_diff(ii,:)=abs(signal_bined.freq_obs_min_max(ii,:)-signal_bined.freq_mean(ii));
    end
end

y_norm = max(signal_bined.val);
x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
[ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on
plot(x_sample_fit-cen_val,ysamp_val./y_norm,'k','LineWidth',1.5)
drawnow
yl=ylim;
plot(x_sample_fit-cen_val,ysamp_ci./y_norm,'color',[1,1,1].*0.5)

curve1 = ysamp_ci(:,1)';
curve2 = ysamp_ci(:,2)';
x1 = (x_sample_fit-cen_val)';
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween./y_norm, 'g');
h.FaceColor = [0.31 0.31 0.32].*2;
h.FaceAlpha = 0.5;

errorbar(signal_bined.freq_mean(3:end-1)-cen_val,signal_bined.val(3:end-1)./y_norm,...
    signal_bined.unc_val(3:end-1,1)./y_norm,signal_bined.unc_val(3:end-1,1)./y_norm,...
    signal_bined.freq_obs_min_max_mean_diff(3:end-1,1), signal_bined.freq_obs_min_max_mean_diff(3:end-1,2),...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5);
hold on
plot(signal_bined.freq_mean(2:end-1)-cen_val,signal_bined.val(2:end-1)./y_norm,'o','MarkerSize',5,'MarkerFaceColor',colors_detail(1,:),'MarkerEdgeColor',colors_main(2,:))

ylim(yl)
xlim([min(xdata),max(xdata)]-cen_val)
xlabel('\(f-f_{0,d}\) (MHz)','fontsize',font_size_global,'interpreter','latex')
% show the inital guess
%plot(x_sample_fit,gauss_fun1d(inital_guess,x_sample_fit)*ymultipler)
box on
fprintf('transition frequnency %s\n',string_value_with_unc(predicted_freq+fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2)))
set(gca,'fontsize',font_size_global)
ylabel(ylabel_str,'fontsize',font_size_global-0.8,'interpreter','latex')
xlim([-36.5,36.5])
ax = gca;
set(gca,'TickLabelInterpreter','latex')
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom-0.0036 ax_width ax_height];


fig = gcf;
set(fig,'Position',[1126 491 693 442])
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% print(fig,'C:\Users\kieran\Documents\MATLAB\Forbidden_Transition\figs\direct_scan','-dpdf')
    %%
    stfig('All data')
    plot(xdata,ydata,'x')
    hold on
    plot(probe_freq_bins,probe_freq_bins*0,'s')
        yl=ylim;
    ylim(yl)
    xlim([min(xdata),max(xdata)])
    xlabel('freq-theory (MHz)')
    ylabel(ylabel_str)
%% residuals
if 0
freq_window = [-60,60];
% freq_window = [-3,18];
% freq_window = [3,10];
% freq_window = [-60,-10];
% freq_window = [19,60];
        
freq_mask=xdata<=freq_window(2) & xdata>freq_window(1);
int_pd = data.integrated_pd(data.is_shot_good);
int_pd = int_pd(is_not_oulier);
int_num = data.probe_num(data.is_shot_good);
int_num = int_num(is_not_oulier);
atom_num = data.atom_num(data.is_shot_good);
atom_num = atom_num(is_not_oulier);
ypred_val=predict(fitobject,xdata,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
res = ydata - ypred_val;
sfigure(345);
clf
plot(res)
xlabel('residuals')
ylabel('shot indx')
sfigure(456);
clf
res_mdl = corr_plot(int_num(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('background counts')
sfigure(567);
clf
corr_plot(atom_num(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('total atom counts')
sfigure(678);
clf
corr_plot(xdata(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('freq - theory (MHz)')
sfigure(123)
clf
corr_plot(int_pd(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('integrated pd voltage')
sfigure(1234)
clf
res_mdl_ratio = corr_plot(int_num(freq_mask)./atom_num(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('atom num ratio')
sfigure(8012)
clf
corr_plot(int_num(freq_mask)./int_pd(freq_mask)./atom_num(freq_mask),res(freq_mask),ones(size(res(freq_mask))))
ylabel('res')
xlabel('test')
end
