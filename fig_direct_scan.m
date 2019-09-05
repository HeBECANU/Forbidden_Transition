%generate the direct scan figure
%analysis options
%drift model 'makima' pchip 'spline' 'nearest'?
%plot drift model?
plot_wm_model = 0;
ylabel_str='\textbf{Scattered fraction (arb. units)}';
xlabel_str='\textbf{\boldmath{\(f-f_{0,d}\)} (MHz)}';
%set up the colors to use
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]; %[88,113,219]%[96,144,201]
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=20;

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);
%which model to use

%% Set up configureations for the direct det analysis
config_direct_det

%% import the preanalysied data from the give directories
cli_format_text('','c',3)
cli_format_text('IMPORTING DATA','c',3)
cli_format_text('','c',3)
data = import_data(data_dirs);

%% create drift model
data.offset = wm_drift_model(data,wm_offset,plot_wm_model);

%% Apply ac stark shift
% crete the ac stark shift model
wnlm = ac_stark_shift(0);
% generate the actual shift
ac_gradient = wnlm.Coefficients.Estimate(2);
data.ac_shift = ac_gradient.*data.integrated_pd(data.is_shot_good);

%% Fit model
fitobject = fit_data(data,fun_type);
fit_coeff=fitobject.Coefficients.Estimate;
fit_se=fitobject.Coefficients.SE;

%% Plot final results
% import physical constants
stfig('combined data');
clf
predicted_freq=700939267; %MHz
cen_val =fit_coeff(2);

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
xlabel(xlabel_str,'fontsize',font_size_global,'interpreter','latex')
box on
set(gca,'fontsize',font_size_global)
ylabel(ylabel_str,'fontsize',font_size_global-0.8,'interpreter','latex')
xlim([-36.5,36.5])
ax = gca;
set(gca,'TickLabelInterpreter','latex','FontWeight','bold')
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