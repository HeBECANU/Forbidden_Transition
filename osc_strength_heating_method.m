%oscillator strength via the hating method

hebec_constants
kb = const.kb;
mhe = const.mhe;
h = const.h;
c = const.c;
alpha = 3*8*pi*h*kb/c^2;
%trap freq
wx=2*pi*40;
wy=2*pi*420;
wz=2*pi*40;
predicted_freq=700939267; %MHz

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


%import heating data
load('X:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden427_overnight_heating_method\out\20190724T095142\data_results.mat')%20190802T130932,20190718T110943,20190830T144824,20190724T095142
Ti = out_data.data.signal.msr.val(:,2);
Tf = out_data.data.signal.msr.val(:,3);
T=(Ti+Tf)./2;
dT_dt = out_data.data.cal.calibrated_signal.val(:);%out_data.data.signal.msr.val(:,1);
f = out_data.data.signal.msr.freq;

Ry= sqrt(2*kb.*T./(mhe*wy^2));
Rz= sqrt(2*kb.*T./(mhe*wz^2));

Ep = 4/3*const.h.*f;%average kinetic energy added by each photon
%convert pd voltage to power

P_ratio = 1- const.h;

sfigure(1001);
plot(Ry.*Rz)
ylabel('Ry*Rz')

%%
convert = 1;%0.9044.*3.2616e-2;
gauss_fun1d = @(b,x) b(1).*exp(-((x-b(3)).^2)./(2*b(2).^2));
loren_fun1d = @(b,x) b(1)./((x-b(3)).^2+b(2).^2);%+b(4);
coeff_names={'amp','sig','mu'};

fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
% 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
inital_guess=[4e-4,0.1,4];
fitobject=fitnlm(f,dT_dt./(convert),...
    gauss_fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo)
inital_guess_l=[4e-3,10,4];
fitobject_l=fitnlm(f,dT_dt./(convert),...
    loren_fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo)
fit_coeff_l=fitobject_l.Coefficients.Estimate;
sfigure(1002);
clf
fit_coeff=fitobject.Coefficients.Estimate;
fit_se=fitobject.Coefficients.SE;
x_sample_fit=col_vec(linspace(min(f),max(f),1e3));
% [ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on
% plot(x_sample_fit,ysamp_val,'r')
% drawnow
% yl=ylim;
% plot(x_sample_fit,ysamp_ci,'color',[1,1,1].*0.5)
[ysamp_val,ysamp_ci]=predict(fitobject_l,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
plot(x_sample_fit,ysamp_val,'b')
plot(x_sample_fit,ysamp_ci,'color',[1,1,1].*0.2)
scatter(f,dT_dt./(convert))
xlabel('freq')
ylabel('signal')
Ep = 3.888*10^(-10);
nanmean(Ry.*Rz)*9.48*10^(-9)*sqrt(pi)*fit_coeff(1)*sqrt(2)*fit_coeff(2)*10^6*1/Ep;
nanmean(Ry.*Rz)*9.48*10^(-9)*pi*fit_coeff_l(1)/fit_coeff_l(2)*10^6*1/Ep;
nanmean(Ry.*Rz)*3.576*10^(-8)*pi*fit_coeff_l(1)/fit_coeff_l(2)*10^6*1/Ep;
Int_1 = (predicted_freq+fit_coeff_l(3))*1e6.*pi*fit_coeff_l(1)/fit_coeff_l(2)
%%
stfig('heating data')
clf
cen_val = fit_coeff_l(3);
ylabel_str='Change in heating rate (nK/s)';
xdata = f;
ydata= dT_dt.*1e9;
    probe_freq_bins =[linspace(min(xdata),max(xdata),8)];
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
hold on
      
x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
    [ysamp_val,ysamp_ci]=predict(fitobject_l,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    curve1 = (ysamp_ci(:,1).*1e9)';
curve2 = (ysamp_ci(:,2).*1e9)';
x1 = (x_sample_fit-cen_val)';
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, 'g');
h.FaceColor = [0.31 0.31 0.32].*2.5;
h.FaceAlpha = 0.5;
    hold on
    plot(x_sample_fit-cen_val,ysamp_val.*1e9,'k','LineWidth',1.5)
    %drawnow
    yl=ylim;
    plot(x_sample_fit-cen_val,ysamp_ci.*1e9,'color',[1,1,1].*0.5)
    ylim([0,2.8])
    xlim([min(xdata),max(xdata)]-cen_val)
    xlabel('\(f-f_{0,h}\) (MHz)','fontsize',14,'interpreter','latex')
       % show the inital guess
    %plot(x_sample_fit,gauss_fun1d(inital_guess,x_sample_fit)*ymultipler)
     box on
    fprintf('transition frequnency %s\n',string_value_with_unc(predicted_freq+fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2)))
   set(gca,'fontsize',font_size_global)
    ylabel(ylabel_str,'fontsize',19.5,'interpreter','latex')
   xlim([-7,13])
  errorbar(signal_bined.freq_mean-cen_val,signal_bined.val,...
        signal_bined.unc_val(:,1),signal_bined.unc_val(:,1),...
         signal_bined.freq_obs_min_max_mean_diff(:,1), signal_bined.freq_obs_min_max_mean_diff(:,2),...
        'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
         'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5);
hold on    
plot(signal_bined.freq_mean-cen_val,signal_bined.val,'o','MarkerSize',5,'MarkerFaceColor',colors_detail(1,:),'MarkerEdgeColor',colors_main(2,:))
set(gca,'TickLabelInterpreter','latex')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


fig = gcf;
set(fig,'Position',[1126 491 693 442])
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'C:\Users\kieran\Documents\MATLAB\Forbidden_Transition\figs\heating_scan','-dpdf')
