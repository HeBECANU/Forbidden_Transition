function offset = wm_drift_model(time,wm_offset,plot_wm_model)
%generates drift model for wm
offset = 2.*interp1(wm_offset(:,1),wm_offset(:,2),time,'pchip');


%plot drift model
if plot_wm_model
    t_temp = linspace(min(wm_offset(:,1)).*0.99995,max(wm_offset(:,1)).*1.0002,4000);
    offset_mdl =  2.*interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'pchip');%
    wm_err = 4.1.*ones(size(t_temp));
    wm_m = 2*mean(wm_offset(:,2)).*ones(size(t_temp));
    stfig('wm drift model')
    clf
    hold on
    %plot(t_temp/1e9,offset_mdl,'linewidth',1.5)
    unit=60*60*24;
    
    curve1 = (wm_m+wm_err);
curve2 = (wm_m-wm_err);
x1 = ((t_temp-min(wm_offset(:,1)))./unit);
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, 'g');
h.FaceColor = [0.31 0.31 0.32].*2;
h.FaceAlpha = 0.4;
    
    plot((t_temp-min(wm_offset(:,1)))./unit,wm_m+wm_err,'k','linewidth',1.5)
    plot((t_temp-min(wm_offset(:,1)))./unit,wm_m,'r','linewidth',1.5)
    plot((t_temp-min(wm_offset(:,1)))./unit,wm_m-wm_err,'k','linewidth',1.5)
    scatter((wm_offset(:,1)-min(wm_offset(:,1)))/unit,2.*wm_offset(:,2),75,'kx')
    ylabel('Wavemeter offset (MHz)','interpreter','latex')
    xlabel('Time after initial calibration (days)','interpreter','latex')
    set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',15)
xlim([-0.1,13])
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
    box on
    t_period = 3*30*3600;%seconds ie one hour
n_period = floor(t_period*length(t_temp)/range(t_temp));
offset_std = [];
% for ii = 1:(length(offset_mdl)-n_period)
%     offset_std = [offset_std,std(offset(ii:(ii+n_period)))];
% end
% sfigure(45)
% plot(t_temp(1:end-n_period),offset_std)
% mean(offset_std)
% std(offset_std)
end
end
