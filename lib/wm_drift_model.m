function offset = wm_drift_model(time,wm_offset,plot_wm_model)
%generates drift model for wm
offset = 2.*interp1(wm_offset(:,1),wm_offset(:,2),time,'pchip');


%plot drift model
if plot_wm_model
    t_temp = linspace(min(wm_offset(:,1)),max(wm_offset(:,1)),4000);
    offset_mdl =  2.*interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'pchip');%
    stfig('wm drift model')
    clf
    hold on
    plot(t_temp/1e9,offset_mdl,'linewidth',1.5)
    scatter(wm_offset(:,1)/1e9,2.*wm_offset(:,2),75,'kx')
    ylabel('Wavemetre offset (MHz)','interpreter','latex')
    xlabel('Posix time (\(10^9\)s)','interpreter','latex')
    set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',15)
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
