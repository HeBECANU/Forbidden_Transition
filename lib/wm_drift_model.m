function offset = wm_drift_model(time,wm_offset,plot_wm_model)
%generates drift model for wm
offset = 2.*interp1(wm_offset(:,1),wm_offset(:,2),time,'pchip');

 set(gca, 'FontName', 'cmr10')
  set(gca,  'defaultAxesFontName', 'cmr10')
   set(gca,  'defaultTextFontName', 'cmr10')


 
%plot drift model
if plot_wm_model
    time_padd=1e5;
    t_temp = linspace(min(wm_offset(:,1))-time_padd,max(wm_offset(:,1))+time_padd,4000);
    offset_mdl =  2.*interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'pchip');%
    wm_err = 4.1.*ones(size(t_temp));
    wm_m = 2*mean(wm_offset(:,2)).*ones(size(t_temp));
    stfig('wm drift model')
    clf
    hold on
    %plot(t_temp/1e9,offset_mdl,'linewidth',1.5)
    time_unit=60*60*24;
    
    curve1 = (wm_m+wm_err);
curve2 = (wm_m-wm_err);
x1 = ((t_temp-min(wm_offset(:,1)))./time_unit);
% inBetween = [curve1, fliplr(curve2)];
% h = fill(x2, inBetween, 'g');
% h.FaceColor = [0.31 0.31 0.32].*2.5;
% h.FaceAlpha = 0.4;
    color_shaded=[0.89 0.9 0.9];
    patch([x1, fliplr(x1)],...
            [curve1, curve2], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    
    plot((t_temp-min(wm_offset(:,1)))./time_unit,wm_m+wm_err,'k','linewidth',1.5)
    plot((t_temp-min(wm_offset(:,1)))./time_unit,wm_m,'r','linewidth',1.5)
    plot((t_temp-min(wm_offset(:,1)))./time_unit,wm_m-wm_err,'k','linewidth',1.5)
    plot((wm_offset(:,1)-min(wm_offset(:,1)))/time_unit,2.*wm_offset(:,2),'ko',...
        'MarkerSize',8,...
        'LineWidth',1.5,...
        'MarkerEdgeColor',[0.3,0.3,0.8],...
        'MarkerFaceColor',[0.45,0.45,1])
    ylabel('Wavemeter offset (MHz)','interpreter','latex')
    xlabel('Time after initial calibration (days)','interpreter','latex')
    set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',15)
xlim([-0.5,13])
ax = gca;
% ax.XAxis.TickLabelFormat= '\\textbf{%g}';
% ax.YAxis.TickLabelFormat= '\\textbf{%g}';
set(gca,'TickLength',[0.02, 0.02])
xticks(0:4:12)
yticks(-8:2:2)
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
