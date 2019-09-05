function offset = wm_drift_model(data,wm_offset,plot_wm_model)
%generates drift model for wm
offset = 2.*interp1(wm_offset(:,1),wm_offset(:,2),data.time,'pchip');


%plot drift model
if plot_wm_model
    t_temp = linspace(min(wm_offset(:,1)),max(wm_offset(:,1)),4000);
    offset_mdl =  2.*interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'pchip');%
    stfig('wm drift model')
    clf
    hold on
    plot(t_temp,offset_mdl)
    scatter(wm_offset(:,1),2.*wm_offset(:,2),'kx')
    ylabel('wavemetre offset (MHz)')
    xlabel('Posix time (s)')
    box on
end
end
