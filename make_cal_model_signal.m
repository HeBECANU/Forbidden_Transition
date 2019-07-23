function out=make_cal_model(anal_opts_cal,data)
%creates a model of what the trap freq of the calibration shots does over time
% TODO:
% - [ ] make sure the time is using the probe time not the labview trig time


%create a mask of shots that are calibrations and are good fits
temp_cal=col_vec(data.mcp_tdc.probe.calibration);
temp_cal(isnan(temp_cal))=0;
cal_dat_mask=col_vec(data.mcp_tdc.ok.all) & temp_cal;

if sum(cal_dat_mask)<3
    error('not enough calibration data pts')
end

temp_cal=col_vec(data.mcp_tdc.probe.calibration);
temp_cal(isnan(temp_cal))=1;
probe_dat_mask=col_vec(data.mcp_tdc.ok.all) & ~temp_cal;
time_probe=col_vec(data.mcp_tdc.time_create_write(probe_dat_mask,2));


%create a time vector
time_cal=col_vec(data.mcp_tdc.time_create_write(cal_dat_mask,2));
time_start_cal=time_cal(1); %shift the times to start at zero for better interp perfromance
time_cal=time_cal-time_start_cal;
%create a vector of the reconstructed calibration frequencies
signal_cal=col_vec(data.signal.raw.val(cal_dat_mask));

%interpolate the input data to subsample by a factor of 100
time_samp_interp=col_vec(linspace(min(time_cal),max(time_cal),size(time_cal,1)*30));
%signal_interp_raw=col_vec(interp1(time_cal,signal_cal,time_samp_interp,'pchip'));
signal_interp_raw=col_vec(interp1(time_cal,signal_cal,time_samp_interp,'linear'));

%% do some filtering of the data
%smooth the interpolated data
%BMH 20190523 change - moved to using gausfilt
%yinterp_smooth=gaussfilt(time_samp_interp,trap_freq_interp_raw,anal_opts_cal.smooth_time);
smoothing_window_len=round(5*anal_opts_cal.smooth_time/mean(diff(time_samp_interp)));
yinterp_smooth=smoothdata(signal_interp_raw,'sgolay',smoothing_window_len);


%% create a model based on interpolating the smoothed/filtered interpolated data
out.signal_cal_drift_model=@(x) interp1(time_samp_interp,yinterp_smooth,x-time_start_cal,'pchip','extrap'); %'linear'
out.num_shots=sum(cal_dat_mask);
out.cal_mask=cal_dat_mask;

%calculate the residuals between the model and the calibration data
model_resid=out.signal_cal_drift_model(time_cal+time_start_cal)-signal_cal;

if mean(model_resid)>std(model_resid)
    error('error mean of residuals is not within 1sd')
end
    
out.unc=nanstd(model_resid);
% shot noise limited atom counting
mean_cal_shot_unc=mean(data.signal.raw.unc(cal_dat_mask));

fprintf('%s:residual std %f vs model mean uncert %f \n',...
     mfilename,nanstd(model_resid),mean_cal_shot_unc)


if anal_opts_cal.plot
    %%
    hour_in_s=60*60;
    x_samp=linspace(min(time_cal)-30,max(time_cal)+30,size(time_cal,1)*1e2);
    stfig('osc cal model','add_stack',1);
    clf
    subplot(3,1,1)
    errorbar(data.mcp_tdc.shot_num(cal_dat_mask),...
        data.signal.raw.val(cal_dat_mask),data.signal.raw.unc(cal_dat_mask)...
        ,'.k-','capsize',0,'MarkerSize',10,'linewidth',1.0)%,'r.',
    hold on
        errorbar(data.mcp_tdc.shot_num(probe_dat_mask),...
        data.signal.raw.val(probe_dat_mask),data.signal.raw.unc(probe_dat_mask)...
        ,'.b-','capsize',0,'MarkerSize',10,'linewidth',1.0)%,'r.',
    hold off
    legend('calibration data','probe data')
    xlabel('Shot Number')
    title('Input Data')
    ylabel('Measured Trap Freq (Hz)')
    
    %plot the calibration model
    subplot(3,1,2)
    plot(time_cal/hour_in_s,signal_cal,'x','MarkerSize',20)
    hold on
    plot(time_samp_interp/hour_in_s,signal_interp_raw,'m')
    plot(x_samp/hour_in_s,out.signal_cal_drift_model(x_samp+time_start_cal),'k')
    legend('data','interp data','calibration model' )
    hold off
    title('Trap Freq Calibration Model')
    xlabel('experiment time (h)')
    ylabel('no probe trap freq')
    %set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    %plot the calibration model
    
    subplot(3,1,3)
    errorbar(time_cal/hour_in_s,...
       model_resid,data.signal.raw.unc(cal_dat_mask)...
        ,'.k-','capsize',0,'MarkerSize',10,'linewidth',1.0)%,'r.',
    xl=[min(time_cal),max(time_cal)]./hour_in_s;
    line(xl,[1,1]*mean(model_resid),'color','g')
    line(xl,[1,1]*mean(model_resid)+std(model_resid),'color','r')
    line(xl,[1,1]*mean(model_resid)-std(model_resid),'color','r')
    xlabel('experiment time (h)')
    ylabel('residuals')
    %set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    %%
end

%% now we produce a calibrated signal value
calibrated_signal.val=data.signal.raw.val*nan; %initalize
calibrated_signal.unc=data.signal.raw.val*nan; %initalize

calibrated_signal.val(probe_dat_mask)=(data.signal.raw.val(probe_dat_mask)-out.signal_cal_drift_model(time_probe));%./data.ai_log.pd.integrated(probe_dat_mask);
calibrated_signal.unc(probe_dat_mask)=data.signal.raw.unc(probe_dat_mask);

out.calibrated_signal=calibrated_signal;


end


% code for ploting atom number    
%     subplot(2,1,2)
%     num_t =data.mcp_tdc.time_create_write(cal_dat_mask,2)-data.mcp_tdc.time_create_write(1,2);
%     plot(num_t/(60*60),data.mcp_tdc.num_counts(cal_dat_mask))
%     saveas(gcf,[anal_opts_cal.global.out_dir,plot_name,'.png'])
%     saveas(gcf,[anal_opts_cal.global.out_dir,plot_name,'.fig'])
%     title('Hit count trend')
%     xlabel('Time (h)')
%     ylabel('counts')
