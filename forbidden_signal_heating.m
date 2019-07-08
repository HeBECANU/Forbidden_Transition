function output=forbidden_signal_heating(data,anal_opts)
% todo
% make a calibration model and subtract off the difference
global const


%%
%anal_opts.global.fall_velocity=const.g0*anal_opts.global.fall_time;

thermal_fit=[];
plot_therm_fit=1; 
num_shots=size(tdc_masked_shots,2);
mean_first_few=nan(num_shots,2);
idx_range_to_avg=[1:40];

single_shot_temp_fig=stfig('single shot temp. curve');

fprintf('fitting temp. rise in shots %04u:%04u',num_shots,0)
for shot_idx=1:size(al_pulse_data.num_counts,1)
    pulse_time=al_pulse_data.time_cen;
    pulse_time=pulse_time-pulse_time(1);
    
    temperature_from_time=al_pulse_data.fit.temperature.val(shot_idx,:);
    temperature_uncert=al_pulse_data.fit.temperature.unc(shot_idx,:);
    
    %temperature_from_time=(anal_opts.global.fall_velocity*al_pulse_data.pos.std(shot_idx,:,1)/anal_opts.global.fall_time).^2 *const.mhe/const.kb;
    % if we assume a normal distribution we can estimate the error in the estimated std
    % std(sigma)= \[Sigma] Sqrt[1 - c4^2]
    %std_in_std=data.mcp_tdc.al_pulses.pos.std(5,:,1).*sqrt(1-normal_correction_c4(data.mcp_tdc.al_pulses.num_counts(shot_number,:)).^2);
    %approx way
    %frac_std_in_std=1./sqrt(2*(al_pulse_data.num_counts(shot_idx,:)-1));
    % the power makes the fractional uncert a factor of 2 more
    %temperature_uncert=frac_std_in_std.*2.*temperature_from_time;

    xdata=col_vec(pulse_time);
    ydata=col_vec(temperature_from_time);
    yunc=col_vec(temperature_uncert);
    
    [mean_first_few(shot_idx,1),mean_first_few(shot_idx,2)]=unc_wmean(ydata(idx_range_to_avg),yunc(idx_range_to_avg));
    
    if al_mcp_data.all_ok(shot_idx) && sum(~isnan(yunc) & ~isnan(ydata) & ~isnan(xdata))>10
% %         modelfun= @(b,x) b(1)+x.*b(2)+x.^2.*b(3);
% %         beta0=[0.2,0.2,0];
% %         cof_names={'offset','lin','square'};

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
            temp_multipler_plot=1e-6;
            stfig(single_shot_temp_fig);
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
            pause(0.001)
        end
      
    else
         thermal_fit.RMSE(shot_idx)=nan;
         thermal_fit.fit_coeff.SE(shot_idx,:)=[nan,nan];
         thermal_fit.fit_coeff.Estimate(shot_idx,:)=thermal_fit.fit_coeff.SE(shot_idx,:);
    end
    
    if mod(shot_idx,10)==0,fprintf('\b\b\b\b%04u',shot_idx),end   
end
fprintf('...Done\n') 


%%
predicted_freq=700939267; %MHz %from ma boi drake

is_cal=data.mcp_tdc.probe.calibration;
%is_cal(isnan(is_cal))=1;
is_cal(isnan(is_cal))=0;

is_freq_good=data.wm_log.proc.probe.freq.act.std<5 &...
    (data.wm_log.proc.probe.freq.set-data.wm_log.proc.probe.freq.act.mean)<5 &...
    ~is_cal;

probe_freq=data.wm_log.proc.probe.freq.act.mean*2;%freq in blue
probe_freq=probe_freq-predicted_freq;

signal_unbinned.val=cat(2,thermal_fit.fit_coeff.Estimate(is_freq_good,2),thermal_fit.fit_coeff.Estimate(is_freq_good,1),mean_first_few(is_freq_good,1));
signal_unbinned.names={'fit grad','fit inital term','mean of temp for first few'};
signal_unbinned.ystr={'heating rate (nk/s)','temp. (uk)','temp. (uk)'};
signal_unbinned.ymult=[1e9,1e7,1e7];
signal_unbinned.freq=probe_freq(is_freq_good);


signal_bined=[];
probe_freq_bins=col_vec(linspace(-80,80,60));
%probe_freq_bins=col_vec(linspace(-25,20,9));
iimax=numel(probe_freq_bins)-1;
for ii=1:iimax
    signal_bined.freq_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=signal_unbinned.freq<probe_freq_bins(ii+1) & signal_unbinned.freq>probe_freq_bins(ii);
    signal_bined.val(ii,:)=nanmean(signal_unbinned.val(bin_mask,:),1);
    sum(bin_mask)
    std(signal_unbinned.val(bin_mask))
    signal_bined.unc_val(ii,:)=nanstd(signal_unbinned.val(bin_mask,:),[],1)./sqrt(sum(bin_mask));
    
    signal_bined.freq(ii)=nanmean(probe_freq_bins(ii:ii+1));
    signal_bined.freq_lims_diff(ii,:)=abs(signal_bined.freq_lims(ii,:)-signal_bined.freq(ii));
end

stfig('counts vs probe freq')
clf

tot_plots=3;
signal_idx=1;
for signal_idx=1:size(signal_unbinned.val,2)
ymultipler=signal_unbinned.ymult(signal_idx);
ylabel_str=signal_unbinned.ystr(signal_idx);
subplot(2,tot_plots,0+signal_idx)
plot(signal_unbinned.freq,signal_unbinned.val(:,signal_idx)*ymultipler,'x')
xlabel('freq-theory (MHz)')
%ylabel('heating rate (nk/s)')
ylabel(ylabel_str)
title(signal_unbinned.names{signal_idx})
xlim([min(probe_freq_bins),max(probe_freq_bins)])
xl=xlim;
subplot(2,tot_plots,tot_plots+signal_idx)
errorbarxy(col_vec(signal_bined.freq),signal_bined.val(:,signal_idx)*ymultipler,...
    signal_bined.freq_lims_diff,signal_bined.unc_val(:,signal_idx)*ymultipler)
xlim(xl)
xlabel('freq-theory (MHz)')
ylabel(ylabel_str)
title(signal_unbinned.names{signal_idx})
end



%%






%%



output=[];



end
