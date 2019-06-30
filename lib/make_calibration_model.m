function calib  = make_calibration_model(data,opts)    
    
    header({0,'Making calibration model'})
    calib = data.sync.msr;
    %signal = data.sync.msr.N_atoms;
%data.sync.cal.N_atoms

    % Plotting the calibration model
        
    plot_raw_T = data.sync.msr.tdc_time-data.sync.start_time;
    plot_raw_X = data.sync.msr.probe_set;
    plot_raw_Y = data.signal.raw;
    
    
    plot_cal_T = data.sync.cal.tdc_time-data.sync.start_time;
    plot_cal_X = data.sync.cal.probe_set;
    plot_cal_Y = data.signal.cal;
    

%     calib = interp1(data.sync.cal.tdc_time,data.sync.cal.N_atoms,data.sync.msr.tdc_time);     
%     plot_mdl_Y = signal./calib;
%     
    % Cull out shots with too few counts based on interpolation of calibration
    mdl_mask = interp1(data.sync.cal.tdc_time,data.sync.cal.N_atoms,data.sync.msr.tdc_time);
    out_mask = ~isnan(data.signal.cal);
    calib.mdl = interp1(data.sync.cal.tdc_time(out_mask),data.signal.cal(out_mask),data.sync.msr.tdc_time);     
    calib.signal = -calib.mdl+data.signal.raw;
    calib.low_count_mask = mdl_mask > opts.check.min_counts;
    calib.probe_set = data.sync.msr.probe_set;
%     calib.signal = calib.signal(calib.low_count_mask);

    
%     sync_data.msr.calib = data.cal.mdl;
%     sync_data.msr.diff = data.sync.msr.calib - signal;
%     sync_data.msr.normdiff = 1-(data.sync.msr.calib - signal)./data.sync.msr.calib;
%     sync_data.msr.ratio = signal./data.sync.msr.calib;
    
    plot_mdl_Y = calib.signal(calib.low_count_mask);
    plot_mdl_X = data.sync.msr.probe_set(calib.low_count_mask);
    
    f1=sfigure(5000);
    clf;
    subplot(3,1,1)
    plot(plot_raw_T,plot_raw_Y,'.')
    hold on
    plot(plot_cal_T,plot_cal_Y,'.')
    plot(plot_raw_T,calib.mdl)
%     plot(plot_raw_T,calib.mdl)
    xlabel('Time elapsed')
    ylabel('Signal')
    legend('Measurement shots','calibration shots','Model')
    title('Raw data')

    
    subplot(3,1,2)
    shot_num_temp = 1:length(calib.signal);
    plot(calib.signal,'.')
    hold on
    plot(shot_num_temp(~calib.low_count_mask),calib.signal(~calib.low_count_mask),'rx')
    xlabel('Shot number')
    ylabel('Signal')
    title('Calibrated time series')    
    legend('Raw data','N<min')
    
    
    subplot(3,1,3)

    plot(plot_mdl_X-mean(plot_mdl_X),plot_mdl_Y,'.')

    hold on

    nbins=50;
    X_bin=linspace(min(plot_mdl_X),max(plot_mdl_X),nbins);
    X_bin_width=(-min(plot_mdl_X)+max(plot_mdl_X))/nbins;
    Y_val = zeros(nbins,2);
    for ii = 1:nbins
       bin_cen = X_bin(ii);
       indx = abs(plot_mdl_X-bin_cen)<X_bin_width/2;
       Y_val(ii,1) = mean(plot_mdl_Y(indx));
       Y_val(ii,2) = std(plot_mdl_Y(indx))./sqrt(sum(indx));
    end
    errorbar(X_bin-mean(plot_mdl_X),Y_val(:,1),Y_val(:,2),'.');
    yy = gaussfilt(X_bin,Y_val(:,1),7);
    plot(X_bin-mean(plot_mdl_X),yy)
    xlabel(sprintf('f-%u (MHz)',mean(plot_mdl_X)))
    ylabel('Signal')

    title('Calibrated spectra')    
    legend('Raw data','Binned data','Smoothed Response')
    grid on
    filename2 = fullfile(opts.out_dir,sprintf('%s_diagnostic',mfilename));
    saveas(f1,[filename2,'.fig']);
    saveas(f1,[filename2,'.png'])
    header({1,'Done.'})
end