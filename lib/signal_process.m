function sig = signal_process(data,opts)
%% Total atom number
sig.raw = data.sync.msr.N_atoms;
sig.cal = data.sync.cal.N_atoms;
%% Atom Number in a given window
% temp = [];
% xmin = -19e-3;
% xmax = 15e-3;
% ymin = -30e-3;
% ymax = 30e-3;
% t_window_min = 6.4;
% t_window_max = 6.9;
% 
% for c_data = data.tdc.counts_txy
%     tdc_data = cell2mat(c_data);
%     mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
%         & tdc_data(:,2)>xmin & tdc_data(:,2)<xmax ...
%         & tdc_data(:,3)>ymin & tdc_data(:,3)<ymax ;
%     data_windowed = tdc_data(mask,:);
%     atom_num = size(data_windowed,1); %atom number in given window
%     temp = [temp,atom_num];%num_dif./atom_num
% end
% sig.cal = temp(data.sync.cal.mask);
% sig.raw = temp(data.sync.msr.mask);
%% Trap frequency fit
% data.tdc.all_ok = ones(size(data.tdc.shot_num,2),1);
% data.tdc.al_pulses=bin_al_pulses(opts.atom_laser,data);
% data.osc_fit=fit_trap_freq(opts.osc_fit,data);
% data.osc_fit.trap_freq_recons=nan*data.osc_fit.ok.did_fits;
% mask=data.osc_fit.ok.did_fits;
% data.osc_fit.trap_freq_recons(mask)=3*(1/opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(mask,2,1);
% 
% data.osc_fit.trap_freq_recons(isoutlier(data.osc_fit.trap_freq_recons)) = nan;
% 
% sig.cal = data.osc_fit.trap_freq_recons(data.sync.cal.mask).';
% sig.raw = data.osc_fit.trap_freq_recons(data.sync.msr.mask).';
%% Simple fits (variance in out coupled distrabution)
% temp = [];
% 
% xmin = -25e-3;
% xmax = 25e-3;
% ymin = -30e-3;
% ymax = 35e-3;
% xybins = 500;
% t_window_min = 6.4;
% t_window_max = 9.5;
% 
% for c_data = data.tdc.counts_txy
%     tdc_data = cell2mat(c_data);
%     mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
%         & tdc_data(:,2)>xmin & tdc_data(:,2)<xmax ...
%         & tdc_data(:,3)>ymin & tdc_data(:,3)<ymax ;
%     data_windowed = tdc_data(mask,:);
%     X_bin=linspace(xmin,xmax,xybins);
%     X_bin_width=(xmax-xmin)/xybins;
%     [X1d_counts,X1d_edges]=histcounts(data_windowed(:,2),X_bin);
%     X1d_centers=mean([X1d_edges(1:end-1);X1d_edges(2:end)]);
%     X1d_counts=1e-3*X1d_counts/(X_bin_width);%why 1e-3 factor?
% 
%     amp_guess=max(X1d_counts);
%     mu_guess=sum(X1d_centers.*X1d_counts)/sum(X1d_counts); %compute the weighted mean
%     sig_guess=sum((X1d_centers-mu_guess).^2.*X1d_counts)/sum(X1d_counts); %compute the mean square weighted deviation
%     
% %     temp = [temp,sig_guess];
%     temp = [temp,sig_guess];
%     %pause(0.5)
% end
% sig.cal = temp(data.sync.cal.mask);
% sig.raw = temp(data.sync.msr.mask);
%% Rectangular number difference
% temp = [];
% 
% xmin = -19e-3;
% xmax = 15e-3;
% ymin = -30e-3;
% ymax = 30e-3;
% t_window_min = 0.5;
% t_window_max = 6.3;
% 
% for c_data = data.tdc.counts_txy
%     tdc_data = cell2mat(c_data);
%     mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
%         & tdc_data(:,2)>xmin & tdc_data(:,2)<xmax ...
%         & tdc_data(:,3)>ymin & tdc_data(:,3)<ymax ;
%     data_windowed = tdc_data(mask,:);
%     num_dif = sum(data_windowed(:,2)>1.2e-3) - sum(data_windowed(:,2)<1.2e-3);%x cut
% %     num_dif = sum(data_windowed(:,3)>-0e-3) - sum(data_windowed(:,3)<-0e-3);
% %     size(data_windowed,1) %atom number in given window
%     temp = [temp,num_dif];
% 
% end
% sig.cal = temp(data.sync.cal.mask);
% sig.raw = temp(data.sync.msr.mask);
%% Circular number difference
temp = [];
cen = [1.5,4.9].*1e-3;
r = 30.5e-3;
theta = pi/4;%the angle about which to partiction
t_window_min = 0.5;
t_window_max = 6.3;

for c_data = data.tdc.counts_txy
    tdc_data = cell2mat(c_data);
    mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
        & sqrt((tdc_data(:,2)-cen(1)).^2+(tdc_data(:,3)-cen(2)).^2)<r;
    data_windowed = tdc_data(mask,:);
    angs = atan(data_windowed(:,3)./data_windowed(:,2));
    atom_num = size(data_windowed,1); %atom number in given window
    num_dif = 2.*sum(angs>theta & angs<theta+pi) - atom_num;%x cut
%     num_dif = sum(data_windowed(:,3)>-0e-3) - sum(data_windowed(:,3)<-0e-3);
    temp = [temp,num_dif];%num_dif./atom_num

end
sig.cal = temp(data.sync.cal.mask);
sig.raw = temp(data.sync.msr.mask);
end