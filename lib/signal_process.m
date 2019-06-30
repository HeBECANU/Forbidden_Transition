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
temp = [];

xmin = -25e-3;
xmax = 25e-3;
ymin = -30e-3;
ymax = 35e-3;
xybins = 500;
t_window_min = 6.4;
t_window_max = 9.5;

for c_data = data.tdc.counts_txy
    tdc_data = cell2mat(c_data);
    mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
        & tdc_data(:,2)>xmin & tdc_data(:,2)<xmax ...
        & tdc_data(:,3)>ymin & tdc_data(:,3)<ymax ;
    data_windowed = tdc_data(mask,:);
    X_bin=linspace(xmin,xmax,xybins);
    X_bin_width=(xmax-xmin)/xybins;
    [X1d_counts,X1d_edges]=histcounts(data_windowed(:,2),X_bin);
    X1d_centers=mean([X1d_edges(1:end-1);X1d_edges(2:end)]);
    X1d_counts=1e-3*X1d_counts/(X_bin_width);%why 1e-3 factor?

    amp_guess=max(X1d_counts);
    mu_guess=sum(X1d_centers.*X1d_counts)/sum(X1d_counts); %compute the weighted mean
    sig_guess=sum((X1d_centers-mu_guess).^2.*X1d_counts)/sum(X1d_counts); %compute the mean square weighted deviation
    
%     temp = [temp,sig_guess];
    temp = [temp,sig_guess];
    %pause(0.5)
end
sig.cal = temp(data.sync.cal.mask);
sig.raw = temp(data.sync.msr.mask);
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
% temp = [];
% cen = [1.5,4.9].*1e-3;%Centering values, obtained from image of probe beam %1.212,5.253
% r =24.0e-3;
% theta = pi/2;%the angle about which to partiction between 0 and pi from x-axis
% t_window_min = 0.5;
% t_window_max = 6.3;
% 
% for c_data = data.tdc.counts_txy
%     tdc_data = cell2mat(c_data);
%     mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
%         & sqrt((tdc_data(:,2)-cen(1)).^2+(tdc_data(:,3)-cen(2)).^2)<r;
%     data_windowed = tdc_data(mask,:);
%     y_data = data_windowed(:,3)-cen(2);
%     x_data = data_windowed(:,2)-cen(1);
%     angs = cart2pol(x_data,y_data);
%     atom_num = size(data_windowed,1); %atom number in given window
%     num_dif = 2.*sum(angs<theta & angs>theta-pi) - atom_num;%x cut
% %     num_dif = sum(data_windowed(:,3)>-0e-3) - sum(data_windowed(:,3)<-0e-3);
%     temp = [temp,num_dif];%num_dif./atom_num
% 
% end
% sig.cal = temp(data.sync.cal.mask);
% sig.raw = temp(data.sync.msr.mask);
%%
shot_lim = length(data.sync.msr.mask);
temp_msr = [];
temp_cal = [];
temp_no = [];
cen = [1.5,4.9].*1e-3;%4.9
r = 24.0e-3;
t_window_min = 0.5;
t_window_max = 6.3;%6.3;

for  c_data = data.tdc.counts_txy(data.sync.msr.mask(1:shot_lim))
    tdc_data = cell2mat(c_data);
    mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
        & sqrt((tdc_data(:,2)-cen(1)).^2+(tdc_data(:,3)-cen(2)).^2)<r;
    data_windowed = tdc_data(mask,:);
    temp_msr = [temp_msr;data_windowed];
end
for  c_data = data.tdc.counts_txy(data.sync.cal.mask(1:shot_lim))
    tdc_data = cell2mat(c_data);
    mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
        & sqrt((tdc_data(:,2)-cen(1)).^2+(tdc_data(:,3)-cen(2)).^2)<r;
    data_windowed = tdc_data(mask,:);
    temp_cal = [temp_cal;data_windowed];
end
% Bin the data:
pts = linspace(-40e-3, 40e-3, 101);
N_m = histcounts2(temp_msr(:,3), temp_msr(:,2), pts, pts);
N_c = histcounts2(temp_cal(:,3), temp_cal(:,2), pts, pts);
figure(44)
subplot(2,1,1)
imagesc(pts.*1e3, pts.*1e3, N_m)
title('measure avg')
xlabel('x (mm)')
ylabel('y (mm)')
axis equal;
colorbar
set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
subplot(2,1,2)
imagesc(pts.*1e3, pts.*1e3, N_c)
title('calibration avg')
xlabel('x (mm)')
ylabel('y (mm)')
axis equal;
colorbar
set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
% imagesc(pts.*1e3, pts.*1e3, N_m./size(temp_msr,1)-N_c./size(temp_cal,1))
figure(67)
imagesc(pts.*1e3, pts.*1e3, N_m./sum(data.sync.msr.mask)-N_c./sum(data.sync.cal.mask))
title('hist diff')
xlabel('x (mm)')
ylabel('y (mm)')
axis equal;
colorbar
set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
if sum(data.sync.no_atoms.mask)>0
    for  c_data = data.tdc.counts_txy(data.sync.no_atoms.mask(1:shot_lim))
        tdc_data = cell2mat(c_data);
        mask=tdc_data(:,1)>t_window_min & tdc_data(:,1)<t_window_max...
            & sqrt((tdc_data(:,2)-cen(1)).^2+(tdc_data(:,3)-cen(2)).^2)<r;
        data_windowed = tdc_data(mask,:);
        temp_no = [temp_no;data_windowed];
    end
    figure(608)
    N_no = histcounts2(temp_no(:,3), temp_no(:,2), pts, pts);
    imagesc(pts.*1e3, pts.*1e3,N_no)
    title('No atoms avg')
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal;
    colorbar
    set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
    figure(545)
    subplot(3,1,1)
    imagesc(pts.*1e3, pts.*1e3, N_m./sum(data.sync.msr.mask)-N_no./sum(data.sync.no_atoms.mask))
    title('msr minus background')
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal;
    colorbar
    set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
    subplot(3,1,2)
    imagesc(pts.*1e3, pts.*1e3, N_c./sum(data.sync.cal.mask)-N_no./sum(data.sync.no_atoms.mask))
    title('cal minus background')
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal;
    colorbar
    set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
    subplot(3,1,3)
    imagesc(pts.*1e3, pts.*1e3, N_m./sum(data.sync.msr.mask)-N_c./sum(data.sync.cal.mask)-N_no./sum(data.sync.no_atoms.mask))
    title('diff minus background')
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis equal;
    colorbar
    set(gca, 'XLim', pts([1 end]).*1e3, 'YLim', pts([1 end]).*1e3, 'YDir', 'normal');
end
end