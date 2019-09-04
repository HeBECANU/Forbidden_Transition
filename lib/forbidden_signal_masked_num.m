function output=forbidden_signal_masked_num(data,import_opts)
do_2d_plot=false;

num_shots=numel(data.mcp_tdc.counts_txy);
if do_2d_plot

    not_empty_shots=~cellfun(@isempty,data.mcp_tdc.counts_txy);
    all_txy=cat(1,data.mcp_tdc.counts_txy{not_empty_shots});

    
    all_txy=masktxy_square(all_txy,import_opts.square_mask);

    all_txy=masktxy_2d_circle(all_txy,import_opts.circ_mask);

    dyn_range_pow=0.2;


    num_bins=import_opts.plot.nbins;
    edge_x=linspace(min(import_opts.plot.lim.x),max(import_opts.plot.lim.x),num_bins);
    edge_y=linspace(min(import_opts.plot.lim.y),max(import_opts.plot.lim.y),num_bins);

    spatial_blur=import_opts.plot.blur;
    bin_area=(range(import_opts.plot.lim.x)/num_bins)*(range(import_opts.plot.lim.x)/num_bins);
    [counts,centers]=hist3(all_txy(:,2:3),'edges',{edge_x,edge_y});
    counts=counts/bin_area;
    counts=counts/num_shots;

    %imagesc seems to plot the wrong way round so we transpose here

    if import_opts.plot.cmp_dyn_range
        counts=counts.^dyn_range_pow;
    end
    if  ~spatial_blur==0
        counts=imgaussfilt(counts,spatial_blur);
    end
    stfig('counts during probe')
    imagesc(10^3*centers{1},10^3*centers{2},transpose(counts))
    colormap(viridis())
    set(gca,'Ydir','normal')
    set(gcf,'Color',[1 1 1]);
    title('Spatial Dist. TOP')
    xlabel('X(mm)')
    ylabel('Y(mm)')
    h=colorbar;
    if import_opts.plot.cmp_dyn_range
    xlabel(h,sprintf('Count Density^{%.2f} (m^{-2})',dyn_range_pow))
    else
    xlabel(h,'Count Density (m^{-2})')
    end

end

%% mask each shot

total_num=[];
quadrant_num=[];
for ii=1:num_shots
    txy_shot=data.mcp_tdc.counts_txy{ii};
    if ~isempty(txy_shot)
    txy_shot=masktxy_square(txy_shot,import_opts.square_mask);
    txy_shot=masktxy_2d_circle(txy_shot,import_opts.circ_mask);
    total_num(ii)=numel(txy_shot);
    quadrant_num(ii,:)=count_quadrants(txy_shot);
    else
        warning('empty shot')
    end
end

%% plot as a function of wavelength

predicted_freq=700939267; %MHz


num_diff=[];
num_diff.x=sum(quadrant_num(:,1:2),2)-sum(quadrant_num(:,3:4),2);
num_diff.y=sum(quadrant_num(:,[1,3]),2)-sum(quadrant_num(:,[2,4]),2);
num_diff.quad=sum(quadrant_num(:,[1,4]),2)-sum(quadrant_num(:,[2,3]),2);

is_cal=data.mcp_tdc.probe.calibration;
%is_cal(isnan(is_cal))=1;
is_cal(isnan(is_cal))=0;

is_freq_good=data.wm_log.proc.probe.freq.act.std<5 &...
    (data.wm_log.proc.probe.freq.set-data.wm_log.proc.probe.freq.act.mean)<5 &...
    ~is_cal;

probe_freq=data.wm_log.proc.probe.freq.act.mean*2;%freq in blue
probe_freq=probe_freq-predicted_freq;

signal_unbinned.val=cat(2,col_vec(total_num(is_freq_good)),...
    num_diff.x(is_freq_good),num_diff.y(is_freq_good),num_diff.quad(is_freq_good));
signal_unbinned.names={'tot num','x axis diff','y axis diff','clover diff'};
signal_unbinned.freq=probe_freq(is_freq_good);


signal_bined=[];
probe_freq_bins=col_vec(linspace(-70,70,60));
%probe_freq_bins=col_vec(linspace(-25,20,9));
iimax=numel(probe_freq_bins)-1;
for ii=1:iimax
    signal_bined.freq_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=signal_unbinned.freq<probe_freq_bins(ii+1) & signal_unbinned.freq>probe_freq_bins(ii);
    signal_bined.val(ii,:)=mean(signal_unbinned.val(bin_mask,:),1);
    sum(bin_mask)
    std(signal_unbinned.val(bin_mask))
    signal_bined.unc_val(ii,:)=std(signal_unbinned.val(bin_mask,:),[],1)./sqrt(sum(bin_mask));
    
    signal_bined.freq(ii)=mean(probe_freq_bins(ii:ii+1));
    signal_bined.freq_lims_diff(ii,:)=abs(signal_bined.freq_lims(ii,:)-signal_bined.freq(ii));
end

stfig('counts vs probe freq')
clf
tot_plots=4;
signal_idx=1;
subplot(2,tot_plots,0+signal_idx)
plot(signal_unbinned.freq,signal_unbinned.val(:,signal_idx),'x')
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})

xl=xlim;
subplot(2,tot_plots,tot_plots+signal_idx)
errorbarxy(col_vec(signal_bined.freq),signal_bined.val(:,signal_idx),...
    signal_bined.freq_lims_diff,signal_bined.unc_val(:,signal_idx))
xlim(xl)
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})

signal_idx=2;
subplot(2,tot_plots,0+signal_idx)
plot(signal_unbinned.freq,signal_unbinned.val(:,signal_idx),'x')
xl=xlim;
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})
subplot(2,tot_plots,tot_plots+signal_idx)
errorbarxy(col_vec(signal_bined.freq),signal_bined.val(:,signal_idx),...
    signal_bined.freq_lims_diff,signal_bined.unc_val(:,signal_idx))
xlim(xl)
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})

signal_idx=3;
subplot(2,tot_plots,0+signal_idx)
plot(signal_unbinned.freq,signal_unbinned.val(:,signal_idx),'x')
xl=xlim;
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})
subplot(2,tot_plots,tot_plots+signal_idx)
errorbarxy(col_vec(signal_bined.freq),signal_bined.val(:,signal_idx),...
    signal_bined.freq_lims_diff,signal_bined.unc_val(:,signal_idx))
xlim(xl)
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})

signal_idx=4;
subplot(2,tot_plots,0+signal_idx)
plot(signal_unbinned.freq,signal_unbinned.val(:,signal_idx),'x')
xl=xlim;
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})
subplot(2,tot_plots,tot_plots+signal_idx)
errorbarxy(col_vec(signal_bined.freq),signal_bined.val(:,signal_idx),...
    signal_bined.freq_lims_diff,signal_bined.unc_val(:,signal_idx))
xlim(xl)
xlabel('freq-theory (MHz)')
ylabel('counts')
title(signal_unbinned.names{signal_idx})


%%



output=[];



end
