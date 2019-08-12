function out_struct=k_space_det_eff(in_opts)

global const

num_scattered=in_opts.numsim;

f23s1_33s1=f2wl(427.7e-9);
%1803 p2 transition
%https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.119.263002
f23s1_23p2=276732186526.2e3;
f23p2_33s1=f23s1_33s1-f23s1_23p2;
recoil_vel1=freq_to_recoil_vel_he(f23s1_33s1); %absorb 427nm photon
recoil_vel2=freq_to_recoil_vel_he(f23p2_33s1); %emit 706nm photon
recoil_vel3=freq_to_recoil_vel_he(f23s1_23p2);%emit 1083nm photon

%knife_vel=recoil_vel1*1;
knife_freq=in_opts.rf_knife_freq;
% 1/2 m v^2 = h f
% sqrt(h f 2 /m ) =v_knife
knife_vel=sqrt(const.h*knife_freq*2/const.mhe) ;
fprintf('rf knife velocity %f mm/s \n',knife_vel*1e3)
det_radius=70e-3/2;

%in the limit of a infinitesimal size detector then you want the recoil velocity at the peak of the radial count density
%vs k norm which is at about 1.2*recoil_vel1
% for our detector the optimal is ~1 recoil

%%


%% make a spont. scatt. halo of radius recoil_vel1 centered at (1,0,0)*recoil_vel1
% this is for a single photon decay
% k_scattered=randn(num_scattered,3);
% k_scattered=k_scattered./repmat(vecnorm(k_scattered,2,2),1,3); 
% k_scattered=k_scattered.*recoil_vel1;
% k_scattered(:,1)=k_scattered(:,1)+recoil_vel1;
% scatter3(k_scattered(:,1),k_scattered(:,2),k_scattered(:,3),'k.')
% xlabel('x')
% ylabel('y')
% zlabel('z')

%% try to simulate a 2 photon decay
%absorbtion photon

num_photons=3;
% photon index,atom index, cartesian index 
k_photons=zeros(num_photons,num_scattered,3);

k_photons(1,:,1)=recoil_vel1; % displace in the x dirn for the first photon
% for the decay photons pick a random direction by generating a norm dist random number in each cartesian axis
rand_dir= randn(num_scattered,3); % random direction
% normalize each to be on the unit sphere
rand_dir=rand_dir./repmat(vecnorm(rand_dir,2,2),1,3); 
% then multiply by the recoil velocity
k_photons(2,:,:)=rand_dir.*recoil_vel2;
% repeat for the 2nd emitted photon
rand_dir= randn(num_scattered,3);
rand_dir=rand_dir./repmat(vecnorm(rand_dir,2,2),1,3); 
k_photons(3,:,:)=rand_dir.*recoil_vel3;

% sum up the k space
k_scattered=squeeze(sum(k_photons,1));


%% find the distribution of the k space radi
dist_from_origin=vecnorm(k_scattered,2,2);


if in_opts.do_plots
    
    %find the mean kenetic energy
    fprintf('mean kinetic energy                %e J , %e * mass of He \n',mean((1/2)*(dist_from_origin.^2)).*[const.mhe,1])
    fprintf('kinetic energy of one 427nm photon %e J , %e * mass of He \n',mean((1/2)*(recoil_vel1.^2)).*[const.mhe,1])
    fprintf('ratio %e \n',mean((dist_from_origin.^2))/(recoil_vel1.^2))
    fprintf('ratio from sum of ke %f \n',mean((dist_from_origin.^2))/(((recoil_vel1^2)+(recoil_vel2^2)+(recoil_vel3^2)))) 

    
    %%
    
    out=smooth_hist(dist_from_origin,'sigma',1e-4,'lims',[0,sum([recoil_vel1,recoil_vel2,recoil_vel3])]);
    %% plot the k_space
    stfig('k space dist. inital');
    subplot(2,2,1)
    scatter3(k_scattered(:,1),k_scattered(:,2),k_scattered(:,3),'k.')
    xlabel('x')
    ylabel('y')
    zlabel('z')



    stfig('k space dist. inital');
    subplot(2,2,2)
    plot(out.bin.centers,out.count_rate.smooth)
    xlabel('k norm')
    ylabel('radial count density')
    hold on
    yl=ylim;
    line([1,1]*knife_vel,yl,'color','k')
    line([1,1]*recoil_vel1*2,yl,'color','r')
    hold off
end
% should normalize by the shel size at each radial bin to give the count density
%% cdf


if in_opts.do_plots
    cum_density=cumsum(out.counts.raw)./num_scattered;
    stfig('k space dist. inital');
    subplot(2,2,3)
    plot(out.bin.centers,cum_density)
    xlabel('radial distance from origin')
    ylabel('cumulative count fraction')
    hold on
    line([1,1]*knife_vel,[0,1.1],'color','k')
    line([1,1]*recoil_vel1*2,[0,1.1],'color','r')
    ylim([0,1.1])
    hold off
end


%% try to find the radial distribution about the center
meank=mean(k_scattered);
dist_from_mean=vecnorm(k_scattered-repmat(meank,size(k_scattered,1),1),2,2);
if in_opts.do_plots
    
    out=smooth_hist(dist_from_mean,'sigma',1e-5,'lims',[min(dist_from_mean),max(dist_from_mean)]);
    stfig('k space dist. inital');
    subplot(2,2,4)
    plot(out.bin.centers,out.count_rate.smooth)
    xlabel('distance from center')
    ylabel('radial count density')
    hold on
    yl=ylim;
    line([1,1]*(recoil_vel2-recoil_vel3),yl,'color','k')
    line([1,1]*recoil_vel1,yl,'color','r')
    hold off
end

% should normalize by the shel size at each radial bin to give the count density

%% use the rf kife to lense in momentum space
k_scatt_norm=vecnorm(k_scattered,2,2);
greater_than_knife_mask=k_scatt_norm>knife_vel;
num_outcoupled=sum(greater_than_knife_mask);
if num_outcoupled==0
    warning('no atoms outcoupled')
end

out_struct.outcoupled.fraction=num_outcoupled/num_scattered;
out_struct.outcoupled.num=num_outcoupled;

fprintf('outcoupled fraction %f\n',num_outcoupled/num_scattered)

k_outcoupled=k_scattered(greater_than_knife_mask,:);

if in_opts.do_plots
    stfig('k space dist. after knife');
    subplot(2,2,1)
    scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')
    title('sudo k space')
end

%% remove the velocity that is given up to the trap
k_outcoupled_norm=vecnorm(k_outcoupled,2,2);
k_out_unit_vec=k_outcoupled./repmat(k_outcoupled_norm,1,3);
k_outcoupled=k_outcoupled-k_out_unit_vec*knife_vel;
if in_opts.do_plots
    stfig('k space dist. after knife');
    subplot(2,2,2)
    scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')
    title('out k space')
end

%% find the distirbution on the detector
%z(t)=x(0)+z'(0)*t-1/2 g t^2
%0=fall_dist+z'(0)*t-1/2 g t^2
% use the quadratic formula
%t=(z'(0)+sqrt((z'(0))^2+2*g*fall_dist))/g
% do it properly later
detector_fall_distance=0.8517; % check this value
fall_time=(k_outcoupled(:,3)+sqrt((k_outcoupled(:,3)).^2+2*const.g0*detector_fall_distance))/const.g0;
xy_det=k_outcoupled(:,1:2).*repmat(fall_time,1,2);

%%
radial_det_distance=vecnorm(xy_det,2,2);
detected_mask=radial_det_distance<det_radius;
num_det_counts=sum(detected_mask);
out_struct.detected.number=num_det_counts;
out_struct.detected.fraction=num_det_counts/num_scattered;
fprintf('detected count fraction %f\n',out_struct.detected.fraction)
if in_opts.do_plots
    %% plot the dectected and undetected counts
    stfig('detector dist. after knife');
    subplot(2,2,1)
    scatter(xy_det(detected_mask,1),xy_det(detected_mask,2),'.k')
    hold on
    scatter(xy_det(~detected_mask,1),xy_det(~detected_mask,2),'.b')
    hold off


    %%
    dyn_range_pow=0.5;
    spatial_blur=1;
    % x_edges=col_vec(linspace(min(min(xy_det(:,1))),max(xy_det(:,1)),1e2));
    % y_edges=col_vec(linspace(min(min(xy_det(:,2))),max(xy_det(:,2)),1e2));
    x_edges=col_vec(linspace(-det_radius,det_radius,1e2));
    y_edges=col_vec(linspace(-det_radius,det_radius,1e2));
    bin_area=diff(x_edges(1:2))*diff(y_edges(1:2));
    [counts,centers]=hist3(xy_det,'edges',{x_edges,y_edges});
    counts=counts/bin_area;

    %imagesc seems to plot the wrong way round so we transpose here

    if dyn_range_pow~=1
        counts=counts.^dyn_range_pow;
    end
    if  ~spatial_blur==0
        counts=imgaussfilt(counts,spatial_blur);
    end

    stfig('detector dist. after knife');
    subplot(2,2,2)
    imagesc(10^3*centers{1},10^3*centers{2},transpose(counts))
    colormap(viridis)
    set(gca,'Ydir','normal')
    set(gcf,'Color',[1 1 1]);
    title('Spatial Dist. TOP')
    xlabel('X(mm)')
    ylabel('Y(mm)')
    c = colorbar;
    c.Label.String=sprintf('count density^{%.2f}',dyn_range_pow);
end
        
%% find the distribution of the detector radi
if in_opts.do_plots
    
    out=smooth_hist(radial_det_distance,'sigma',1e-4,'lims',[0,max(radial_det_distance)]);
    
    stfig('detector dist. after knife');
    subplot(2,2,3)
    plot(out.bin.centers*1e3,out.count_rate.smooth)
    xlabel('radial distance cen of det (mm)')
    ylabel('radial count density')
    hold on
    yl=ylim;
    line([1,1]*det_radius*1e3,yl,'color','k')
    hold off
    % should normalize by the shel size at each radial bin to give the count density
    %% cdf
    cum_density=cumsum(out.counts.smooth)./num_scattered;
    stfig('detector dist. after knife');
    subplot(2,2,4)
    plot(out.bin.centers*1e3,cum_density)
    xlabel('radial distance cen of det (mm)')
    ylabel('cumulative count fraction')
    hold on
    yl=ylim;
    line([1,1]*det_radius*1e3,yl,'color','k')
    hold off
end



end