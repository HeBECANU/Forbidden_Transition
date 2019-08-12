% wm drift model
%format unix time, offset
wm_offset = [];

data_dirs = {};

data.signal = [];
data.freq = [];
data.time = [];
data.atom_num = [];
data.probe_num = [];
data.integrated_pd =[];

for loop_idx=1:length(data_dirs)
    current_dir = data_dirs{loop_idx};
    fprintf('importing data from \n %s \n',current_dir)
    if ~strcmp(current_dir(end),'\')
        current_dir = [current_dir,'\'];
    end
    out_dirs=dir(fullfile(current_dir,'out'));
    out_dirs=out_dirs(3:end);
    out_dirs=out_dirs(cat(1,out_dirs.isdir));
    if size(out_dirs,1)==0
        warning(sprintf('dir \n %s \n does not contain any out dirs',current_dir)), 
    else 
        % convert the folder name (iso time) to posix time
        time_posix=cellfun(@(x) posixtime(datetime(datenum(x,'yyyymmddTHHMMSS'),'ConvertFrom','datenum')),{out_dirs.name});
        [~,sort_idx]=sort(time_posix,'descend');
        out_dirs=out_dirs(sort_idx);
        looking_for_data_dir=1;
        folder_index=1;
        %runs through all the out put dirs for a given run and looks for saved data, if none is there
        %skips that data dir
        while looking_for_data_dir
            try
                out_instance_folder_path=fullfile(current_dir,'out',out_dirs(folder_index).name,'done.txt');
                if (exist(out_instance_folder_path,'file') || ...
                    exist(out_instance_folder_path,'file')) && ...
                    exist(out_instance_folder_path,'file')
                    looking_for_data_dir=0;
                else
                    folder_index=folder_index+1;
                    if folder_index>numel(out_dirs) %if beyon the end of the folder list return nan;
                         looking_for_data_dir=0;
                         folder_index=nan;
                         warning('did not find a valid output direcory for folder %s',current_dir)
                    end
                end
                %~and(isfile([current_dir,'out\',most_recent_dir.name,'\main_data.mat']),isfile([current_dir,'out\',most_recent_dir.name,'\drift_data.mat']))
                %offset = offset + 1;
                %most_recent_dir=out_dirs(end-offset,1);
                %check = drift_data.avg_coef; %check if it has the avg coefs update

            catch e
                fprintf('\n dir: %s didnt work \n',current_dir)
                msgText = getReport(e)
                continue
            end
        end
        if ~isnan(folder_index)

            load(fullfile(current_dir,'out',out_dirs(folder_index).name,'data_results.mat'))
            % now do some serious data plumbing
            %append to main structure
            
            data.signal = cat(1,data.signal,);
            data.freq = cat(1,data.freq,);
            data.time = cat(1,data.time,);
            data.atom_num = cat(1,data.atom_num,);
            data.probe_num = cat(1,data.probe_num,);
            data.integrated_pd =cat(1,data.integrated_pd,);
            if ~isequal(size(drift_data_compiled.to.val),size(drift_data_compiled.wp.qwp))
                error('things are not the same size') 
            end

        end
    end
end
%create drift model
offset = interp1(wm_offset(:,1),wm_offset(:,2),data.time,'spline');

gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4); 
coeff_names={'amp','mu','sig','offset'};

xdata=data.freq - offset;
ydata=data.signal;
num_bins = 15;

probe_freq_bins = linspace(min(xdata),max(xdata),num_bins)

iimax=numel(probe_freq_bins)-1;
signal_bined.freq_std=nan(iimax,1);
signal_bined.val=nan(iimax,1);
signal_bined.unc_val=nan(iimax,1);
signal_bined.freq_mean=nan(iimax,1);
signal_bined.freq_obs_min_max_mean_diff=nan(iimax,2);
for ii=1:iimax
    signal_bined.freq_bin_lims(ii,:)=[probe_freq_bins(ii),probe_freq_bins(ii+1)];
    bin_mask=xdata<=probe_freq_bins(ii+1) & xdata>probe_freq_bins(ii);
     signal_bined.freq_bin_cen(ii)=nanmean(probe_freq_bins(ii:ii+1));
    if sum(bin_mask)==0
        warning('no elements')
        signal_bined.num_bin(ii)=0;
    else
        signal_bined.num_bin(ii)=sum(bin_mask);
        signal_bined.val(ii,:)=nanmean(ydata(bin_mask,:),1);
        signal_bined.unc_val(ii,:)=nanstd(ydata(bin_mask,:),[],1)./sqrt(sum(bin_mask));
        signal_bined.freq_mean(ii)=nanmean(xdata(bin_mask));
        signal_bined.freq_std(ii)=nanstd(xdata(bin_mask));
        signal_bined.freq_obs_min_max(ii,:)=[min(xdata(bin_mask)),max(xdata(bin_mask))];
        signal_bined.freq_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_bin_lims_mean_diff(ii,:)=abs(signal_bined.freq_bin_lims(ii,:)-signal_bined.freq_mean(ii));
        signal_bined.freq_obs_min_max_mean_diff(ii,:)=abs(signal_bined.freq_obs_min_max(ii,:)-signal_bined.freq_mean(ii));
     end
end

plot(signal_bined.freq_mean,signal_bined.val,'o','MarkerSize',5,'MarkerFaceColor',colors_detail(1,:))
    errorbar(signal_bined.freq_mean,signal_bined.val,...
        signal_bined.unc_val(:,1),signal_bined.unc_val(:,1),...
         signal_bined.freq_obs_min_max_mean_diff(:,1), signal_bined.freq_obs_min_max_mean_diff(:,2),...
        'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
         'MarkerFaceColor',colors_detail(3,:),'LineWidth',2.5);

    
    amp_guess=max(ydata);
    ydata_shifted =ydata-min(ydata);
    mu_guess=wmean(xdata,ydata_shifted); %compute the weighted mean
    %sig_guess=sqrt(nansum((xdata-mu_guess).^2.*ydata_shifted)/nansum(ydata_shifted)); %compute the mean square weighted deviation
    sig_guess=10;
    fo = statset('TolFun',10^-6,...
        'TolX',1e-4,...
        'MaxIter',1e4,...
        'UseParallel',1);
    % 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
    inital_guess=[amp_guess,mu_guess,sig_guess,0];
    fitobject=fitnlm(xdata,ydata,...
        gauss_fun1d,...
         inital_guess,...
        'CoefficientNames',coeff_names,'Options',fo);
    fit_coeff=fitobject.Coefficients.Estimate;
    fit_se=fitobject.Coefficients.SE;
    x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
    [ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    hold on
    plot(x_sample_fit,ysamp_val,'r')
    drawnow
    yl=ylim;
    plot(x_sample_fit,ysamp_ci,'color',[1,1,1].*0.5)
    ylim(yl)
    % show the inital guess
    %plot(x_sample_fit,gauss_fun1d(inital_guess,x_sample_fit)*ymultipler)
    fitobject
    amp_str=string_value_with_unc(fitobject.Coefficients.Estimate(1),fitobject.Coefficients.SE(1),'b');
    cen_str=string_value_with_unc(fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2),'b');
    width_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject.Coefficients.SE(3),'b');
    offset_str=string_value_with_unc(abs(fitobject.Coefficients.Estimate(3)),fitobject.Coefficients.SE(3),'b');
    width_units='MHz';
    offset_units='counts';
    amp_units='counts';
    str=sprintf('Gauss fit \n   Cen    %s %s \n   Width %s %s \n   Amp   %s %s \n   Offset %s %s',...
        cen_str,width_units,width_str,width_units,amp_str,amp_units,offset_str,offset_units);
    text(0.01,0.9,str,'Units','normalized'); 
    
    fprintf('transition frequnency %s\n',string_value_with_unc(predicted_freq+fitobject.Coefficients.Estimate(2),fitobject.Coefficients.SE(2)))
                    