function data=import_data(data_dirs)
%extract data from analysis directories
data.signal = [];
data.freq = [];
data.time = [];
data.atom_num = [];
data.probe_num = [];
data.integrated_pd = [];
data.is_shot_good = [];
data.run_time = [];
data.avg_pd = [];

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
                out_instance_folder_path=fullfile(current_dir,'out',out_dirs(folder_index).name,'data_results.mat');
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
            
            data.signal = cat(1,data.signal,out_data.signal.cal.calibrated_signal.val);
            data.freq = cat(1,data.freq,out_data.freq);
            data.time = cat(1,data.time,out_data.time);
            data.atom_num = cat(1,data.atom_num,out_data.atom_num');
            data.probe_num = cat(1,data.probe_num,out_data.probe_num');
            data.integrated_pd = cat(1,data.integrated_pd,out_data.integrated_pd);
            data.is_shot_good = cat(1,data.is_shot_good,out_data.is_shot_good);
            data.run_time = [data.run_time,mean(out_data.time)];
            data.avg_pd = [data.avg_pd,mean(out_data.integrated_pd(out_data.is_shot_good))];
%             if ~isequal(size(drift_data_compiled.to.val),size(drift_data_compiled.wp.qwp))
%                 error('things are not the same size') 
%             end

        end
    end
end
%ensure is shot good is a logical array
data.is_shot_good = logical(data.is_shot_good);
end