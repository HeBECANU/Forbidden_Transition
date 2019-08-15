% calculates the second (and most complicated) integral required for the A
% values
const = (pi*2*kb.*nanmean(T)/mhe)^(3/2)*1/(wx*wy*wz);
R_r = 2*kb*nanmean(T)/mhe*(1/wx^2+1/wy^2+1/wz^2);
%% import camera data

cam_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190715_camera_measurments';
% Get a list of all files and folders in this folder.
files = dir(cam_dir);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
subDirs = char(subFolders.folder);
Foldername = char(subFolders.name);
I_field = zeros(1080,1920,size(subDirs,1)-2);
pos = zeros(size(subDirs,1)-2,1);
x = (-960:959).*3e-6;
y = (-540:539).*3e-6;
y_indx = 700:830;
x_indx = 
[x_I,y_I] = meshgrid(x,y);
for ii = 3:size(subDirs,1)
    current_folder = [subDirs(ii,:),'\',Foldername(ii,:)];
    %cd(current_folder)
    pos(ii-2) = str2double(Foldername(ii,7:11))*1e-3;
    cam_files = dir(current_folder);
    png_files = char(cam_files.name);
    I_tot = double(zeros(1080,1920));
    for jj =3:(size(png_files,1)-1)
        png_dir = [current_folder,'\',png_files(jj,:)];
        imdata = imread(png_dir);
        I = double(rgb2gray(imdata));
        I_tot = I_tot + I;
    end
    I_tot = I_tot./(size(png_files,1)-3);
    I_tot = I_tot - nanmean(nanmean(I_tot(:,800:900))); % remove background
    I_norm = trapz(y(y_indx),trapz(x(800:1200),I_tot(700:830,800:1200),2));
    I_tot = I_tot./I_norm;
    sfigure(560);
    surf(x(800:1200),y(y_indx),I_tot(y_indx,800:1200))
    shading flat
    I_field(:,:,ii-2) = I_tot;
end
