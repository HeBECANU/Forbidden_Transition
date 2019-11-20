% calculates the second (and most complicated) integral required for the A
% values
beta = (pi*2*kb.*nanmean(T)/mhe)^(3/2)*1/(wx*wy*wz);
d_beta = beta.*sqrt((nanstd(T)/nanmean(T)).^2+(dwx/wx).^2+(dwy/wy).^2+(dwz/wz).^2);
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
y_indx = 700:900;
x_indx = 800:1300;
[x_I,y_I] = meshgrid(x(x_indx),y(y_indx));
int_vec = zeros(size(subDirs,1)-2,1);
int_vec_2 = zeros(size(subDirs,1)-2,1);
int_vec_unc = zeros(size(subDirs,1)-2,1);
% loads data in from the images themselves
% for ii = 3:size(subDirs,1)
%     current_folder = [subDirs(ii,:),'\',Foldername(ii,:)];
%     %cd(current_folder)
%     pos(ii-2) = str2double(Foldername(ii,7:11))*1e-3;
%     cam_files = dir(current_folder);
%     png_files = char(cam_files.name);
%     I_tot = double(zeros(1080,1920));
%     for jj =3:(size(png_files,1)-1)
%         png_dir = [current_folder,'\',png_files(jj,:)];
%         imdata = imread(png_dir);
%         I = double(rgb2gray(imdata));
%         I_tot = I_tot + I;
%     end
%     I_tot = I_tot./(size(png_files,1)-3);
%     I_tot = I_tot - nanmean(nanmean(I_tot(:,800:900))); % remove background
%     I_tot(abs(I_tot)<0.5) = 0;
%     I_norm = trapz(y(y_indx),trapz(x(x_indx),I_tot(y_indx,x_indx),2));
%     I_tot = I_tot./I_norm;
%     sfigure(560);
%     surf(x(x_indx),y(y_indx),I_tot(y_indx,x_indx))
%     shading flat
%     I_field(:,:,ii-2) = I_tot;
%     [M,I]= max(I_tot);
%     [m,ind]=max(M);
%     y0 = y(I(ind));
%     z0 = x(ind);
%     prod_val = I_tot(y_indx,x_indx).*exp(-(x_I-z0).^2./R_z.^2).*exp(-(y_I-y0).^2./R_y.^2);%value of the integral over this plane
%     int_val = trapz(y(y_indx),trapz(x(x_indx),prod_val,2));
%     int_vec(ii-2) = int_val;
% end
%%
load('Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190715_camera_measurments\probe_image_data.mat','I_field','pos')
for ii = 1:10
    I_tot = I_field(:,:,ii);
    [M,I]= max(I_tot);
    [m,ind]=max(M);
    y0 = y(I(ind));
    z0 = x(ind);
    prod_val = I_tot(y_indx,x_indx).*exp(-(x_I-z0).^2./R_z.^2).*exp(-(y_I-y0).^2./R_y.^2);%value of the integral over this plane
    prod_val_unc = prod_val.^2.*((d_beta/beta)^2+(x_I-z0).^4.*(2*R_z^(-3)*dR_z)^2+(y_I-y0).^4.*(2*R_y^(-3)*dR_y)^2);%value of the integral over this plane
    prod_val_2 = prod_val.^2;
    int_val = trapz(y(y_indx),trapz(x(x_indx),prod_val,2));
    int_val_2 = trapz(y(y_indx),trapz(x(x_indx),prod_val_2,2));
    int_val_unc = trapz(y(y_indx),trapz(x(x_indx),prod_val_unc,2));
    int_vec(ii) = int_val;
    int_vec_2(ii) = int_val_2;
    int_vec_unc(ii) = int_val_unc;
end
x0 = 0.013;
d_x0 = 0.0005;
int_2 = 1/beta.*trapz(pos,int_vec.*exp(-(pos-x0).^2./R_x.^2));
int_2_err = int_2*sqrt((2*R_x^(-1)*dR_x)^2+(2*R_y^(-1)*dR_y)^2+(2*R_z^(-1)*dR_z)^2);%first guess at uncertainty
% int_2_err_int = sqrt(1/beta.^2.*trapz(pos,int_vec_unc.*(exp(-(pos-x0).^2./R_x.^2)).^2+int_vec_2.*(2*R_x^(-3)*dR_x)^2.*(pos-x0).^4.*(exp(-(pos-x0).^2./R_x.^2)).^2))
%sfigure(5556)
%plot(pos,int_vec)%.*const have to multiply final val by const
