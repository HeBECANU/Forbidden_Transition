t_period = 30*3600;%seconds ie one hour
n_period = floor(t_period*length(t_temp)/range(t_temp));
offset_std = [];
for ii = 1:(length(offset_mdl)-n_period)
    offset_std = [offset_std,std(offset(ii:(ii+n_period)))];
end
sfigure(45)
plot(t_temp(1:end-n_period),offset_std)
mean(offset_std)
std(offset_std)