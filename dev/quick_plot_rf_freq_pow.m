%quick plot of rf knife & power dep


rf_knife_dat= [];
rf_knife_dat.freq=col_vec([1.85,1.95,2.05,2.55,2.00,2.45,2.15,2.25,nan]);
%re analysis
rf_knife_dat.res_num=col_vec([,521.49,]);
% rf_knife_dat.res_num=col_vec([824.6,510.2,443,-7,264.6,-55,257,-74,-666]);
rf_knife_dat.detuned_num=col_vec([548,205,196,-141,168.2,-215,119,77,-312]);
rf_knife_dat.diff_num=rf_knife_dat.res_num-rf_knife_dat.detuned_num;
rf_knife_dat.freq(isnan(rf_knife_dat.freq))=4
stfig('rf knife dep')
plot(rf_knife_dat.freq,rf_knife_dat.diff_num,'x')




%%

opt_power_dat=[];
opt_power_dat.power=col_vec([1.1,0.27,0.46,0.75]);
opt_power_dat.res_num=col_vec([443,17.6,-36,141]);
opt_power_dat.detuned_num=col_vec([196,-7,35,111]);
opt_power_dat.diff_num=opt_power_dat.res_num-opt_power_dat.detuned_num;
plot(opt_power_dat.power,opt_power_dat.diff_num,'x')



%% David taking power dep with scans
stfig('probe power dep')

opt_power_dat=[];
opt_power_dat.power=col_vec([1.1,0.2,0.6,1]);
opt_power_dat.res_num=col_vec([190,100,260,280]);
opt_power_dat.res_num_unc=col_vec([20,40,80,70]);
opt_power_dat.detuned_num=col_vec([2.7,4,7,40]);
opt_power_dat.detuned_num_unc=col_vec([4,0,0,30]);
opt_power_dat.diff_num=opt_power_dat.res_num-opt_power_dat.detuned_num;
opt_power_dat.diff_num_unc=sqrt(opt_power_dat.res_num_unc.^2+opt_power_dat.detuned_num_unc.^2);
errorbar(opt_power_dat.power,opt_power_dat.diff_num,opt_power_dat.diff_num_unc,'x')
xlim([0,1.3])
yl=ylim;
ylim([0,yl(2)])