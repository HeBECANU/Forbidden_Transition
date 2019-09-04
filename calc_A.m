%A value calculator
clear all
%constants
clear const
hebec_constants
kb = const.kb;
mhe = const.mhe;
h = const.h;
c = const.c;

%trap freq
wx=2*pi*53.5;
dwx=2*pi*0.6;
wy=2*pi*426.56;
dwy=2*pi*0.05;
wz=2*pi*430.3;
dwz=2*pi*0.06;
v0=700939267*1e6; %MHz

%% import heating data
load('Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden427_overnight_heating_method\out\20190718T110943\data_results.mat')
Ti = out_data.data.signal.msr.val(:,2);
Tf = out_data.data.signal.msr.val(:,3);
T=(Ti+Tf)./2;
dT_dt = out_data.data.cal.calibrated_signal.val(:);%out_data.data.signal.msr.val(:,1);
f = out_data.data.signal.msr.freq;
% P = out_data.data.ai_log.pd.mean(~out_data.data.cal.cal_mask).*3.2616e-2; %power in jules
% P(isnan(P)) = 3.2616e-2;
P = ones(1102,1).*3.2616e-2;
%% set up dist variables
R_x = sqrt(2*kb*nanmean(T)/(mhe*wx^2));
dR_x = R_x*sqrt((nanstd(T)/nanmean(T)).^2+(dwx/wx).^2);
R_y = sqrt(2*kb*nanmean(T)/(mhe*wy^2));
dR_y = R_y*sqrt((nanstd(T)/nanmean(T)).^2+(dwy/wy).^2);
R_z = sqrt(2*kb*nanmean(T)/(mhe*wz^2));
dR_z = R_z*sqrt((nanstd(T)/nanmean(T)).^2+(dwz/wz).^2);
R_avg = (R_x*R_y*R_z)^(1/3);
%% proportionality constant
alpha = 3*8*pi*h*kb*v0^2/c^2;

%% Energy per photon

energy_per_photon;

%% Integral 1

first_integral;

%% Integral 2

second_integral;

%% Print out final value
A=alpha/E_p*int_1/int_2;
d_A = A.*sqrt((int_1_err/int_1).^2+(int_2_err/int_2).^2+(E_p_err/E_p).^2);
fprintf(['A_{ul} = (%1.2f',177,'%1.1f)*10^{-9} \n'],A*1e9,d_A*1e9)
       