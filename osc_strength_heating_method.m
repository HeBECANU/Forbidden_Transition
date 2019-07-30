%oscillator strength via the hating method

hebec_constants

%trap freq
wx=;
wy=;
wz=;

%import heating data
load('Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190716_forbidden427_overnight_heating_method\out\20190718T110943\data_results.mat')

dT = out_data.data.signal.msr.val(:,1);
f = out_data.data.signal.msr.freq;

Ep = 4/3*const.h.*f;%average kinetic energy added by each photon
%convert pd voltage to power

P_ratio = 1- const.h
