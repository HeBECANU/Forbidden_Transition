%calculates frequency center of transition (including all shifts)
clear all
%analysis options
%drift model 'makima' pchip 'spline' 'nearest'?
%plot drift model?
plot_wm_model = 0;
fun_type = 'loren';
%which model to use
%% Set up configureations for the direct det analysis
config_direct_det

%% import the preanalysied data from the give directories
cli_format_text('','c',3)
cli_format_text('IMPORTING DATA','c',3)
cli_format_text('','c',3)
data = import_data(data_dirs);

%% create drift model
data.offset = wm_drift_model(data.time,wm_offset,plot_wm_model);

%% Apply ac stark shift
% crete the ac stark shift model
wnlm = ac_stark_shift(0);
% generate the actual shift
ac_gradient = wnlm.Coefficients.Estimate(2);
data.ac_shift = ac_gradient.*data.integrated_pd(data.is_shot_good);
ac_err = wnlm.Coefficients.SE(1);

%% Fit model
fitobject = fit_data(data,fun_type);
fit_coeff=fitobject.Coefficients.Estimate;
fit_se=fitobject.Coefficients.SE;

%% Calculate remaining shifts
% import physical constants
hebec_constants
predicted_freq=700939267; %MHz
% recoil shift
recoil_shift = 1e6*const.h/(2*const.mhe)*(predicted_freq/const.c)^2; %can use predicted due to megahurts unc on terahurts value
% Zeeman shift
Zeeman_shift = -1.7154; %calculated from the B-field
Zeeman_err = 0.003;
% WM unc
wm_err = 4.0;
% Cell shifts
cs_cell_ac_shift = -1.88;
cs_cell_ac_err = 0.43;
% RF stark shift
[no_RF_predict,no_RF_ci]= predict(wnlm,1.620427285097443e+01);
RF_shift = no_RF_predict-5.41580230847063e+00;
RF_err = sqrt(range(no_RF_ci).^2+0.5.^2);
total_shift = recoil_shift+Zeeman_shift+RF_shift+cs_cell_ac_shift
% Calculate final value
freq_val = predicted_freq+fit_coeff(2)-total_shift;
freq_err = (fit_se(2)^2+RF_err^2+Zeeman_err^2+ac_err^2+wm_err^2+cs_cell_ac_err^2)^(1/2)
%% Do shifts for heating freq
freq_heating = 700939271.9;
pd_int_heating = 18.9915;
time_heating = 1563243507.46945;
% Statistical heating error
stat_err = 0.6;
% AC stark shift
ac_shift_heating = ac_gradient.*18.9915;
% WM shift
wm_shift = wm_drift_model(time_heating,wm_offset,plot_wm_model);
total_shift_heating = recoil_shift+Zeeman_shift+cs_cell_ac_shift+ac_shift_heating+wm_shift
% Calculate final value
freq_val_heating = freq_heating-total_shift_heating;
freq_err_heating = (stat_err^2+Zeeman_err^2+ac_err^2+wm_err^2+cs_cell_ac_err^2)^(1/2)
%% Calculate that other frequency
other_freq = freq_val-276736495.6246;
%% Write out final anwser

cli_format_text('','c',3)
cli_format_text('FINAL RESULTS','c',3)
cli_format_text('','c',3)
fprintf('\n transition frequnency (direct) %s\n',string_value_with_unc(freq_val,freq_err))
fprintf('\n transition frequnency (heating) %s\n',string_value_with_unc(freq_val_heating,freq_err_heating))
fprintf('\n Constrained frequnency %s\n',string_value_with_unc(other_freq,freq_err))
