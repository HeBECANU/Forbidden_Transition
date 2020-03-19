%% Calculate the excited state lifetime from the linewidth of our transition.

clear all
%analysis options
%drift model 'makima' pchip 'spline' 'nearest'?
%plot drift model?
plot_wm_model = 1;
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
% data.offset = zeros(size(data.offset));
%% Apply ac stark shift
% crete the ac stark shift model
wnlm = ac_stark_shift(0);
% generate the actual shift
ac_gradient = wnlm.Coefficients.Estimate(2);
% ac_gradient = 0;
data.ac_shift = ac_gradient.*data.integrated_pd(data.is_shot_good);
ac_err = wnlm.Coefficients.SE(1);

%% Fit model
fitobject = fit_data(data,fun_type);
fit_coeff=fitobject.Coefficients.Estimate;
fit_se=fitobject.Coefficients.SE;

%% Remove broadening effects

%% Calculate state lifetime

%% Write out final anwser

cli_format_text('','c',3)
cli_format_text('FINAL RESULTS','c',3)
cli_format_text('','c',3)
fprintf('\n transition frequnency (direct) %s\n',string_value_with_unc(freq_val,freq_err))
fprintf('\n transition frequnency (heating) %s\n',string_value_with_unc(freq_val_heating,freq_err_heating))
fprintf('\n Constrained frequnency %s\n',string_value_with_unc(other_freq,freq_err))