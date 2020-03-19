%% Calculate the excited state lifetime from the linewidth of our transition.

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
% data.offset = zeros(size(data.offset));
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

%% Remove broadening effects
measured_linewidth = fit_coeff(3)*2;
measured_linewidth_err = fit_se(3)*2;
laser_effective_linewidth = std(wm_offset(:,2))*2*sqrt(2*log(2));
laser_effective_linewidth_err = 0.8;
syms x
eqn = measured_linewidth == 0.54*x+sqrt(0.22*x^2+laser_effective_linewidth^2); %voigt linewidth
true_linewidth = double(solve(eqn));

%% Calculate state lifetime
tau = ceil(1/(2*pi*true_linewidth)*1e3); %life time in nearest nanoseconds
tau_err = tau/2*sqrt((measured_linewidth_err/measured_linewidth).^2+(laser_effective_linewidth_err/laser_effective_linewidth).^2);
%% Write out final anwser

cli_format_text('','c',3)
cli_format_text('FINAL RESULTS','c',3)
cli_format_text('','c',3)
fprintf('\n 3^3S_1 state lifetime (direct) %s ns\n',string_value_with_unc(tau,tau_err))