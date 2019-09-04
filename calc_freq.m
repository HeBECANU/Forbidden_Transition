%calculates frequency center of transition (including all shifts)

%analysis options
%drift model 'makima' pchip 'spline' 'nearest'?
%plot drift model?
plot_wm_model = 0;
%which model to use

% wm drift model
%format: unix time, offset
wm_offset = [
    1.563278643754e9,-0.71111273765563965
    1.563183199886e9,-1.1711589097976685
    1.563115846138e9,-1.7073591947555542
    1.563016150706e9,-0.57892698049545288
    1.562861326011e9,-3.2293498516082764
    1.562724583967e9,-1.4483662247657776
    1.5637490316e9,-1.85798978805542
];

%data directories we wish to use for the final graph and value
data_dirs = {'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190713_forbidden427_direct_det_narrow\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190710_forbidden427_direct_det_narrow_dither_on\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190710_forbidden427_direct_det\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden_Rf_2.00\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_weekend_run\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden_Rf_2.05\'
    'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190715_forbidden_427_narrow_scan_missed\'
    };

%% import the preanalysied data from the give directories
cli_format_text('','c',3)
cli_format_text('IMPORTING DATA','c',3)
cli_format_text('','c',3)
data = import_data(data_dirs);

%% create drift model
offset = 2.*interp1(wm_offset(:,1),wm_offset(:,2),data.time,'pchip');
t_temp = linspace(min(wm_offset(:,1)),max(wm_offset(:,1)),4000);
offset_mdl =  2.*interp1(wm_offset(:,1),wm_offset(:,2),t_temp,'pchip');% 
%plot drift model
if plot_wm_model
stfig('wm drift model')
clf
hold on
plot(t_temp,offset_mdl)
scatter(wm_offset(:,1),2.*wm_offset(:,2),'kx')
ylabel('wavemetre offset (MHz)')
xlabel('Posix time (s)')
box on
end

%% Apply ac stark shift
% crete the ac stark shift model
wnlm = ac_stark_shift(0);
% generate the actual shift
ac_gradient = wnlm.Coefficients.Estimate(2);
ac_shift = ac_gradient.*data.integrated_pd(data.is_shot_good);%2.95062819924799e-01.

%% Fit model
gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4);
loren_fun1d = @(b,x) b(1)./((x-b(2)).^2+b(3).^2) + b(4);
fun1d = loren_fun1d;
coeff_names={'amp','mu','sig','offset'};

xdata=data.freq(data.is_shot_good) - offset(data.is_shot_good) - ac_shift;
ydata=data.signal(data.is_shot_good);

%cut outliers
pd_tol = 5;
sigma_from_mean=(ydata-nanmean(ydata,1))/nanstd(ydata,1);
is_not_oulier=sigma_from_mean<7 & data.integrated_pd(data.is_shot_good)>pd_tol;
ydata=ydata(is_not_oulier);
xdata=xdata(is_not_oulier);

%intial guess
amp_guess=4.556661697251370e+01;
mu_guess=1;
sig_guess=2.866375907039084e+00;
offset_guess=-3.340314780403577e-01;
fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
inital_guess=[amp_guess,mu_guess,sig_guess,offset_guess];
fitobject=fitnlm(xdata,ydata,...
    fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo);
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
% RF stark shift
[no_RF_predict,no_RF_ci]= predict(wnlm,1.620427285097443e+01);
RF_shift = no_RF_predict-5.41580230847063e+00;
RF_err = sqrt(range(no_RF_ci).^2+0.5.^2);
total_shift = recoil_shift+Zeeman_shift+RF_shift;
% Calculate final value
freq_val = predicted_freq+fit_coeff(2)-total_shift;
freq_err = (fit_se(2)^2+RF_err^2+Zeeman_err^2)^(1/2);
%% Write out final anwser

cli_format_text('','c',3)
cli_format_text('FINAL RESULTS','c',3)
cli_format_text('','c',3)
fprintf('\n transition frequnency %s\n',string_value_with_unc(freq_val,freq_err))
