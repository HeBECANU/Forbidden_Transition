%oscillator strength via the hating method

hebec_constants
kb = const.kb;
mhe = const.mhe;
h = const.h;
c = const.c;
alpha = 3*8*pi*h*kb/c^2;
%trap freq
wx=2*pi*40;
wy=2*pi*420;
wz=2*pi*40;
predicted_freq=700939267; %MHz


%import heating data
load('Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190716_forbidden427_overnight_heating_method\out\20190718T110943\data_results.mat')
Ti = out_data.data.signal.msr.val(:,2);
Tf = out_data.data.signal.msr.val(:,3);
T=(Ti+Tf)./2;
dT_dt = out_data.data.cal.calibrated_signal.val(:);%out_data.data.signal.msr.val(:,1);
f = out_data.data.signal.msr.freq;

Ry= sqrt(2*kb.*T./(mhe*wy^2));
Rz= sqrt(2*kb.*T./(mhe*wz^2));

Ep = 4/3*const.h.*f;%average kinetic energy added by each photon
%convert pd voltage to power

P_ratio = 1- const.h;

sfigure(1001);
plot(Ry.*Rz)
ylabel('Ry*Rz')

%%
convert = 1;%0.9044.*3.2616e-2;
gauss_fun1d = @(b,x) b(1).*exp(-((x-b(3)).^2)./(2*b(2).^2));
loren_fun1d = @(b,x) b(1)./((x-b(3)).^2+b(2).^2)+b(4);
coeff_names={'amp','sig','mu','offset'};

fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
% 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
inital_guess=[4e-4,0.1,4,0];
fitobject=fitnlm(f,dT_dt./(convert),...
    gauss_fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo)
inital_guess_l=[4e-3,10,4];
fitobject_l=fitnlm(f,dT_dt./(convert),...
    loren_fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo)
fit_coeff_l=fitobject_l.Coefficients.Estimate;
sfigure(1002);
clf
fit_coeff=fitobject.Coefficients.Estimate;
fit_se=fitobject.Coefficients.SE;
x_sample_fit=col_vec(linspace(min(f),max(f),1e3));
% [ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on
% plot(x_sample_fit,ysamp_val,'r')
% drawnow
% yl=ylim;
plot(x_sample_fit,ysamp_ci,'color',[1,1,1].*0.5)
[ysamp_val,ysamp_ci]=predict(fitobject_l,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
plot(x_sample_fit,ysamp_val,'b')
plot(x_sample_fit,ysamp_ci,'color',[1,1,1].*0.2)
scatter(f,dT_dt./(convert))
xlabel('freq')
ylabel('signal')
Ep = 3.888*10^(-10);
nanmean(Ry.*Rz)*9.48*10^(-9)*sqrt(pi)*fit_coeff(1)*sqrt(2)*fit_coeff(2)*10^6*1/Ep;
nanmean(Ry.*Rz)*9.48*10^(-9)*pi*fit_coeff_l(1)/fit_coeff_l(2)*10^6*1/Ep;
nanmean(Ry.*Rz)*3.576*10^(-8)*pi*fit_coeff_l(1)/fit_coeff_l(2)*10^6*1/Ep;
Int_1 = (predicted_freq+fit_coeff_l(3))*1e6.*pi*fit_coeff_l(1)/fit_coeff_l(2)