function fitobject = fit_data(data,fun_type)
gauss_fun1d = @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2))+b(4);
loren_fun1d = @(b,x) b(1)./((x-b(2)).^2+b(3).^2) + b(4);
if strcmp(fun_type,'gauss')
    fun1d = gauss_fun1d;
else
    fun1d = loren_fun1d;
end
coeff_names={'amp','mu','sig','offset'};

xdata=data.freq(data.is_shot_good) - data.offset(data.is_shot_good) - data.ac_shift;
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
end
