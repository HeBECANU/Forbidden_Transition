function wnlm = ac_stark_shift(plot_model)
%measure of ac stark shift

%center of runs for differing average intensity
%weekened
data=[%6.26160218413720e+00,1.96828982727647e-01,1.446229756998192e+01
%rf 2.05
7.65113886773931e+00,0.259,2.018938031264871e+01
%forbidden427_direct_det
7.36535432156470e+00,0.301,2.054670601647565e+01
%narrow
7.48072351365071e+00,0.163,1.954504829263426e+01
%dither
%8.99252361715919e+00,0.1625,1.909103065832209e+01
%rf 2.00
%8.01096300825301e+00,0.0944,1.917097098838119e+01
%miss
7.17312778805293e+00,0.889,1.893204889446240e+01
%bad atom num
6.07657729703271e+00,0.4366,1.579524060634340e+01
];
data_no_RF=[%heating overnight
4.56197793958740e+00,0.6372,1.899150254186126e+01.*22/25
%no RF
5.41580230847063e+00,0.49,1.620427285097443e+01
];
x=data(:,3);
y=data(:,1);
x_RF=data_no_RF(:,3);
y_RF=data_no_RF(:,1);
w=1./data(:,2).^2;
modelFun = @(b,x) b(1)+b(2).*x;
start=[1,1];
wnlm = fitnlm(x,y,modelFun,start,'Weight',w);
xx=linspace(0,21)';
if plot_model
sfigure(8080)
errorbar(x,y,data(:,2),'ko');
hold on
line(xx,predict(wnlm,xx),'color','b')

errorbar(x_RF,y_RF,data_no_RF(:,2),'bs');
%line(xx,predict(wnlm,xx),'color','bs')
hold off
%with ac shift 700939267.98±0.07
%with new ac shift transition frequnency 700939267.68±0.08
%without 700939276.25
end
%%
%no RF 5.41580230847063e+00
%compare to 
[no_RF_predict,no_RF_ci]= predict(wnlm,1.620427285097443e+01);
RF_shift = no_RF_predict-5.41580230847063e+00;
RF_err = sqrt(range(no_RF_ci).^2+0.5.^2);
end