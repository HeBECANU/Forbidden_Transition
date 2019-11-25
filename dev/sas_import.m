%read the txt file
clear all
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\ultrawide\';
% fname = '2p_cs_log_20190721T202650.txt';
% fname = '2p_cs_log_20190721T203527.txt';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\width_dep\70MHz\';
% fname = '2p_cs_log_20190721T210118.txt';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\width_dep\10MHz\';
% fname = '2p_cs_log_20190721T205830.txt';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\width_dep\40MHz\';
% fname = '2p_cs_log_20190721T205352.txt';
data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\scan_speed\40pts_40MHz_0.1s\';
fname = '2p_cs_log_20190721T221643.txt';
data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\scan_speed\40pts_40MHz_0.05s\';
fname = '2p_cs_log_20190721T221832.txt';
data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\scan_speed\40pts_40MHz_0.4s\';
fname = '2p_cs_log_20190721T215558.txt';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\power_dep\540uw\';
% fname = '2p_cs_log_20190721T205352.txt';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\power_dep\63uw\';
% fname = '2p_cs_log_20190721T214251.txt';
% data_dir = 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190721_wm_cal\power_dep\730uw\';
% fname = '2p_cs_log_20190721T211406.txt';

%%load the data
path=strcat(data_dir,fname);
fid = fopen(path,'r');
raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
fclose(fid);
lines = split(raw_line,'}}}');
lines = lines(1:end-1);
f = [];
pmt_volt = [];
for ii = 1:(length(lines))
    lines{ii} = [lines{ii},'}}}'];
    sas_dat=jsondecode(lines{ii});
    f = [f,sas_dat.parameters.set_freq'];
    pmt_volt = [pmt_volt,sas_dat.parameters.pmt_voltage_mean'];
    
end
%%
theory_cen = 351721835.04;
% fun1d = @(b,x) b(1).*exp(-((x-b(3)).^2)./(2*b(2).^2))+b(4)+b(5).*x+b(6).*x.^2;
fun1d = @(b,x) b(1)./((x-b(3)).^2+b(2).^2)+b(4)+b(5).*x+b(6).*x.^2;
coeff_names={'amp','sig','mu','offset','lin','quad'};
% fun1d = @(b,x) 129.3994./((x-b(1)).^2+10.9358.^2)+0.3608-0.0047.*x-0.0001.*x.^2;
% coeff_names={'mu'};

fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
% inital_guess=[1.1069800120393503,21.807871898594758,1.3452651113107825,0.37296995919095949,-0.006370908089692269,0];
inital_guess=[129.3994   10.9358   -1.8044    0.3608   -0.0047   -0.0001];
% inital_guess=[-1.8044];
fitobject=fitnlm(f-theory_cen,pmt_volt,...
    fun1d,...
    inital_guess,...
    'CoefficientNames',coeff_names,'Options',fo)
stfig('Saturated Absorption Spectroscopy');
% clf
hold on
x_sample_fit=col_vec(linspace(min(f),max(f),1e3))-theory_cen;
[ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
plot(x_sample_fit,ysamp_val,'b')


scatter(f-theory_cen,pmt_volt,'*')
xlabel('relative freq')
ylabel('pmt voltage')

%% stored data
%scan size parameters
data = [10, 714.2441   11.9910   -0.5511   -3.5091    0.0072    0.0208
    40, 90.7936    9.7445   -1.5828    0.4956   -0.0071   -0.0002
    70, 129.3994   10.9358  -1.8044    0.3608   -0.0047   -0.0001];
%scan size with controlled background
data = [10,-1.337,0.1195
    40,-1.9807,0.089
    70,-1.8296,0.0996];
% ac shift
data = [540,-1.5828,0.074104
    63,1.2396,0.045945
    730,-2.1566,0.0888];
% data bryce took
data_b = [540 8.75 -1.57 0.07
    730 9.80 -2.14 0.08
    390 8.00 -0.67 0.07
    968 11.00 -2.8 0.01
    246 7.3 	-0.18 0.05
    130 6.9 	0.49 0.03
    98 6.5 	0.94 0.05
    63 6.4 	1.22 0.07
    540 8.75 -1.01 0.08
    2090 	15.00 -6.2 0.02];
x=data_b(:,1);
y=data_b(:,3);
w=1./data_b(:,4).^2;
modelFun = @(b,x) b(1)+b(2).*x;%+b(3).*x.^2;
start=[1,1];
wnlm = fitnlm(x,y,modelFun,start,'Weight',w)
xx=linspace(0,2500)';
stfig('bryce power dep data')
clf
hold on
box on
errorbar(data_b(:,1),data_b(:,3)-0.5989,4*log(4).*data_b(:,4),'k.')
plot(xx,predict(wnlm,xx)-0.5989,'k--','linewidth',1.5)
xlabel('Applied laser power ($\mu$W)','interpreter','latex')
ylabel('Center frequency - \(f_r\) (MHz)','interpreter','latex')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',15)
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';