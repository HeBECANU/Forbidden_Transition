clear all
anal_opts=[]; %reset the options (would be good to clear all variables except the loop config
anal_opts.tdc_import.dir='Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden427_overnight_heating_method\individual_shots\';
% anal_opts.tdc_import.dir= 'Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden_Rf_2.25\';
font_size_global = 20;
anal_opts.tdc_import.save_cache_in_data_dir=true;
tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
tlim=[0,inf];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];
anal_opts.dld_aquire=22;
anal_opts.trig_dld=20.3;
anal_opts.probe.t0=20;
anal_opts.probe.duration=6;


cli_format_text('','c',2)
cli_format_text('STARTING ANALYSIS','c',2)
cli_format_text('','c',2)

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;
load('Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden427_overnight_heating_method\out\20190724T095142\data_results.mat')
data.mcp_tdc = import_mcp_tdc_data(anal_opts.tdc_import);
%%
tlims = [0, 24.5];
t_sig = 0.001;
num_bins = 100000;
counts_txy = data.mcp_tdc.counts_txy{1};
[z t] = histcounts(counts_txy(:,1),num_bins);
stfig('time of flight for heating')
clf
t = (t(1:end-1)+t(2:end))./2;
data_gauss = gaussfilt(t,z,t_sig);
% data_gauss = z;
plot(t,data_gauss./anal_opts.global.qe.*96/288,'k-')
xlim(tlims)
xlabel('Arrivial time (s)','fontsize',font_size_global,'interpreter','latex')
ylabel('Count rate (kHz)','fontsize',font_size_global,'interpreter','latex')
set(gca,'fontsize',font_size_global)
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
set(gca,'linewidth',1.5)
box on
ylim([1 125])
%%
stfig('Individual Pulse')
clf
plot(t,data_gauss./anal_opts.global.qe,'k-')
xlim([1.735 1.78])
xlabel('Arrivial time (s)','fontsize',font_size_global,'interpreter','latex')
ylabel('Count rate (kHz)','fontsize',font_size_global,'interpreter','latex')
set(gca,'fontsize',font_size_global)
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
set(gca,'linewidth',1.5)
box on
ylim([1 max(data_gauss./anal_opts.global.qe)*1.1])
gauss_fun1d = @(b,x) b(1).*exp(-((x-1.756).^2)./(2*b(3).^2))+b(4); 
b=out_data.data.al_pulses.fit.Estimate(37,1,:);
hold on 
plot(t,gauss_fun1d(b,t))