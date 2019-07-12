% sim rf knife method for forbidden transition

% find this .m file's path, this must be in the project root dir
this_folder = fileparts(fileparts(which(mfilename)));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));%add all subfolders to the path to find genpath_exclude
path_to_genpath=fileparts(which('genpath_exclude'));
path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(this_folder)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(this_folder,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(this_folder,'dev'),'\.'))
%addpath(genpath_exclude(fullfile(this_folder,'bin'),'\.'))

hebec_constants
global const

%%

in_opts.numsim=1e6;
in_opts.do_plots=0;
in_opts.rf_knife_freq=300e3;

out=k_space_det_eff(in_opts)


%%
outs={};
in_opts=[];
in_opts.numsim=1e5;
in_opts.do_plots=0;
in_opts.rf_knife_freq=300e3;
rf_knife_freqs=col_vec(linspace(10,1000,500)*1e3);
iimax=numel(rf_knife_freqs);
for ii=1:iimax
    in_opts.rf_knife_freq=rf_knife_freqs(ii);
    outs{ii}=k_space_det_eff(in_opts);
end


det_frac_vec=col_vec(cellfun(@(x) x.detected.fraction,outs));

%%
stfig('rf knife freq depedence')
plot(rf_knife_freqs*1e-3,det_frac_vec,'xk')
xlabel('rf freq (khz)')
ylabel('detected fraction of scattered')



%% thermal leakage
%one consideration is how much thermal will leak out from having the rf knife on
% one d Maxwell–Boltzmann
% sqrt(m/(2*pi*k*T))*exp(-m v^2 / (2*k*T))
% which if we integrate fomr vmax to inf and -vmax to -inf gives
% see mathematica document thermal_leak_rate
% Erfc[(m vmax)/(Sqrt[2] Sqrt[k m T])]
% erfc(sqrt(const.mhe)*knife_vel/(sqrt(2*const.kb*1e-6)))

knife_freq=in_opts.rf_knife_freq;
% 1/2 m v^2 = h f
% sqrt(h f 2 /m ) =v_knife
knife_vel=sqrt(const.h*knife_freq*2/const.mhe) ;

leak_rate=@(therm_temp,v_knife) erfc(sqrt(const.mhe)*v_knife/(sqrt(2*const.kb*therm_temp)));
leak_per_atom=leak_rate(1e-6,knife_vel)
total_leak=leak_per_atom*1e6



%% another way the signal could apear is if the rf leakage is changeing based on the temperature
% the derivative of the thermal leak rate with temp. is
leak_rate_dt=@(therm_temp,v_knife) exp(-((const.mhe*v_knife^2)./(2*const.kb*therm_temp)))*const.kb*const.mhe^2*v_knife/...
    (sqrt(2*pi)*(const.kb*const.mhe*therm_temp)^(3/2));

num_atoms=1e6;
delt_temp=5e-9;
starting_temp=1e-6;
iimax=numel(rf_knife_freqs);
therm_leak.increase=[];
therm_leak.mean=[];
for ii=1:iimax
    knife_vel=sqrt(const.h*rf_knife_freqs(ii)*2/const.mhe) ;
    therm_leak.increase(ii)=leak_rate_dt(starting_temp,knife_vel)*num_atoms*delt_temp;
    therm_leak.mean(ii)=leak_rate(starting_temp,knife_vel)*num_atoms;
end

%%

stfig('thermal increase in outcoupling')
subplot(2,1,1)
plot(rf_knife_freqs*1e-3,therm_leak.increase)
xlabel('rf knife freq above bottom (kHz)')
ylabel('atom number increase (arb units)')
subplot(2,1,2)
plot(rf_knife_freqs*1e-3,therm_leak.mean)
xlabel('rf knife freq above bottom (kHz)')
ylabel('atom number (arb units)')

%%




%% now to turn velocity into mhz above trap bottom

% \ub =1.399 625 MHz/G
% \delta E = - \mu · B 
% \mu= -\mu_B * g^{J} /\hbar
% freq_drive=  \mu_B * g^{2s}_J =2.8MHz/G

%E=1/2 m v^2 
% B= E/\mu 
% freq = mu_B * g^{2s}_J*E/\mu 
% freq = mu_B * g^{2s}_J*E *\hbar / \mu_B * g^{J}
% freq = E /h



