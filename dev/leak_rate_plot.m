%% thermal leakage
%one consideration is how much thermal will leak out from having the rf knife on
% one d Maxwell–Boltzmann
% sqrt(m/(2*pi*k*T))*exp(-m v^2 / (2*k*T))
% which if we integrate fomr vmax to inf and -vmax to -inf gives
% see mathematica document thermal_leak_rate
% Erfc[(m vmax)/(Sqrt[2] Sqrt[k m T])]
% erfc(sqrt(const.mhe)*knife_vel/(sqrt(2*const.kb*1e-6)))

f23s1_33s1=f2wl(427.7e-9);
leak_rate_fun=@(therm_temp,v_knife) erfc(sqrt(const.mhe)*v_knife/(sqrt(2*const.kb*therm_temp)));
recoil_vel1=freq_to_recoil_vel_he(f23s1_33s1); %absorb 427nm photon

knife_factor=linspace(0.8,2,1e3);
leakage=arrayfun(@(x) leak_rate_fun(1e-6,x*recoil_vel1),knife_factor);

stfig('leak rate')
plot(knife_factor,leakage,'x')
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')

