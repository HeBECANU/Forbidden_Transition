%calculates the energy added to the system per photon scattered

lambda_706 = 706.53e-9;
lambda_1083 = 1083.3e-9;

E_427 = 1/(2*mhe)*(h*v0/c)^2;
E_706 = 1/(2*mhe)*(h/lambda_706)^2;
E_1083 = 1/(2*mhe)*(h/lambda_1083)^2;

E_p = E_427+E_706+E_1083;
v_p = 2*sqrt(E_p/mhe);

detection_efficeny = 0.01;
N = 2e5/detection_efficeny;

sig = pi*const.ahe_scat^2;
n_avg = N/((R_x*R_y*R_z));
v_avg = 4*sqrt(kb*nanmean(T)/(pi*mhe));
scat_rate = n_avg*v_avg*sig;
avg_time = 4*R_avg/v_p;%time to really leave the cloud
avg_scatters = scat_rate*avg_time;

%monte carlo sim
avg_scatters = 0.55;
s_err = 0.5;

%final value
E_p = 0.2407.*E_p+(1-0.2407).*avg_scatters.*E_p;
E_p_err = (1-0.2407).*avg_scatters.*E_p.*s_err./avg_scatters;
