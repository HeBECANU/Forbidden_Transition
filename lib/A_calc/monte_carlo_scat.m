%monte carlo sim of particle leaving the thermal cloud
lambda_706 = 706.53e-9;
lambda_1083 = 1083.3e-9;

E_427 = 1/(2*mhe)*(h*v0/c)^2;
v_427 = 2*sqrt(E_427/mhe);
E_706 = 1/(2*mhe)*(h/lambda_706)^2;
v_706 = 2*sqrt(E_706/mhe);
E_1083 = 1/(2*mhe)*(h/lambda_1083)^2;
v_1083 = 2*sqrt(E_1083/mhe);

v_avg = 4*sqrt(kb*nanmean(T)/(pi*mhe));

no_trials = 2000;%number of trials
detection_efficeny = 0.01;
N = 5e5/detection_efficeny;

sig = pi*const.ahe_scat^2;

x= R_x.*randn(no_trials,1);
y= R_y.*randn(no_trials,1);
z= R_z.*randn(no_trials,1);

r = [x,y,z];

phi = rand(no_trials,2).*pi;
theta = rand(no_trials,2).*2*pi;

k_x=ones(no_trials,1).*v_427+cos(theta(:,1)).*sin(phi(:,1)).*v_706+...
    +cos(theta(:,2)).*sin(phi(:,2)).*v_1083;
k_y=sin(theta(:,1)).*sin(phi(:,1)).*v_706+...
    +sin(theta(:,2)).*sin(phi(:,2)).*v_1083;
k_z=cos(phi(:,1)).*v_706+...
    cos(phi(:,2)).*v_1083;

k = [k_x,k_y,k_z];

s_rate = @(x) v_avg*sig*N/(pi^(3/2)*R_x*R_y*R_z).*exp(-(x(:,1)./R_x).^2)...
    .*exp(-(x(:,2)./R_y).^2).*exp(-(x(:,3)./R_z).^2);
s = zeros(no_trials,1);

t_step = 0.01;
step_num = 10000;
for ii = 1:step_num
    s = s+s_rate(r);
    r = r+k.*t_step;
end

sum((s.*t_step)>1)/no_trials
mean((s.*t_step))
mean((s./no_trials))
