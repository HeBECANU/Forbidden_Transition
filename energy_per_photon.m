%calculates the energy added to the system per photon scattered

lambda_706 = 706.53e-9;

E_427 = 1/(2*mhe)*(h*predicted_freq*1e6/c)^2;
E_706 = 1/(2*mhe)*(h/lambda_706)^2;
%sphere convolutions
theta = linspace(0,2*pi,3e1);
phi = linspace(-pi/2,pi/2,3e1);
points = [];
for ii = theta
    for jj = phi
        x1 = E_427/2.*[cos(ii).*sin(jj),sin(ii).*sin(jj),cos(jj)]+[E_427/2,0,0];
        for kk = theta
            for hh = phi
                x2 = E_706/2.*[cos(kk).*sin(hh),sin(kk).*sin(hh),cos(hh)];
                x = x1+x2;
                points = [points;norm(x)];
            end
        end
    end
end
Ep = mean(points)
%energy added by mj = +1 states