% sim rf knife method for forbidden transition
fprintf('== Simulating outcoupling ==\n')

% Physical parameters
num_scattered=1e4;
recoil_k1=2*pi/(427); %nm ^-1
recoil_k2=2*pi/(706);
recoil_k3=2*pi/(1083);
k_detector = 0.5*recoil_k1; 

%Scan stuff
k_fracs = linspace(0,2,100);
sample_fracs = zeros(size(k_fracs));
num_tests = length(k_fracs);
% Misc settings
verbose = 0;

for kidx = 1:num_tests
    % set independent variable
    knife_vel=k_fracs(kidx)*recoil_k1;
    
    % Create two-photon scattering distribution in kspace
    k_scatter = recoil_k1*repmat([1,0,0],num_scattered,1);
    seed_decay1 = randn(num_scattered,3);
    k_decay1 = recoil_k2*seed_decay1./vecnorm(seed_decay1,2,2);
    seed_decay2 = randn(num_scattered,3);
    k_decay2 = recoil_k3*seed_decay2./vecnorm(seed_decay2,2,2);
    k_scattered = k_scatter + k_decay1 + k_decay2;
    
    % % Select the atoms that are addressed by the rf knife
    k_scatt_norm=vecnorm(k_scattered,2,2);
    greater_than_knife_mask=k_scatt_norm>knife_vel;
    k_outcoupled=k_scattered(greater_than_knife_mask,:);
    % % remove the velocity that is given up to the trap
    k_outcoupled_norm=vecnorm(k_outcoupled,2,2);
    k_out_unit_vec=k_outcoupled./k_outcoupled_norm;
    k_final=k_outcoupled-k_out_unit_vec*knife_vel;
    k_final_norm = vecnorm(k_final,2,2);
    % Clip by detector field of view
    on_detector_mask = k_final_norm < k_detector;
    sample_frac = sum(on_detector_mask)/num_scattered;
    
    sample_fracs(kidx) = sample_frac;
    
    %Graphical output if desired
    if verbose
        fprintf('Captured %.2f pc of scattered atoms\n',100*sample_frac)
    end
    if verbose > 1
        figure(1)
        clf;
        subplot(3,1,1)
        scatter3(k_scattered(:,1),k_scattered(:,2),k_scattered(:,3),'k.')
        xlabel('x')
        ylabel('y')
        zlabel('z')

        subplot(3,1,2)
        scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')

        subplot(3,1,3)
        scatter3(k_final(:,1),k_final(:,2),k_final(:,3),'k.')
    end
end
figure(2)
plot(k_fracs,sample_fracs,'k.')
hold on
xlabel('k_{RF}/k_{1}')
ylabel('Captured fraction of scattered atoms')
title('Collection efficiency of RF outcoupling 3-photon scatter')