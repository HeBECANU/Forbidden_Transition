% dependance on RF knife point
% using a gaussian fit with no offset
%RF knife height, cen, width, amp, atom num, integrated pd
data.val = [1.75,6.1,1.8,21,1.713928e+04,1.729e+01
            0,3.8,2.1,6.8,1.0558e+05,1.62e+01
            2.0,6.56,2.63,7.3,1.2836e+05,1.917e+01
            2.05,6.0,2.5,7.8,1.3923468e+05,2.0189e+01
            2.45,5.8,2.0,3.4,1.3017e+05,1.7766e+01
            1.85,6.5,2.1,24,9.3908e+04,1.758e+01
            ];
data.unc = [1.75,0.4,0.3,3,2.403e+03,8.502
            0,0.4,0.4,0.11,2.0334e+04,7.97
            2.0,0.12,0.12,0.3,1.169e+04,9.4
            2.05,0.3,0.3,0.8,4.21e+03,1.0
            2.45,0.2,0.2,0.3,9.2596e+03,8.69
            1.85,0.6,0.5,6,6.3478e+03,8.68
            ];
        
counts = [1.75,400,70;
    1.85,1100,400
stfig('RF dep')
subplot(2,1,1)
errorbar(data.val(:,1),data.val(:,4),data.unc(:,4),'kx')
xlabel('RF knife height')
ylabel('Signal')
subplot(2,1,2)
count_unc = data.val(:,4).*data.val(:,5).*data.val(:,6).*3.2616e-5.*sqrt((data.unc(:,4)./data.val(:,4)).^2+(data.unc(:,5)./data.val(:,5)).^2);
errorbar(data.val(:,1),data.val(:,4).*data.val(:,5).*data.val(:,6).*3.2616e-5,count_unc,'kx')
xlabel('RF knife height')
ylabel('\Delta Counts')
stfig('RF cen dep')
errorbar(data.val(:,1),data.val(:,2),data.unc(:,2),'kx')