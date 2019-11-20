%energy level diagram

levels = {'1^1S_0','2^3S_1','2^1S_0','2^3P_2','2^3P_1','2^3P_0','2^1P_1','3^3S_1','2^3P'};
    
energy = [0,19.819,20.61577,20.96408703,20.96409651,20.96421899,21.218023,22.718467];
energy_adj = [18.8,19.819,20.4,20.7,21.12,21.85,22.09,22.718467,21.0433];
state_pos = [1,2,1,3.9,3.9,3.9,1,2,3];
lin_width = 0.25;
y_shift = 0.25;
x_shift = 0.1;
font_size=18.9;
stfig('energy diagram')
clf
for ii = 1:length(state_pos)
    level_str = levels{ii};
    pos = state_pos(ii);
    energy_current = energy_adj(ii);
    plot(pos+[-lin_width,lin_width],energy_current.*[1,1],'k-','LineWidth',2.0)
    if ii<7 && ii>3
        text(pos-x_shift,energy_current-y_shift,['\(',level_str,'\)'],'interpreter','latex','FontSize',font_size-3)
    else
        text(pos-x_shift,energy_current-y_shift,['\(',level_str,'\)'],'interpreter','latex','FontSize',font_size)
    end
    hold on
end
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
%add the arrows
% quiver(state_pos(2),energy_adj(2),state_pos(8)-state_pos(2),energy_adj(8)-energy_adj(2),1,'filled','color',[0.4940 0.1840 0.5560])
% quiver(state_pos(2),energy_adj(2),state_pos(9)-state_pos(2),energy_adj(9)-energy_adj(2),0,'filled')
% quiver(state_pos(9),energy_adj(9),state_pos(8)-state_pos(9),energy_adj(8)-energy_adj(9),0,'filled')
% quiver(state_pos(1),energy_adj(1),state_pos(2)-state_pos(1),energy_adj(2)-energy_adj(1),0,'filled')
plot([state_pos(9)+lin_width,state_pos(4)-lin_width],[energy_adj(9),energy_adj(4)],'k-.')
plot([state_pos(9)+lin_width,state_pos(5)-lin_width],[energy_adj(9),energy_adj(5)],'k-.')
plot([state_pos(9)+lin_width,state_pos(6)-lin_width],[energy_adj(9),energy_adj(6)],'k-.')
box off

set(gca,'Visible','off')
%scatter(state_pos,energy_adj,'kx')