% load in heating data
load('Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\20190716_forbidden427_overnight_heating_method\out\20190724T095142\data_results.mat')
shot = 7;
y = out_data.data.al_pulses.fit.temperature.val(shot,:).*1e6;
dy = out_data.data.al_pulses.fit.temperature.unc(shot,:).*1e6;
x = out_data.data.al_pulses.pos.mean(shot,:,1);

%set colours
colors_main=[[88,113,219];[88,113,219]./1.2;[88,113,219]./1.3]; %[88,113,219]%[96,144,201]
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=17;

% bin data points and fit
fun1d = @(b,x) b(1).*x+b(2);
fo = statset('TolFun',10^-6,...
    'TolX',1e-4,...
    'MaxIter',1e4,...
    'UseParallel',1);
% 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
inital_guess=[1,1];
fitobject=fitnlm(x,y,...
    fun1d,...
    inital_guess,...
    'Options',fo)


colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);

stfig('heating demo')
clf
ylabel_str='Temperature (\(\mu\)K)';
x_sample_fit=col_vec(linspace(-5,29,1e4));
[ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on


curve1 = ysamp_ci(:,1)';
curve2 = ysamp_ci(:,2)';
x1 = (x_sample_fit)';
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, 'g');
h.FaceColor = [0.31 0.31 0.32]*2;
h.FaceAlpha = 0.5;

plot(x_sample_fit,ysamp_val,'k','LineWidth',1.5)
drawnow
yl=ylim;
plot(x_sample_fit,ysamp_ci,'color',[1,1,1].*0.5)

 errorbar(x,y,...
        dy,...
        'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
         'MarkerFaceColor',colors_detail(3,:),'LineWidth',2.5);
hold on    
xlabel('time of pulse arrival (s)','fontsize',font_size_global,'interpreter','latex')
ylabel(ylabel_str,'fontsize',font_size_global,'interpreter','latex')
set(gca,'fontsize',font_size_global)
xlim([0,25])
box on

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height-0.01];


fig = gcf;
set(fig,'Position',[948 545 738 420])
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,'C:\Users\kieran\Documents\MATLAB\Forbidden_Transition\figs\heating_fig','-dpdf')

