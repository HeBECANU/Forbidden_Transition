% load in heating data
load('Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190716_forbidden427_overnight_heating_method\out\20190724T095142\data_results.mat')

%set colours
colors_main=[[233,87,0];[33,188,44];[0,165,166]];
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=17;

% bin data points and fit

predicted_freq=700939267; %MHz

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
ylabel_str='';

 errorbar(x,y,...
        dy,...
        'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
         'MarkerFaceColor',colors_detail(3,:),'LineWidth',2.5);
hold on    
xlabel('time of pulse arrival (s)','fontsize',font_size_global,'interpreter','latex')
ylabel(ylabel_str,'fontsize',font_size_global,'interpreter','latex')
set(gca,'fontsize',font_size_global)