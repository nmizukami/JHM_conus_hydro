%MTCLIM cloudy-sky transmittance sensitivity
close all
clear all
clc

DTR = 2:0.1:20;
delta = -1.8:0.3:1.8;
    
B     = 0.031+0.201*exp(-0.185*DTR);
Tfmax = 1-0.9*exp(-B.*DTR.^1.5);
for i=1:length(delta)    
    B_del(:,i)     = 0.031+0.201*exp(-0.185*(DTR+delta(i)));
    Tfmax_del(:,i) = 1-0.9*exp(-B_del(:,i)'.*(DTR+delta(i)).^1.5);
end

%Color setting
mycolors2  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/temp_diff_18lev.rgb',3,6);
mycolors2  = (1/255).*mycolors2;
mycolors2(end-2:end,:)=[];
mycolors2(1:3,:)=[];

%plot
axes
set(gca, 'ColorOrder', mycolors2);
hold all;
for i=1:length(delta)
    plot(DTR,Tfmax_del(:,i),'--','LineWidth',2);
end
for i=1:length(delta)
    plot(DTR,Tfmax_del(:,i)*0.75,'--','LineWidth',2); 
end
plot(DTR,Tfmax,'-k','LineWidth',3);
plot(DTR,Tfmax*0.75,'-k','LineWidth',3);
ylabel('T_{f,max}, [-]');
xlabel('DTR, [\circC]');
set(gca,'XLim',[2 20]);
set(gca,'XTick',2:2:20);
%Title

%Some text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.845,0.625,'\DeltaDTR [\circC]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.625,0.800,'P = 0 [mm]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
text(0.650,0.625,'P > 0 [mm]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
%Colorbar
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-1.92 1.92]); colormap(mycolors2);
colorbar('vertical','position',[0.815,0.125,0.030,0.475],...    
   'YLim',[-1.92 1.92],'YTick',-1.8:0.6:1.8,'FontSize',10);

%save figure
figfile='./figure/paper/Paper_m11_fig1_Cloudy-sky_T_sensitivity';
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-r400',figfile);  
