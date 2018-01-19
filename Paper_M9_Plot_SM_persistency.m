%{ 
Comparison of 3month SW, LW, RH climatology between Maurer vs 4 SDs
         plot_type   temp_res.     climatology  spatial res.   variables
 -------------------------------------------------------------------------------
 Fig 1-1 map         3month         yes         grid           sw
 Fig 2-1 map         3month         yes         grid           lw
 Fig 3-1 map         3month         yes         grid           rh 

This script does not use CLM history file. Instead, use forcing dataset
(netCDF prepared in preprocess)
%}

close all
clear
clc

tic
%% Definition
% -------------------------------------------------------------------------
% Forcing product list (in orde of increasing resolution)
% ------------------------------------------------------------------------- 
prdct{1,1}='BCCA12K';
prdct{2,1}='BCSD12K';
prdct{3,1}='BCSDdisag12K';
prdct{4,1}='BCSAR12K';
prdct{5,1}='MAURER12K';
% -------------------------------------------------------------------------
% Grid configuration option 
% -------------------------------------------------------------------------
region{1,1} = 'NE';
region{2,1} = 'MA';
region{3,1} = 'SA';
region{4,1} = 'GL';
region{5,1} = 'OH';
region{6,1} = 'TN';
region{7,1} = 'UM';
region{8,1} = 'LM';
region{9,1} = 'RR';
region{10,1} = 'MR';
region{11,1} = 'AR';
region{12,1} = 'GUL';
region{13,1} = 'RIO';
region{14,1} = 'UCO';
region{15,1} = 'LCO';
region{16,1} = 'GB';
region{17,1} = 'PN';
region{18,1} = 'CA';
% -------------------------------------------------------------------------
% Model
% -------------------------------------------------------------------------
model{1} = 'CLM';
model{2} = 'VIC';
model{3} = 'PRMS';
% -------------------------------------------------------------------------
% resolution
% -------------------------------------------------------------------------
res = 12;
% -------------------------------------------------------------------------
% Starting and ending time period for analysis
% -------------------------------------------------------------------------
byr=1980; bmon=10; bday=1;
%byr=2001; bmon=10; bday=1;
%eyr=2008; emon=9; eday=30;
eyr=1999; emon=9; eday=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Variables that you don't need to change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Read conus domain grid 
% -------------------------------------------------------------------------
us_domain_grid = ['/d3/mizukami/domain_huc/domain_conus_' num2str(res) 'k.nc'];
lon2d = netcdf2mat(us_domain_grid,'x'); lon2d(lon2d == -999) = NaN;
lat2d = netcdf2mat(us_domain_grid,'y');
lon2d = 360+lon2d;
lon1d = netcdf2mat(us_domain_grid,'lon');
lat1d = netcdf2mat(us_domain_grid,'lat');
ele   = netcdf2mat(us_domain_grid,'ele');
huc2g = netcdf2mat(us_domain_grid,'huc2'); huc2g(huc2g == -999) = NaN;
xlon_num = length(lon1d);
ylat_num = length(lat1d);
clear us_domain_grid
% -------------------------------------------------------------------------
% Read shapefile for western region HUC2
% -------------------------------------------------------------------------
us_huc2_shp = '/d3/mizukami/domain_huc/shapefile/HUC02_conus.shp';
S = shaperead(us_huc2_shp);
clear us_huc2_shp
% -------------------------------------------------------------------------
% water year array 
% -------------------------------------------------------------------------
wyr = byr+1:1:eyr;

%Import color map for surface map        
mycolors  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/rainbow.gp',3,10);
mycolors  = (1/255).*mycolors;
mycolors1  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/BlRe.rgb',3,2);
mycolors1  = (1/255).*mycolors1;
mycolors2  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/temp_diff_18lev.rgb',3,6);
mycolors2  = (1/255).*mycolors2;
mycolors3  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/MPL_Blues.rgb',3,2);
% line color order
clr4=[0 1 1;1 0 0;0 1 0;0 0 1];
clr5=[0 0 1;0.65 0.16 0.16;1 0 0;0 1 0;0 0 0];
%% Getting daily clm datm variable from various product files
tlag = 0:5:180;
% -------------------------------------------------------------------------
% Memory allocation
% -------------------------------------------------------------------------
% CONUS [lat x lon x tlag x prod]          
rho_us     = zeros(ylat_num,xlon_num,length(tlag),length(prdct),length(model));
%Create simulation area mask
sim_mask1 = repmat(lon2d,[1,1,length(tlag),length(prdct),length(model)]);
sim_mask1(~isnan(sim_mask1)) = 1;

sim_mask2 = repmat(lon2d,[1,1,length(prdct),length(model)]);
sim_mask2(~isnan(sim_mask2)) = 1;

for r = 1:length(region)
  % -------------------------------------------------------------------------
  % Grid configuration option 
  % -------------------------------------------------------------------------
  if strcmp(region{r},'NE') 
      huc=1;
  elseif strcmp(region{r},'MA') 
      huc=2;
  elseif strcmp(region{r},'SA') 
      huc=3;
  elseif strcmp(region{r},'GL') 
      huc=4;
  elseif strcmp(region{r},'OH') 
      huc=5;
  elseif strcmp(region{r},'TN') 
      huc=6;
  elseif strcmp(region{r},'UM') 
      huc=7;
  elseif strcmp(region{r},'LM') 
      huc=8;
  elseif strcmp(region{r},'RR') 
      huc=9;
  elseif strcmp(region{r},'MR') 
      huc=10;
  elseif strcmp(region{r},'AR')
      huc=11;
  elseif strcmp(region{r},'GUL')
      huc=12;
  elseif strcmp(region{r},'RIO')
      huc=13;    
  elseif strcmp(region{r},'UCO')
      huc=14;
  elseif strcmp(region{r},'LCO')
      huc=15;
  elseif strcmp(region{r},'GB')
      huc=16;    
  elseif strcmp(region{r},'PN')
      huc=17;         
  elseif strcmp(region{r},'CA')
      huc=18;     
  end
  % -------------------------------------------------------------------------
  % Read region domain grid 
  % -------------------------------------------------------------------------
  region_domain_grid = ['/d3/mizukami/domain_huc/domain_' region{r} '_' num2str(res) 'k.nc'];
  lon_region_1d = netcdf2mat(region_domain_grid,'lon');
  lat_region_1d = netcdf2mat(region_domain_grid,'lat');
  xlon_region_num = length(lon_region_1d);
  ylat_region_num = length(lat_region_1d);
  clear region_domain_grid
  
  % get indices for region domain from western domain array
  i1 = find( lat1d==lat_region_1d(1) );
  i2 = find( lat1d==lat_region_1d(end) );               
  j1 = find( lon1d==lon_region_1d(1) );
  j2 = find( lon1d==lon_region_1d(end) );
  % -------------------------------------------------------------------------
  %  Reading NetCDF
  % -------------------------------------------------------------------------
  for p=1:length(prdct)
      for m = 1:length(model)
          % -------------------------------------------------------------------------
          % Directory where forcing data is located
          % -------------------------------------------------------------------------
          if strcmp(model{m},'CLM')
              main_indir ='/d3/mizukami/CLM_OUTPUT';
          elseif strcmp (model{m},'VIC')
              main_indir ='/d3/mizukami/VIC_OUTPUT/netcdf/processed';
          elseif strcmp (model{m},'PRMS')
              main_indir ='/d3/mizukami/PRMS_OUTPUT';
          end
          % grid box
          ncname = [main_indir '/' prdct{p} '_'  region{r} '/calib/SM_autocorr.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'rho','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          rho_us(i1:i2,j1:j2,:,p,m) = var0+rho_us(i1:i2,j1:j2,:,p,m);
      end
  end
end
rho_us      = rho_us.*sim_mask1;
%%
edecayt_us  = zeros(ylat_num,xlon_num,length(prdct),length(model));
r1=ones(ylat_num,xlon_num,length(prdct),length(model));
for lag=1:length(tlag)
    r1( squeeze( rho_us(:,:,lag,:,:) ) < exp(-1) | r1==0 )=0;
    edecayt_us = edecayt_us+r1;
    
end
edecayt_us=edecayt_us*5;
% -------------------------------------------------------------------------
% Make plot 
% -------------------------------------------------------------------------
%% plot 1 - 1 month auto correlation
figure('Units','inches','position',[1 4 6.00 6.75],'Visible','on','Color',[1 1 1]);
%M02
for i = 1:3
  subplot(5,3,i)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  pcolorm(lat2d, lon2d, edecayt_us(:,:,5,i));
  caxis([0 180]);
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  colormap(mycolors);freezeColors;
  framem off
end
%BCCA-M02
for i = 1:3
  subplot(5,3,i+3)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  pcolorm(lat2d, lon2d, edecayt_us(:,:,1,i)-edecayt_us(:,:,5,i));
  caxis([-20 20]);
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem off
end
%BCSDd-M02
for i = 1:3
  subplot(5,3,i+6)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  framem off
  pcolorm(lat2d, lon2d, edecayt_us(:,:,2,i)-edecayt_us(:,:,5,i));
  caxis([-20 20]);

  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem off
end
%BCSDm-M02
for i = 1:3
  subplot(5,3,i+9)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  pcolorm(lat2d, lon2d, edecayt_us(:,:,3,i)-edecayt_us(:,:,5,i));
  caxis([-20 20]);
  plotm([S.Y],[S.X], 'k');
  tightmap; colormap(mycolors2); freezeColors;
  framem off
end
%AR-MA02
for i = 1:3
  subplot(5,3,i+12)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  pcolorm(lat2d, lon2d, edecayt_us(:,:,4,i)-edecayt_us(:,:,5,i));
  caxis([-20 20]);
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove = 0.2;
hscale=1.2;
vscale=1.2;
hspace=0.10;
vspace=0.02;
%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)-hspace*hmove  pos(1,2)-vspace*2.2  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)-hspace*hmove  pos(2,2)-vspace*2.2  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)-hspace*hmove  pos(3,2)-vspace*2.2  pos(3,3)*hscale pos(3,4)*vscale])

set(h(4), 'Position',[ pos(4,1)-hspace*hmove  pos(4,2)-vspace*2.2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)-hspace*hmove  pos(5,2)-vspace*2.2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)-hspace*hmove  pos(6,2)-vspace*2.2  pos(6,3)*hscale pos(6,4)*vscale])

set(h(7), 'Position',[ pos(7,1)-hspace*hmove  pos(7,2)-vspace*2.2  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)-hspace*hmove  pos(8,2)-vspace*2.2  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)-hspace*hmove  pos(9,2)-vspace*2.2  pos(9,3)*hscale pos(9,4)*vscale])

set(h(10),'Position',[pos(10,1)-hspace*hmove pos(10,2)-vspace*2.2 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)-hspace*hmove pos(11,2)-vspace*2.2 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)-hspace*hmove pos(12,2)-vspace*2.2 pos(12,3)*hscale pos(12,4)*vscale])

set(h(13),'Position',[pos(13,1)-hspace*hmove pos(13,2)+vspace*0.5 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)-hspace*hmove pos(14,2)+vspace*0.5 pos(14,3)*hscale pos(14,4)*vscale])
set(h(15),'Position',[pos(15,1)-hspace*hmove pos(15,2)+vspace*0.5 pos(15,3)*hscale pos(15,4)*vscale])

%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.05,0.880,'M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.675,'BCCA-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.500,'BCSDd-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.325,'BCSDm-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.150,'AR-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.23,0.965,'CLM','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
text(0.50,0.965,'VIC','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
text(0.77,0.965,'PRMS','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);

%CLM  difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-20 20]); colormap(mycolors2);
colorbar('horizontal','position',[0.145,0.030,0.14,0.0225],...
   'XLim',[-20 20],'XTick',-20:20:20,'FontSize',10);
cbfreeze;
%VIC difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-20 20]); colormap(mycolors2);
colorbar('horizontal','position',[0.430,0.030,0.14,0.0225],...
   'XLim',[-20 20],'XTick',-20:20:20,'FontSize',10);
cbfreeze;

%PRMS difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-20 20]); colormap(mycolors2);
colorbar('horizontal','position',[0.725,0.030,0.14,0.0225],...
   'XLim',[-20 20],'XTick',-20:20:20,'FontSize',10);
cbfreeze;

%Maurer PR colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([20 180]); colormap(mycolors);
colorbar('horizontal','position',[0.145,0.775,0.14,0.0225],...
    'XLim',[20 180],'XTick',20:160:180,'FontSize',10);
cbfreeze;

%Maurer RO colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([20 180]); colormap(mycolors);
colorbar('horizontal','position',[0.430,0.775,0.14,0.0225],...
    'XLim',[20 180],'XTick',20:160:180,'FontSize',10);
cbfreeze;

%Maurer ET colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([20 180]); colormap(mycolors);
colorbar('horizontal','position',[0.725,0.775,0.14,0.0225],...
    'XLim',[20 180],'XTick',20:160:180,'FontSize',10);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m9_fig1_edecayt_diff_WYR' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);

%% scatter plot
% hist3 will bin the data
xidecayt = 0:5:190;
yidecayt = xidecayt;
%memory allocation
edecayt_1d = ones(ylat_num*xlon_num,2,length(model)-1,length(prdct))*NaN;
Nedecayt = ones(length(xidecayt),length(yidecayt),length(model)-1,length(prdct))*NaN;

for m = 1:length(model)-1
    for p = 1:length(prdct)
        edecayt_1d(:,:,m,p) = [reshape( edecayt_us(:,:,p,m+1), ylat_num*xlon_num,1) reshape( edecayt_us(:,:,p,1), ylat_num*xlon_num,1)];
        Nedecayt(:,:,m,p)   = hist3(edecayt_1d(:,:,m,p),{xidecayt yidecayt});
    end
end
[Xedecayt,Yedecayt] = meshgrid(xidecayt,yidecayt);

figure('Units','inches','position',[1 4 6.75 3.75],'Visible','on');
%Plots for decayt
for p=1:5
    subplot(2,5,p)
    surface(Xedecayt,Yedecayt,log10(Nedecayt(:,:,1,p)),'EdgeColor','none');caxis([0 3]);
    hold on
    plot3(xidecayt,yidecayt,5*ones(length(xidecayt),1),'k-');
    axis equal
    set(gca,'FontSize',9);
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0 180]);
    set(gca,'XTick',0:180:180);
    set(gca,'YTick',0:180:180);     
%     set(gca,'XTickLabel',{''});
    if p == 1;title('BCCA');
    elseif p == 2;title('BCSDd');
    elseif p == 3;title('BCSDm');    
    elseif p == 4;title('AR');
    elseif p == 5;title('M02');    
    end
end
for p=1:5
    subplot(2,5,p+5)
    surface(Xedecayt,Yedecayt,log10(Nedecayt(:,:,2,p)),'EdgeColor','none');caxis([0 3]);
    hold on
    plot3(xidecayt,yidecayt,5*ones(length(xidecayt),1),'k-');
    axis equal
    set(gca,'FontSize',9);
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0 180]);
    set(gca,'XTick',0:180:180);
    set(gca,'YTick',0:180:180);    
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hspace=0.03;
vspace=0.03;
hscale=1.07;
vscale=1.07;
% %bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hspace*1.65  pos(1,2)+vspace*0.00  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hspace*1.00  pos(2,2)+vspace*0.00  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)-hspace*0.00  pos(3,2)+vspace*0.00  pos(3,3)*hscale pos(3,4)*vscale])
set(h(4), 'Position',[ pos(4,1)-hspace*0.95  pos(4,2)+vspace*0.00  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)-hspace*1.65  pos(5,2)+vspace*0.00  pos(5,3)*hscale pos(5,4)*vscale])

set(h(6), 'Position',[ pos(6,1)+hspace*1.65  pos(6,2)-vspace*1.95  pos(6,3)*hscale pos(6,4)*vscale])
set(h(7), 'Position',[ pos(7,1)+hspace*1.00  pos(7,2)-vspace*1.95  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)-hspace*0.00  pos(8,2)-vspace*1.95  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)-hspace*0.95  pos(9,2)-vspace*1.95  pos(9,3)*hscale pos(9,4)*vscale])
set(h(10), 'Position',[ pos(10,1)-hspace*1.65  pos(10,2)-vspace*1.95 pos(10,3)*hscale pos(10,4)*vscale])

% 
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.525,0.925,'e-fold decay time [day]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0)
text(0.020,0.725,'VIC','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',90);
text(0.525,0.500,'CLM','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
text(0.020,0.275,'PRMS','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',90);
text(0.525,0.075,'CLM','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
text(0.95,0.035,'log_{10}(N)','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
% 
%colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 3]);
colorbar('horizontal','position',[0.80,0.035,0.10,0.025],...
   'XLim',[0 3],'XTick',0:3:3,'FontSize',9);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m9_fig2_edecayt_comparison_model_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);

clear Ndecayt

toc
