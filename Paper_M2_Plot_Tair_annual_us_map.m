% Comparison of annual average of 3 daily temperature climatology between M02 and 4 SDs
%          plot_type   temp_res.     climatology  spatial res.   variables
%  -------------------------------------------------------------------------------
%  Fig 1-1 map         annual        yes          grid           Tmin,Tmax,Tavg

tic

close all
clear
clc
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
model = 'CLM';
% model = 'VIC';
% -------------------------------------------------------------------------
% resolution
% -------------------------------------------------------------------------
res = 12;
% -------------------------------------------------------------------------
% Starting and ending time period for analysis
% -------------------------------------------------------------------------
byr=1980; bmon=10; bday=1;
%byr=2001; bmon=10; bday=1;
% eyr=2008; emon=9; eday=30;
eyr=1999; emon=9; eday=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Variables that cannot be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Directory where forcing data is located
% -------------------------------------------------------------------------
if strcmp(model,'CLM')
    main_indir ='/d3/mizukami/CLM_OUTPUT';
elseif strcmp (model,'VIC')
    main_indir ='/d3/mizukami/VIC_OUTPUT/netcdf/processed';
end
% -------------------------------------------------------------------------
% Read western domain grid 
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
us_huc2_shp = '/d3/mizukami/domain_huc/shapefile/HUC/HUC02_conus.shp';
S = shaperead(us_huc2_shp);
clear us_huc2_shp
% -------------------------------------------------------------------------
% water year array 
% -------------------------------------------------------------------------
wyr = byr+1:1:eyr;
%% Getting daily clm datm variable from various product files
% -------------------------------------------------------------------------
% Memory allocation
% -------------------------------------------------------------------------
% CONUS [lat x lon x region x prod]
pr_us   = zeros(ylat_num,xlon_num,length(prdct));
tavg_us = zeros(ylat_num,xlon_num,length(prdct));
tmin_us = zeros(ylat_num,xlon_num,length(prdct));
tmax_us = zeros(ylat_num,xlon_num,length(prdct));
%Create simulation area mask
sim_mask = lon2d+huc2g;
sim_mask(~isnan(sim_mask) ) = 1;
sim_mask = repmat(sim_mask,[1,1,length(prdct)]);

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
      ncname = [main_indir '/' prdct{p} '_'  region{r} '/calib/Annual_climate_PRCP.nc'];
      [var0,fillvalue]=netcdf2mat(ncname,'PRCP','attname1','_FillValue');
      var0(var0==fillvalue)=0;
      pr_us(i1:i2,j1:j2,p) = var0*60*60*24*365+pr_us(i1:i2,j1:j2,p);
      clear var0 ncname  
      
      ncname = [main_indir '/' prdct{p} '_'  region{r} '/calib/Annual_climate_Tair.nc'];
      [var0,fillvalue]=netcdf2mat(ncname,'Tair','attname1','_FillValue');
      var0(var0==fillvalue)=0;
      tavg_us(i1:i2,j1:j2,p) = var0+tavg_us(i1:i2,j1:j2,p);
      clear var0 ncname  
      
      ncname = [main_indir '/' prdct{p} '_'  region{r} '/calib/Annual_clim_daily_Tair.nc'];
      [var0,fillvalue]=netcdf2mat(ncname,'Tmin','attname1','_FillValue');
      var0(var0==fillvalue)=0;
      tmin_us(i1:i2,j1:j2,p) = var0+tmin_us(i1:i2,j1:j2,p);
      clear var0 ncname
      
      ncname = [main_indir '/' prdct{p} '_'  region{r} '/calib/Annual_clim_daily_Tair.nc'];
      [var0,fillvalue]=netcdf2mat(ncname,'Tmax','attname1','_FillValue');
      var0(var0==fillvalue)=0;
      tmax_us(i1:i2,j1:j2,p) = var0+tmax_us(i1:i2,j1:j2,p);
      clear var0 ncname
  end
end

%Mask out the outside the simulation domain
pr_us   = pr_us.*sim_mask;
tavg_us = tavg_us.*sim_mask; tavg_us=tavg_us-272.15;
tmin_us = tmin_us.*sim_mask;
tmax_us = tmax_us.*sim_mask;

%% Compute stat
fprintf('\nM02 pr max = %f\n',nanmax(reshape(pr_us(:,:,5),ylat_num*xlon_num,1)));
fprintf('M02 pr min = %f\n',nanmin(reshape(pr_us(:,:,5),ylat_num*xlon_num,1)));
fprintf('M02 tmax max = %f\n',nanmax(reshape(tmax_us(:,:,5),ylat_num*xlon_num,1)));
fprintf('M02 tmax min = %f\n',nanmin(reshape(tmax_us(:,:,5),ylat_num*xlon_num,1)));
fprintf('M02 tmin max = %f\n',nanmax(reshape(tmin_us(:,:,5),ylat_num*xlon_num,1)));
fprintf('M02 tmin min = %f\n',nanmin(reshape(tmin_us(:,:,5),ylat_num*xlon_num,1)));

fprintf('\npr bias\n');
fprintf('------------\n');
for p = 1:4
    fprintf('Avg over CONUS for %s = %f\n',prdct{p},nanmean(reshape(pr_us(:,:,p)-pr_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Max over CONUS for %s = %f\n',prdct{p},nanmax(reshape(pr_us(:,:,p)-pr_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Min over CONUS for %s = %f\n',prdct{p},nanmin(reshape(pr_us(:,:,p)-pr_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Std over CONUS for %s = %f\n\n',prdct{p},nanstd(reshape(pr_us(:,:,p)-pr_us(:,:,5),ylat_num*xlon_num,1)));
end
fprintf('\ntmin bias\n');
fprintf('------------\n');
for p = 1:4
    fprintf('Avg over CONUS for %s = %f\n',prdct{p},nanmean(reshape(tmin_us(:,:,p)-tmin_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Max over CONUS for %s = %f\n',prdct{p},nanmax(reshape(tmin_us(:,:,p)-tmin_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Min over CONUS for %s = %f\n',prdct{p},nanmin(reshape(tmin_us(:,:,p)-tmin_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Std over CONUS for %s = %f\n\n',prdct{p},nanstd(reshape(tmin_us(:,:,p)-tmin_us(:,:,5),ylat_num*xlon_num,1)));
end
fprintf('\ntmax bias\n');
fprintf('------------\n');
for p = 1:4
    fprintf('Avg over CONUS for %s = %f\n',prdct{p},nanmean(reshape(tmax_us(:,:,p)-tmax_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Max over CONUS for %s = %f\n',prdct{p},nanmax(reshape(tmax_us(:,:,p)-tmax_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Min over CONUS for %s = %f\n',prdct{p},nanmin(reshape(tmax_us(:,:,p)-tmax_us(:,:,5),ylat_num*xlon_num,1)));
    fprintf('Std over CONUS for %s = %f\n\n',prdct{p},nanstd(reshape(tmax_us(:,:,p)-tmax_us(:,:,5),ylat_num*xlon_num,1)));
end
%%
% -------------------------------------------------------------------------
% Plot setting 
% -------------------------------------------------------------------------
%Import color map for surface map
%mycolors  = colortable('/glade/u/home/mizukami/mtl/myColortables/WhBlGrYeRe.gp',3,10);
mycolors  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/rainbow.gp',3,10);
mycolors  = (1/255).*mycolors;
mycolors1  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/BlRe.rgb',3,2);
mycolors1  = (1/255).*mycolors1;
mycolors2  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/temp_diff_18lev.rgb',3,6);
mycolors2  = (1/255).*mycolors2;
% line color order
% clr=[0.85 0.6 0.05;0.95 0.8 0.55;0 0.75 0;0.95 0.95 0;0.6 0.6 0.6];
clr=[0.90 0.5 0.05;0.95 0.8 0.55;0 0.75 0;0.15 0.15 0.90;0.5 0.5 0.5];
% -------------------------------------------------------------------------
% Make plot 
% -------------------------------------------------------------------------
%% plot 1 - tair difference
figure('Units','inches','position',[1 4 5.75 6.5],'Visible','on','Color',[1 1 1]);
%M02
for i = 1:3
  subplot(5,3,i)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');  
  if i ==1;
      caxis([0 2000]);
      pcolorm(lat2d, lon2d, pr_us(:,:,5)); 
      title('P [mm/yr]','FontSize',9);
  elseif i ==2;
      caxis([-15 20]);
      pcolorm(lat2d, lon2d, tmin_us(:,:,5));
      title('Tmin [\circC]','FontSize',9);
  elseif i ==3;
      caxis([0 30]);
      pcolorm(lat2d, lon2d, tmax_us(:,:,5));
      title('Tmax [\circC]','FontSize',9);
  end    
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  colormap(mycolors);freezeColors; 
  mlabel('off'); plabel('off');
  framem off  
end
%BCCA-M02
for i = 1:3
  subplot(5,3,i+3)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1;
      caxis([-200 200]);
      pcolorm(lat2d, lon2d, pr_us(:,:,1)-pr_us(:,:,5));    
  elseif i ==2;
      caxis([-0.6 0.6]);
      pcolorm(lat2d, lon2d, tmin_us(:,:,1)-tmin_us(:,:,5));
  elseif i ==3;
      caxis([-0.6 0.6]);
      pcolorm(lat2d, lon2d, tmax_us(:,:,1)-tmax_us(:,:,5));
  end  
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  if i == 1;colormap(flipud(mycolors2)); else colormap(mycolors2); end
  freezeColors; mlabel('off'); plabel('off');
  framem off
end
%BCSDd - M02
for i = 1:3
  subplot(5,3,i+6)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1;
      caxis([-200 200]);
      pcolorm(lat2d, lon2d, pr_us(:,:,2)-pr_us(:,:,5));  
  elseif i ==2;
      caxis([-0.6 0.6]);
      pcolorm(lat2d, lon2d, tmin_us(:,:,2)-tmin_us(:,:,5));
  elseif i ==3;
      caxis([-0.6 0.6]);
      pcolorm(lat2d, lon2d, tmax_us(:,:,2)-tmax_us(:,:,5));
  end  
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  if i == 1;colormap(flipud(mycolors2)); else colormap(mycolors2); end
  freezeColors; mlabel('off'); plabel('off');
  framem off
end
%BCSDm - M02
for i = 1:3
  subplot(5,3,i+9)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1;
      caxis([-200 200]);
      pcolorm(lat2d, lon2d, pr_us(:,:,3)-pr_us(:,:,5));
  elseif i ==2;
      caxis([-0.6 0.6]);
      pcolorm(lat2d, lon2d, tmin_us(:,:,3)-tmin_us(:,:,5));
  elseif i ==3;
      caxis([-0.6 0.6]);
      pcolorm(lat2d, lon2d, tmax_us(:,:,3)-tmax_us(:,:,5));
  end
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  if i == 1;colormap(flipud(mycolors2)); else colormap(mycolors2); end
  freezeColors; mlabel('off'); plabel('off');
  framem off 
end
%AR-M02
for i = 1:3
    subplot(5,3,i+12)
    ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
    set(ax2, 'Visible','off','layer','top');
    setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
    if i ==1;
        caxis([-200 200]);
        pcolorm(lat2d, lon2d, pr_us(:,:,4)-pr_us(:,:,5));
    elseif i ==2;
        caxis([-0.6 0.6]);
        pcolorm(lat2d, lon2d, tmin_us(:,:,4)-tmin_us(:,:,5));
    elseif i ==3;
        caxis([-0.6 0.6]);
        pcolorm(lat2d, lon2d, tmax_us(:,:,4)-tmax_us(:,:,5));
    end
    plotm([S.Y],[S.X], 'k');
    tightmap;
    if i == 1;colormap(flipud(mycolors2)); else colormap(mycolors2); end
    freezeColors; mlabel('off'); plabel('off');
    framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 =+0.020;
hmove2 =+0.000;
hmove3 =-0.020;
vmove1 =-0.044;
vmove2 =-0.044;
vmove3 =-0.044;
vmove4 =-0.044;
vmove5 =+0.022;
hscale=1.225;
vscale=1.225;

%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])
%the 4th row fromright to left
set(h(4), 'Position',[ pos(4,1)+hmove1  pos(4,2)+vmove2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove2  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove3  pos(6,2)+vmove2  pos(6,3)*hscale pos(6,4)*vscale])
%The 3rd row fromright to left
set(h(7), 'Position',[ pos(7,1)+hmove1  pos(7,2)+vmove3  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hmove2  pos(8,2)+vmove3  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove3  pos(9,2)+vmove3  pos(9,3)*hscale pos(9,4)*vscale])
%The 2nd row fromright to left
set(h(10),'Position',[pos(10,1)+hmove1 pos(10,2)+vmove4 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)+hmove2 pos(11,2)+vmove4 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)+hmove3 pos(12,2)+vmove4 pos(12,3)*hscale pos(12,4)*vscale])
%top row fromright to left
set(h(13),'Position',[pos(13,1)+hmove1 pos(13,2)+vmove5 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)+hmove2 pos(14,2)+vmove5 pos(14,3)*hscale pos(14,4)*vscale])
set(h(15),'Position',[pos(15,1)+hmove3 pos(15,2)+vmove5 pos(15,3)*hscale pos(15,4)*vscale])
%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.05,0.900,'M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.660,'BCCA-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.490,'BCSDd-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.310,'BCSDm-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.05,0.140,'AR-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);

%M02 pr colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 2000]); colormap(mycolors);
colorbar('horizontal','position',[0.125,0.785,0.185,0.0225],...
   'XLim',[0 2000],'XTick',0:500:2000,'XTickLabel',0:0.5:2,'FontSize',8);
text(0.350,0.770,'x10^4','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
cbfreeze;

%M02 tmin colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-15 25]); colormap(mycolors);
colorbar('horizontal','position',[0.435,0.785,0.185,0.0225],...
   'XLim',[-15 25],'XTick',-15:10:25,'XTickLabel',-15:10:25,'FontSize',8);
cbfreeze;
%M02 tmax colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 30]); colormap(mycolors);
colorbar('horizontal','position',[0.735,0.785,0.185,0.0225],...
    'XLim',[0 30],'XTick',0:7.5:30,'XTickLabel',0:7.5:30,'FontSize',8);
cbfreeze;

% pr difference from M02 colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-200 200]); colormap(flipud(mycolors2));
colorbar('horizontal','position',[0.125,0.030,0.185,0.0225],...
   'XLim',[-202 202],'XTick',-200:100:200,'XTickLabel',-0.2:0.1:0.2,'FontSize',8);
text(0.355,0.015,'x10^4','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
cbfreeze;
% tmin difference from M02 colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-0.6 0.6]); colormap(mycolors2);
colorbar('horizontal','position',[0.435,0.030,0.185,0.0225],...
   'XLim',[-0.6 0.6],'XTick',-0.6:0.3:0.6,'XTickLabel',-0.6:0.3:0.6,'FontSize',8);
cbfreeze;
% tmax difference from M02 colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-0.6 0.6]); colormap(mycolors2);
colorbar('horizontal','position',[0.735,0.030,0.185,0.0225],...
   'XLim',[-0.6 0.6],'XTick',-0.6:0.3:0.6,'XTickLabel',-0.6:0.3:0.6,'FontSize',8);
cbfreeze;

figfile=['./figure/paper/Paper_m2_fig1_annual_sd_data_diff_WYR' num2str(wyr(1)) '-' num2str(wyr(end)) '_' model '_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-r300',figfile);
% print('-depsc2',figfile);
toc