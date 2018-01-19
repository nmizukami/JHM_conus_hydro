% Comparison of annual cycle of Tmax, Tmin, and DTR between two dataset (NLDAS vs. VIC)
% 
% This code maps climatological annual mean flux (precipitation, runoff, and
% ET).   Climatology is period beteween WY1981 through WY2008.
% 
% If climatology over other periods are need, go back to python scripts and recompute 
% 
%          plot_type     temp_res.  time span   climatology  spatial res.           variables
%  ----------------------------------------------------------------------------------------
%  Fig 1-1 map   daily      annual      yes          3 elev band + snotel   Tair
%  Fig 2-1 map   daily      annual      yes          3 elev band + snotel   Tair
%  Fig 3-1 map  daily      annual      yes          3 elev band + snotel   Tair
% 
% This script does not use CLM history file. Instead, use forcing dataset
% (netCDF prepared in preprocess)

tic

close all
clear
clc
%% Definition
% -------------------------------------------------------------------------
% Forcing product list (in orde of increasing resolution)
% -------------------------------------------------------------------------
force{1,1}='BCCA12K';
force{2,1}='BCSD12K';
force{3,1}='BCSDdisag12K';
force{4,1}='BCSAR12K';
force{5,1}='MAURER12K';
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
model{1,1} = 'CLM';
model{2,1} = 'VIC';
model{3,1} = 'PRMS';
% -------------------------------------------------------------------------
% resolution
% -------------------------------------------------------------------------
res = 12;
% -------------------------------------------------------------------------
% Starting and ending time period for analysis
% -------------------------------------------------------------------------
byr=1980; bmon=10; bday=1;
eyr=1999; emon=10; eday=1;
% eyr=2008; emon=9; eday=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Variables that cannot be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Read western domain grid 
% -------------------------------------------------------------------------
us_domain_grid = ['/d3/mizukami/domain_huc/domain_conus_' num2str(res) 'k.nc'];
lon2d = netcdf2mat(us_domain_grid,'x'); lon2d(lon2d == -999) = NaN;
lat2d = netcdf2mat(us_domain_grid,'y');
lon2d = 360+lon2d;
lon1d = netcdf2mat(us_domain_grid,'lon');
lat1d = netcdf2mat(us_domain_grid,'lat');
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
% dimension [lat x lon x region x prod]
pr_us = zeros(ylat_num,xlon_num,length(force),length(model));
et_us = zeros(ylat_num,xlon_num,length(force),length(model));
ro_us = zeros(ylat_num,xlon_num,length(force),length(model));
rr_us = zeros(ylat_num,xlon_num,length(force),length(model));

%Create simulation area mask
sim_mask = lon2d+huc2g;
sim_mask(~isnan(sim_mask) ) = 1;
sim_mask = repmat(sim_mask,[1,1,length(force),length(model)]);

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
  % Read western domain grid 
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
  % Annual climatology value
  % [lat x lon x prod]
  for p=1:length(force)
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
          % Annual PRCP
          ncname = [main_indir '/' force{p} '_' region{r} '/calib/Annual_climate_PRCP.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'PRCP','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          pr_us(i1:i2,j1:j2,p,m) = var0*60*60*24*365+pr_us(i1:i2,j1:j2,p,m);
          clear var0 ncname
          
          % Annual RO
          ncname = [main_indir '/' force{p} '_' region{r} '/calib/Annual_climate_RUNOFF.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'RUNOFF','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          ro_us(i1:i2,j1:j2,p,m) = var0*60*60*24*365+ro_us(i1:i2,j1:j2,p,m);
          clear var0 ncname
          
          % Annual ET
          ncname = [main_indir '/' force{p} '_' region{r} '/calib/Annual_climate_ET.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'ET','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          et_us(i1:i2,j1:j2,p,m) = var0*60*60*24*365+et_us(i1:i2,j1:j2,p,m);
          clear var0 ncname
          
          % Annual Runoff ratio
          ncname = [main_indir '/' force{p} '_' region{r} '/calib/Annual_climate_RR.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'RR','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          rr_us(i1:i2,j1:j2,p,m) = var0+rr_us(i1:i2,j1:j2,p,m);
          clear var0 ncname
      end
  end
end
%% Further manipulation
%Mask out the outside the simulation domain
pr_us = pr_us.*sim_mask;
et_us = et_us.*sim_mask;
ro_us = ro_us.*sim_mask;
rr_us = rr_us.*sim_mask;

%Compute stat
pr_M02_max = round(nanmax( nanmax( pr_us(:,:,5),[],1),[],2));
pr_M02_min = round(nanmin( nanmin( pr_us(:,:,5),[],1),[],2));
et_M02_max = round(nanmax( nanmax( et_us(:,:,5),[],1),[],2));
et_M02_min = round(nanmin( nanmin( et_us(:,:,5),[],1),[],2));
ro_M02_max = round(nanmax( nanmax( ro_us(:,:,5),[],1),[],2));
ro_M02_min = round(nanmin( nanmin( ro_us(:,:,5),[],1),[],2));
rr_M02_max = round(nanmax( nanmax( rr_us(:,:,5),[],1),[],2));
rr_M02_min = round(nanmin( nanmin( rr_us(:,:,5),[],1),[],2));

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
mycolors2 = flipud(mycolors2);

% line color order
clr=[0.90 0.5 0.05;0.95 0.8 0.55;0 0.75 0;0.15 0.15 0.90;0.5 0.5 0.5];
% -------------------------------------------------------------------------
% Make plot 
% -------------------------------------------------------------------------
%% plot 1 - WB difference
cbar_bcca = [-200 200];
cbar = [-100 100];
figure('Units','inches','position',[1 4 7.25 5.75],'Visible','on','Color',[1 1 1]);
%M02
for i = 1:7
  subplot(5,7,i)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1; 
      pcolorm(lat2d, lon2d, pr_us(:,:,5,1));
      caxis([0 2000]);
  elseif i >=2 && i <=4; 
      pcolorm(lat2d, lon2d, et_us(:,:,5,i-1));
      caxis([0 2000]);
      if i==2; title('CLM','FontSize',9); elseif i==3;title('VIC','FontSize',9); elseif i==4; title('PRMS','FontSize',9);end
  elseif i >=5 && i <=7 ; 
      pcolorm(lat2d, lon2d, ro_us(:,:,5,i-4));
      caxis([0 2000]);
      if i==5; title('CLM','FontSize',9); elseif i==6;title('VIC','FontSize',9); elseif i==7; title('PRMS','FontSize',9);end
  end
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  colormap(mycolors);freezeColors;
  framem off
end
%BCCA-M02
for i = 1:7
  subplot(5,7,i+7)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1;
      pcolorm(lat2d, lon2d, pr_us(:,:,1,1)-pr_us(:,:,5,1));
      caxis(cbar_bcca);
  elseif i >=2 && i <=4;
      pcolorm(lat2d, lon2d, et_us(:,:,1,i-1)-et_us(:,:,5,i-1));
      caxis(cbar_bcca);
  elseif i >=5 && i <=7;
      pcolorm(lat2d, lon2d, ro_us(:,:,1,i-4)-ro_us(:,:,5,i-4));
      caxis(cbar_bcca);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem off
end
%BCSDd-M02
for i = 1:7
  subplot(5,7,i+14)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  framem off
  if i ==1;
      pcolorm(lat2d, lon2d, pr_us(:,:,2,1)-pr_us(:,:,5,1));
      caxis(cbar);
  elseif i >=2 && i <=4;
      pcolorm(lat2d, lon2d, et_us(:,:,2,i-1)-et_us(:,:,5,i-1));
      caxis(cbar);
  elseif i >=5 && i <=7;
      pcolorm(lat2d, lon2d, ro_us(:,:,2,i-4)-ro_us(:,:,5,i-4));
      caxis(cbar);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem off
end
%BCSDm-M02
for i = 1:7
  subplot(5,7,i+21)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1;
      pcolorm(lat2d, lon2d, pr_us(:,:,3,1)-pr_us(:,:,5,1));
      caxis(cbar);
  elseif i >=2 && i <=4;
      pcolorm(lat2d, lon2d, et_us(:,:,3,i-1)-et_us(:,:,5,i-1));
      caxis(cbar);
  elseif i >=5 && i <=7;
      pcolorm(lat2d, lon2d, ro_us(:,:,3,i-4)-ro_us(:,:,5,i-4));
      caxis(cbar);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap; colormap(mycolors2); freezeColors;
  framem off
end
%AR-MA02
for i = 1:7
  subplot(5,7,i+28)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i ==1;
      pcolorm(lat2d, lon2d, pr_us(:,:,4,1)-pr_us(:,:,5,1));
      caxis(cbar);
  elseif i >=2 && i <=4;
      pcolorm(lat2d, lon2d, et_us(:,:,4,i-1)-et_us(:,:,5,i-1));
      caxis(cbar);
  elseif i >=5 && i <=7;
      pcolorm(lat2d, lon2d, ro_us(:,:,4,i-4)-ro_us(:,:,5,i-4));
      caxis(cbar);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 =  0.035;
hmove2 =  0.025;
hmove3 =  0.015;
hmove4 = -0.020;
hmove5 = -0.030;
hmove6 = -0.040;
hmove7 = -0.075;
vmove1 = -0.0430;
vmove2 = -0.0465;
vmove3 = -0.0500;
vmove4 = -0.0535;
vmove5 = -0.0200;
hscale=1.375;
vscale=1.375;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])
set(h(4), 'Position',[ pos(4,1)+hmove4  pos(4,2)+vmove1  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove5  pos(5,2)+vmove1  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove6  pos(6,2)+vmove1  pos(6,3)*hscale pos(6,4)*vscale])
set(h(7), 'Position',[ pos(7,1)+hmove7  pos(7,2)+vmove1  pos(7,3)*hscale pos(7,4)*vscale])

set(h(8), 'Position',[ pos(8,1)+hmove1  pos(8,2)+vmove2  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove2  pos(9,2)+vmove2  pos(9,3)*hscale pos(9,4)*vscale])
set(h(10),'Position',[pos(10,1)+hmove3 pos(10,2)+vmove2 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)+hmove4 pos(11,2)+vmove2 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)+hmove5 pos(12,2)+vmove2 pos(12,3)*hscale pos(12,4)*vscale])
set(h(13),'Position',[pos(13,1)+hmove6 pos(13,2)+vmove2 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)+hmove7 pos(14,2)+vmove2 pos(14,3)*hscale pos(14,4)*vscale])

set(h(15),'Position',[pos(15,1)+hmove1 pos(15,2)+vmove3 pos(15,3)*hscale pos(15,4)*vscale])
set(h(16),'Position',[pos(16,1)+hmove2 pos(16,2)+vmove3 pos(16,3)*hscale pos(16,4)*vscale])
set(h(17),'Position',[pos(17,1)+hmove3 pos(17,2)+vmove3 pos(17,3)*hscale pos(17,4)*vscale])
set(h(18),'Position',[pos(18,1)+hmove4 pos(18,2)+vmove3 pos(18,3)*hscale pos(18,4)*vscale])
set(h(19),'Position',[pos(19,1)+hmove5 pos(19,2)+vmove3 pos(19,3)*hscale pos(19,4)*vscale])
set(h(20),'Position',[pos(20,1)+hmove6 pos(20,2)+vmove3 pos(20,3)*hscale pos(20,4)*vscale])
set(h(21),'Position',[pos(21,1)+hmove7 pos(21,2)+vmove3 pos(21,3)*hscale pos(21,4)*vscale])

set(h(22),'Position',[pos(22,1)+hmove1 pos(22,2)+vmove4 pos(22,3)*hscale pos(22,4)*vscale])
set(h(23),'Position',[pos(23,1)+hmove2 pos(23,2)+vmove4 pos(23,3)*hscale pos(23,4)*vscale])
set(h(24),'Position',[pos(24,1)+hmove3 pos(24,2)+vmove4 pos(24,3)*hscale pos(24,4)*vscale])
set(h(25),'Position',[pos(25,1)+hmove4 pos(25,2)+vmove4 pos(25,3)*hscale pos(25,4)*vscale])
set(h(26),'Position',[pos(26,1)+hmove5 pos(26,2)+vmove4 pos(26,3)*hscale pos(26,4)*vscale])
set(h(27),'Position',[pos(27,1)+hmove6 pos(27,2)+vmove4 pos(27,3)*hscale pos(27,4)*vscale])
set(h(28),'Position',[pos(28,1)+hmove7 pos(28,2)+vmove4 pos(28,3)*hscale pos(28,4)*vscale])

set(h(29),'Position',[pos(29,1)+hmove1 pos(29,2)+vmove5 pos(29,3)*hscale pos(29,4)*vscale])
set(h(30),'Position',[pos(30,1)+hmove2 pos(30,2)+vmove5 pos(30,3)*hscale pos(30,4)*vscale])
set(h(31),'Position',[pos(31,1)+hmove3 pos(31,2)+vmove5 pos(31,3)*hscale pos(31,4)*vscale])
set(h(32),'Position',[pos(32,1)+hmove4 pos(32,2)+vmove5 pos(32,3)*hscale pos(32,4)*vscale])
set(h(33),'Position',[pos(33,1)+hmove5 pos(33,2)+vmove5 pos(33,3)*hscale pos(33,4)*vscale])
set(h(34),'Position',[pos(34,1)+hmove6 pos(34,2)+vmove5 pos(34,3)*hscale pos(34,4)*vscale])
set(h(35),'Position',[pos(35,1)+hmove7 pos(35,2)+vmove5 pos(35,3)*hscale pos(35,4)*vscale])

%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.025,0.885,'M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.025,0.665,'BCCA-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.025,0.500,'BCSDd-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.025,0.325,'BCSDm-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.025,0.175,'AR-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.125,0.975,'P [mm]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.375,0.975,'ET [mm]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.775,0.975,'RO [mm]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

% text(0.065,0.0575,num2str(cbar_bcca),'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

%M02 PR colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 2000]); colormap(mycolors);
colorbar('horizontal','position',[0.065,0.760,0.0750,0.0225],...
    'XLim',[0 2000],'XTick',0:1000:2000,'XTickLabel',0:1:2,'FontSize',8);
text(0.1750,0.750,'x10^4','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
cbfreeze;

%M02 ET colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 2000]); colormap(mycolors);
colorbar('horizontal','position',[0.300,0.760,0.1775,0.0225],...
    'XLim',[0 2000],'XTick',0:500:2000,'XTickLabel',0:0.5:2,'FontSize',8);
text(0.5120,0.750,'x10^4','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
cbfreeze;

%M02 RO colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 2000]); colormap(mycolors);
colorbar('horizontal','position',[0.700,0.760,0.1750,0.0225],...
    'XLim',[0 2000],'XTick',0:500:2000,'XTickLabel',0:0.5:2,'FontSize',8);
text(0.9120,0.750,'x10^4','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
cbfreeze;

% %For BCCA bias
% %PR difference colorbar, text on the plots
% axes('Units','Normalized','position',[0 0 1 1],'visible','off');
% caxis(cbar_bcca); colormap(mycolors2);
% colorbar('horizontal','position',[0.065,0.565,0.0775,0.0215],...
%    'XLim',[cbar_bcca(1)-1 cbar_bcca(2)+1],...
%    'XTick',cbar_bcca(1):cbar_bcca(2)-cbar_bcca(1):cbar_bcca(2),...
%    'XTickLabel',cbar_bcca(1):cbar_bcca(2)-cbar_bcca(1):cbar_bcca(2),'FontSize',8);
% cbfreeze;
% 
% %ET difference colorbar, text on the plots
% axes('Units','Normalized','position',[0 0 1 1],'visible','off');
% caxis(cbar_bcca); colormap(mycolors2);
% colorbar('horizontal','position',[0.300,0.565,0.1775,0.0215],...
%    'XLim',[cbar_bcca(1)-1 cbar_bcca(2)+1],...
%    'XTick',cbar_bcca(1):(cbar_bcca(2)-cbar_bcca(1))/2:cbar_bcca(2),...
%    'XTickLabel',cbar_bcca(1):(cbar_bcca(2)-cbar_bcca(1))/2:cbar_bcca(2),'FontSize',8);
% cbfreeze;
% 
% %RO difference colorbar, text on the plots
% axes('Units','Normalized','position',[0 0 1 1],'visible','off');
% caxis(cbar_bcca); colormap(mycolors2);
% colorbar('horizontal','position',[0.700,0.565,0.1775,0.0215],...
%    'XLim',[cbar_bcca(1)-1 cbar_bcca(2)+1],...
%    'XTick',cbar_bcca(1):(cbar_bcca(2)-cbar_bcca(1))/2:cbar_bcca(2),...
%    'XTickLabel',cbar_bcca(1):(cbar_bcca(2)-cbar_bcca(1))/2:cbar_bcca(2),'FontSize',8);
% cbfreeze;

%For other S.D bias
%PR difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis(cbar); colormap(mycolors2);
colorbar('horizontal','position',[0.075,0.0575,0.0775,0.0215],...
   'XLim',[cbar(1)-1 cbar(2)+1],...
   'XTick',cbar(1):cbar(2)-cbar(1):cbar(2),...
   'XTickLabel',cbar(1):cbar(2)-cbar(1):cbar(2),'FontSize',8);
cbfreeze;

%ET difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis(cbar); colormap(mycolors2);
colorbar('horizontal','position',[0.300,0.0575,0.1775,0.0215],...
   'XLim',[cbar(1)-1 cbar(2)+1],...
   'XTick',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),...
   'XTickLabel',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),'FontSize',8);
cbfreeze;

%RO difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis(cbar); colormap(mycolors2);
colorbar('horizontal','position',[0.700,0.0575,0.1775,0.0215],...
   'XLim',[cbar(1)-1 cbar(2)+1],...
   'XTick',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),...
   'XTickLabel',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),'FontSize',8);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m4_fig1_annual_wb_wyr' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
% set(gcf,'PaperPositionMode','auto')
print('-dtiff','-r300',figfile);
% set(gcf, 'Renderer','OpenGL');
% print('-depsc2',figfile);
toc