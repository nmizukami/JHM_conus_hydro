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
force{1,1}='MAURER12K';
force{2,1}='BCCA12K';
force{3,1}='BCSD12K';
force{4,1}='BCSDdisag12K';
force{5,1}='BCSAR12K';
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
%% plot 1 - inter model difference WB difference
cbar = [-200 200];
figure('Units','inches','position',[1 4 5.75 6.5],'Visible','on','Color',[1 1 1]);
c=0;
for f=1:length(force)
    for i = 1:3
        c=c+1;
        subplot(5,3,c)
        ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
        set(ax2, 'Visible','off','layer','top');
        setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
        if i ==1;
            pcolorm(lat2d, lon2d, ro_us(:,:,f,1)-ro_us(:,:,f,2));
            caxis(cbar);
        elseif i ==2 ;
            pcolorm(lat2d, lon2d, ro_us(:,:,f,1)-ro_us(:,:,f,3));
            caxis(cbar);
        elseif i ==3;
            pcolorm(lat2d, lon2d, ro_us(:,:,f,2)-ro_us(:,:,f,3));
            caxis(cbar);
        end
        if c==1; title('CLM-VIC','FontSize',9); 
        elseif c==2;title('CLM-PRMS','FontSize',9); 
        elseif c==3;title('VIC-PRMS','FontSize',9);
        end
        plotm([S.Y],[S.X], 'k');
        tightmap;
        colormap(mycolors2);freezeColors;
        framem off
    end
end
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.005;
hmove2 = -0.015;
hmove3 = -0.035;
vmove1 = -0.0375;
vmove2 = -0.0325;
vmove3 = -0.0275;
vmove4 = -0.0225;
vmove5 = -0.0175;
hscale=1.375;
vscale=1.375;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])

set(h(4), 'Position',[ pos(4,1)+hmove1  pos(4,2)+vmove2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove2  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove3  pos(6,2)+vmove2  pos(6,3)*hscale pos(6,4)*vscale])

set(h(7), 'Position',[ pos(7,1)+hmove1  pos(7,2)+vmove3  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hmove2  pos(8,2)+vmove3  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove3  pos(9,2)+vmove3  pos(9,3)*hscale pos(9,4)*vscale])

set(h(10),'Position',[pos(10,1)+hmove1 pos(10,2)+vmove4 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)+hmove2 pos(11,2)+vmove4 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)+hmove3 pos(12,2)+vmove4 pos(12,3)*hscale pos(12,4)*vscale])

set(h(13),'Position',[pos(13,1)+hmove1 pos(13,2)+vmove5 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)+hmove2 pos(14,2)+vmove5 pos(14,3)*hscale pos(14,4)*vscale])
set(h(15),'Position',[pos(15,1)+hmove3 pos(15,2)+vmove5 pos(15,3)*hscale pos(15,4)*vscale])

%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.055,0.885,'M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.055,0.725,'BCCA','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.055,0.525,'BCSDd','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.055,0.350,'BCSDm','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.055,0.175,'AR','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);

%RO difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis(cbar); colormap(mycolors2);
colorbar('horizontal','position',[0.350,0.0300,0.350,0.0200],...
   'XLim',[cbar(1)-1 cbar(2)+1],...
   'XTick',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),...
   'XTickLabel',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),'FontSize',9);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m12_fig1_intermodel_ro_wyr' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto');
print('-dtiff','-r300',figfile);
% set(gcf, 'Renderer','OpenGL');
% print('-depsc2',figfile);

%% plot 2 - inter model difference ET difference
cbar = [-200 200];
figure('Units','inches','position',[1 4 5.75 6.5],'Visible','on','Color',[1 1 1]);
c=0;
for f=1:length(force)
    for i = 1:3
        c=c+1;
        subplot(5,3,c)
        ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
        set(ax2, 'Visible','off','layer','top');
        setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
        if i ==1;
            pcolorm(lat2d, lon2d, et_us(:,:,f,1)-et_us(:,:,f,2));
            caxis(cbar);
        elseif i ==2 ;
            pcolorm(lat2d, lon2d, et_us(:,:,f,1)-et_us(:,:,f,3));
            caxis(cbar);
        elseif i ==3;
            pcolorm(lat2d, lon2d, et_us(:,:,f,2)-et_us(:,:,f,3));
            caxis(cbar);
        end
        if c==1; title('CLM-VIC','FontSize',9); 
        elseif c==2;title('CLM-PRMS','FontSize',9); 
        elseif c==3;title('VIC-PRMS','FontSize',9);
        end
        plotm([S.Y],[S.X], 'k');
        tightmap;
        colormap(mycolors2);freezeColors;
        framem off
    end
end
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.010;
hmove2 = -0.010;
hmove3 = -0.030;
vmove1 = -0.0400;
vmove2 = -0.0350;
vmove3 = -0.0300;
vmove4 = -0.0250;
vmove5 = -0.0200;
hscale=1.375;
vscale=1.375;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])

set(h(4), 'Position',[ pos(4,1)+hmove1  pos(4,2)+vmove2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove2  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove3  pos(6,2)+vmove2  pos(6,3)*hscale pos(6,4)*vscale])

set(h(7), 'Position',[ pos(7,1)+hmove1  pos(7,2)+vmove3  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hmove2  pos(8,2)+vmove3  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove3  pos(9,2)+vmove3  pos(9,3)*hscale pos(9,4)*vscale])

set(h(10),'Position',[pos(10,1)+hmove1 pos(10,2)+vmove4 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)+hmove2 pos(11,2)+vmove4 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)+hmove3 pos(12,2)+vmove4 pos(12,3)*hscale pos(12,4)*vscale])

set(h(13),'Position',[pos(13,1)+hmove1 pos(13,2)+vmove5 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)+hmove2 pos(14,2)+vmove5 pos(14,3)*hscale pos(14,4)*vscale])
set(h(15),'Position',[pos(15,1)+hmove3 pos(15,2)+vmove5 pos(15,3)*hscale pos(15,4)*vscale])

%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.035,0.885,'BCCA','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.035,0.725,'BCSDd','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.035,0.525,'BCSDm','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.035,0.350,'AR','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.035,0.175,'M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);

%RO difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis(cbar); colormap(mycolors2);
colorbar('horizontal','position',[0.375,0.0525,0.325,0.0215],...
   'XLim',[cbar(1)-1 cbar(2)+1],...
   'XTick',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),...
   'XTickLabel',cbar(1):(cbar(2)-cbar(1))/2:cbar(2),'FontSize',8);
cbfreeze;

%save figure
% figfile=['./figure/paper/Paper_m4_fig1_intermodel_et_wyr' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
% print('-dtiff','-r300',figfile);
% set(gcf, 'Renderer','OpenGL');
% print('-depsc2',figfile);


toc