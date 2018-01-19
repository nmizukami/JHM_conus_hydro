% Comparison of annual water balance (ET and RO) among the simulations 
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
%  Fig 3-1 map   daily      annual      yes          3 elev band + snotel   Tair
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

nhuc=length(region);
%% Getting daily clm datm variable from various product files
% -------------------------------------------------------------------------
% Memory allocation
% -------------------------------------------------------------------------
% dimension [huc x region x prod]
pr_huc = zeros(nhuc,length(force),length(model));
et_huc = zeros(nhuc,length(force),length(model));
ro_huc = zeros(nhuc,length(force),length(model));

for r = 1:nhuc
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
          var0(var0==fillvalue)=NaN;
          pr_huc(r,p,m) = nanmean(nanmean(var0*60*60*24*365,1),2);
          clear var0 ncname
          
          % Annual RO
          ncname = [main_indir '/' force{p} '_' region{r} '/calib/Annual_climate_RUNOFF.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'RUNOFF','attname1','_FillValue');
          var0(var0==fillvalue)=NaN;
          var0(var0>3000/60/60/24/365)=NaN;
          ro_huc(r,p,m) = nanmean(nanmean(var0*60*60*24*365,1),2);
          clear var0 ncname
          
          % Annual ET
          ncname = [main_indir '/' force{p} '_' region{r} '/calib/Annual_climate_ET.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'ET','attname1','_FillValue');
          var0(var0==fillvalue)=NaN;
          var0(var0<0)=NaN;
          et_huc(r,p,m) = nanmean(nanmean(var0*60*60*24*365,1),2);
          clear var0 ncname
      end
  end
end
%% Further manipulation
% -------------------------------------------------------------------------
% Plot setting 
% -------------------------------------------------------------------------
% line color order
clr=[0.0 0.0 1.0;0.5 0.2 0.5;1.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 0.0];
% -------------------------------------------------------------------------
% Make plot 
% -------------------------------------------------------------------------
%% plot 1 - inter model difference WB difference
figure('Units','inches','position',[1 4 5.75 5.75],'Visible','on','Color',[1 1 1]);
huc_list=[5;10;14;17];
for i=1:4
    subplot(2,2,i)
    set(gca, 'Visible','off','layer','top');
    for f=1:length(force)
        plot(0:10:pr_huc(huc_list(i),f,1),pr_huc(huc_list(i),f,1):-10:0,'-','color',clr(f,:));hold on
        text(ro_huc(huc_list(i),f,1),et_huc(huc_list(i),f,1),'c','color',clr(f,:),'FontSize',11);
        text(ro_huc(huc_list(i),f,2),et_huc(huc_list(i),f,2),'v','color',clr(f,:),'FontSize',11);
        text(ro_huc(huc_list(i),f,3),et_huc(huc_list(i),f,3),'p','color',clr(f,:),'FontSize',11);
    end
    if i==1; title('OH','FontSize',10); set(gca,'XLim',[100 500]);set(gca,'YLim',[550 850]);
    elseif i==2;title('AR','FontSize',10);set(gca,'XLim',[0 400]);set(gca,'YLim',[300 600]);
    elseif i==3;title('UCO','FontSize',10);set(gca,'XLim',[0 400]);set(gca,'YLim',[150 450]);
    elseif i==4;title('PN','FontSize',10);set(gca,'XLim',[200 600]);set(gca,'YLim',[200 500]);
    end
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.005;
hmove2 = -0.015;
vmove1 = -0.0315;
vmove2 = -0.0300;

hscale=1.1;
vscale=1.1;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove1  pos(3,2)+vmove2  pos(3,3)*hscale pos(3,4)*vscale])
set(h(4), 'Position',[ pos(4,1)+hmove2  pos(4,2)+vmove2  pos(4,3)*hscale pos(4,4)*vscale])

%legends 1
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr(3,:));
dummy4 = plot(1,1,'-','LineWidth',2,'color',clr(4,:));
dummy5 = plot(1,1,'-','LineWidth',2,'color',clr(5,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off'); set(dummy5,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('BCCA','BCSDd','BCSDm','AR','M02',...
    'Orientation','vertical','Location',[0.850 0.825 0.0075 0.030]);
legend('boxoff');

%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.035,0.475,'Annual ET [mm]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',90);
text(0.55,0.025,'Annual Q [mm]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);


%save figure
figfile=['./figure/paper/Paper_m13_fig1_ET_Q_wyr' num2str(wyr(1)) '-' num2str(wyr(end)) '_huc_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto');
print('-dtiff','-r300',figfile);
% set(gcf, 'Renderer','OpenGL');
% print('-depsc2',figfile);


toc