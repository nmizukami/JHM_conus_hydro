% Comparison of annual wet-day fraction comparison among statistical
% downscaled products
%
%          plot_type     temp_res.  time span   climatology  spatial res.  variables
%  ----------------------------------------------------------------------------------------
%  Fig 1   scatter        daily      annual      yes          grid         Qmax20

tic

close all
clear
clc

prnt=1;

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
% Region
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
% byr=2001; bmon=10; bday=1;
eyr=2008; emon=9; eday=30;
% eyr=1999; emon=9; eday=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Variables that cannot be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Read western domain grid 
% -------------------------------------------------------------------------
us_domain_grid = ['/d3/mizukami/domain_huc/domain_conus_' num2str(res) 'k.nc'];
lon2d = netcdf2mat(us_domain_grid,'x');lon2d(lon2d == -999) = NaN;
lat2d = netcdf2mat(us_domain_grid,'y');
lon2d = 360+lon2d;
lon1d = netcdf2mat(us_domain_grid,'lon');
lat1d = netcdf2mat(us_domain_grid,'lat');
ele   = netcdf2mat(us_domain_grid,'ele');
huc2g = netcdf2mat(us_domain_grid,'huc2'); huc2g(huc2g == -999) = NaN;
nxlon = length(lon1d);
nylat = length(lat1d);
clear us_domain_grid

%Size of dimensions
nforce  = length(force);
nsd     = length(force)-1;
nmodel  = length(model);
nregion = length(region);
% -------------------------------------------------------------------------
% Read shapefile for western region HUC2
% -------------------------------------------------------------------------
huc2_shp = '/d3/mizukami/domain_huc/shapefile/HUC/HUC02_conus.shp';
S = shaperead(huc2_shp);
clear huc2_shp
% -------------------------------------------------------------------------
% water year array 
% -------------------------------------------------------------------------
wyr = byr+1:1:eyr;
%% Get data - CLM I/O and perform runoff process
% memory allocation

% grid [lat x lon x region x prod]
ro_mean_us      = zeros(nylat,nxlon,nforce,nmodel);
ro_var_us       = zeros(nylat,nxlon,nforce,nmodel);
ro7_10yr_min_us = zeros(nylat,nxlon,nforce,nmodel);
ro_20yr_max_us  = zeros(nylat,nxlon,nforce,nmodel);

sim_mask = lon2d+huc2g;
sim_mask(~isnan(sim_mask) ) = 1;
sim_mask = repmat(sim_mask,[1,1,nforce,nmodel]);

%mean areal value for each scale (HUC4, huc8, grid) [force x scale x region x model]
mean_ro_mean_us      =ones(4,nforce,nmodel)*NaN;
mean_ro_var_us       =ones(4,nforce,nmodel)*NaN;
mean_ro7_20yr_min_us =ones(4,nforce,nmodel)*NaN;
mean_ro_20yr_max_us  =ones(4,nforce,nmodel)*NaN;

%HUC2
ID_huc2_us           = [];
ro_mean_huc2_us      = zeros(nregion,nforce,nmodel);
ro_var_huc2_us       = zeros(nregion,nforce,nmodel);
ro_20yr_max_huc2_us  = zeros(nregion,nforce,nmodel);
ro7_10yr_min_huc2_us = zeros(nregion,nforce,nmodel);

%HUC4
ID_huc4_us           = [];
ro_mean_huc4_us      = cell(nregion,nforce,nmodel);
ro_var_huc4_us       = cell(nregion,nforce,nmodel);
ro_20yr_max_huc4_us  = cell(nregion,nforce,nmodel);
ro7_10yr_min_huc4_us = cell(nregion,nforce,nmodel);

%HUC8
ID_huc8_us           = [];
ro_mean_huc8_us      = cell(nregion,nforce,nmodel);
ro_var_huc8_us       = cell(nregion,nforce,nmodel);
ro_20yr_max_huc8_us  = cell(nregion,nforce,nmodel);
ro7_10yr_min_huc8_us = cell(nregion,nforce,nmodel);

for r = 1:nregion
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
    % Read HUC domain grid
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
    for f=1:nforce
        for m=1:nmodel
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
            
            ncname = [main_indir '/' force{f} '_' region{r} '/Extreme_RUNOFF.nc'];
            
            [var0,fillvalue]=netcdf2mat(ncname,'ro_mean','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro_mean_us(i1:i2,j1:j2,f,m) = var0+ro_mean_us(i1:i2,j1:j2,f,m);
            clear var0
            
            [var0,fillvalue]=netcdf2mat(ncname,'ro_var','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro_var_us(i1:i2,j1:j2,f,m) = var0+ro_var_us(i1:i2,j1:j2,f,m);
            clear var0
            
            [var0,fillvalue]=netcdf2mat(ncname,'ro7_10yr_min','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro7_10yr_min_us(i1:i2,j1:j2,f,m) = var0+ro7_10yr_min_us(i1:i2,j1:j2,f,m);
            clear var0
            
            [var0,fillvalue]=netcdf2mat(ncname,'ro_20yr_max','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro_20yr_max_us(i1:i2,j1:j2,f,m) = var0+ro_20yr_max_us(i1:i2,j1:j2,f,m);
            clear var0 ncname
            
            matname = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF_huc2.mat'];
            load(matname,'ro_mean_huc2','ro_var_huc2','ro7_10yr_min_huc2','ro_20yr_max_huc2');
            clear matname
            
            matname = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF_huc4.mat'];
            load(matname,'ro_mean_huc4','ro_var_huc4','ro7_10yr_min_huc4','ro_20yr_max_huc4');
            clear matname
            
            matname = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF_huc8.mat'];
            load(matname,'ro_mean_huc8','ro_var_huc8','ro7_10yr_min_huc8','ro_20yr_max_huc8');
            clear matname
            
            %Concatenate HUC runoff stat to all the regions
            ro_mean_huc2_us(r,f,m)      = ro_mean_huc2';
            ro_var_huc2_us(r,f,m)       = ro_var_huc2';
            ro7_10yr_min_huc2_us(r,f,m) = ro7_10yr_min_huc2';
            ro_20yr_max_huc2_us(r,f,m)  = ro_20yr_max_huc2';
            
            ro_mean_huc4_us{r,f,m}      = ro_mean_huc4;
            ro_var_huc4_us{r,f,m}       = ro_var_huc4;
            ro7_10yr_min_huc4_us{r,f,m} = ro7_10yr_min_huc4;
            ro_20yr_max_huc4_us{r,f,m}  = ro_20yr_max_huc4;
            
            ro_mean_huc8_us{r,f,m}      = ro_mean_huc8;
            ro_var_huc8_us{r,f,m}       = ro_var_huc8;
            ro7_10yr_min_huc8_us{r,f,m} = ro7_10yr_min_huc8;
            ro_20yr_max_huc8_us{r,f,m}  = ro_20yr_max_huc8;
        end
    end
end
%Mask out the outside the simulation domain
ro_mean_us      = ro_mean_us.*sim_mask;
ro_var_us       = ro_var_us.*sim_mask;
ro7_10yr_min_us = ro7_10yr_min_us.*sim_mask;
ro_20yr_max_us  = ro_20yr_max_us.*sim_mask;

%% Compute stat
fprintf('\nCONUS wide statistics\n');
fprintf('\nM02\n');
fprintf('------------\n');
for m=1:3
    fprintf('%s\n',model{m})
    fprintf('M02 Qmax20yr max = %7.2f\n',nanmax(reshape(ro_20yr_max_us(:,:,5,m),nylat*nxlon,1)));
    fprintf('M02 Qmax20yr min = %7.2f\n',nanmin(reshape(ro_20yr_max_us(:,:,5,m),nylat*nxlon,1)));
    fprintf('M02 7Q10 max = %5.2f\n',nanmax(reshape(ro7_10yr_min_us(:,:,5,m),nylat*nxlon,1)));
    fprintf('M02 7Q10 min = %5.2f\n',nanmin(reshape(ro7_10yr_min_us(:,:,5,m),nylat*nxlon,1)));
end
for m=1:3
    fprintf('------------\n');
    fprintf('%s\n',model{m})
    fprintf('------------\n\n');
    fprintf('Qmax20yr bias\n');
    fprintf('------------\n');
    for p = 1:4
        fprintf('Avg over CONUS for %s = %7.2f\n',force{p},nanmean(reshape(ro_20yr_max_us(:,:,p,m)-ro_20yr_max_us(:,:,5,m),nylat*nxlon,1)));
        fprintf('Max over CONUS for %s = %7.2f\n',force{p},nanmax(reshape(ro_20yr_max_us(:,:,p,m)-ro_20yr_max_us(:,:,5,m),nylat*nxlon,1)));
        fprintf('Min over CONUS for %s = %7.2f\n',force{p},nanmin(reshape(ro_20yr_max_us(:,:,p,m)-ro_20yr_max_us(:,:,5,m),nylat*nxlon,1)));
        fprintf('Std over CONUS for %s = %7.2f\n',force{p},nanstd(reshape(ro_20yr_max_us(:,:,p,m)-ro_20yr_max_us(:,:,5,m),nylat*nxlon,1)));
    end
    fprintf('\n 7Q10 bias\n');
    fprintf('------------\n');
    for p = 1:4
        fprintf('Avg over CONUS for %s = %5.2f\n',force{p},nanmean(reshape(ro7_10yr_min_us(:,:,p,m)-ro7_10yr_min_us(:,:,5,m),nylat*nxlon,1)));
        fprintf('Max over CONUS for %s = %5.2f\n',force{p},nanmax(reshape(ro7_10yr_min_us(:,:,p,m)-ro7_10yr_min_us(:,:,5,m),nylat*nxlon,1)));
        fprintf('Min over CONUS for %s = %5.2f\n',force{p},nanmin(reshape(ro7_10yr_min_us(:,:,p,m)-ro7_10yr_min_us(:,:,5,m),nylat*nxlon,1)));
        fprintf('Std over CONUS for %s = %5.2f\n\n',force{p},nanstd(reshape(ro7_10yr_min_us(:,:,p,m)-ro7_10yr_min_us(:,:,5,m),nylat*nxlon,1)));
    end
end
% -------------------------------------------------------------------------
% Make plot
% -------------------------------------------------------------------------
%% joint histogram for inter-model extreme flow comparison
% hist3 will bin the data
%Qmax20yr
xiro_20yr_max = 0:2:180;
yiro_20yr_max = xiro_20yr_max;
[Xro_20yr_max,Yro_20yr_max] = meshgrid(xiro_20yr_max,yiro_20yr_max);
%7Q10
xiro7_10yr_min = 0:0.01:1.6;
yiro7_10yr_min = xiro7_10yr_min;
[Xro7_10yr_min,Yro7_10yr_min] = meshgrid(xiro7_10yr_min,yiro7_10yr_min);

%mean
xiro_mean = 0:0.1:10;
yiro_mean = xiro_mean;
[Xro_mean,Yro_mean] = meshgrid(xiro_mean,yiro_mean);

%Qmax20yr
%memory allocation
ro_20yr_max_model  = ones(nylat*nxlon,2,nforce,nchoosek(nmodel,2))*NaN;
Nro_20yr_max_model = ones(length(xiro_20yr_max),length(yiro_20yr_max),nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            ro_20yr_max_model(:,:,f,c)  = [reshape( ro_20yr_max_us(:,:,f,m2), nylat*nxlon,1) reshape( ro_20yr_max_us(:,:,f,m1), nylat*nxlon,1)];
            Nro_20yr_max_model(:,:,f,c) = hist3(ro_20yr_max_model(:,:,f,c),{xiro_20yr_max yiro_20yr_max});
        end
    end
end

%7Q10
%memory allocation
ro7_10yr_min_model  = ones(nylat*nxlon,2,nforce,nchoosek(nmodel,2))*NaN;
Nro7_10yr_min_model = ones(length(xiro7_10yr_min),length(yiro7_10yr_min),nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            ro7_10yr_min_model(:,:,f,c)  = [reshape( ro7_10yr_min_us(:,:,f,m2), nylat*nxlon,1) reshape( ro7_10yr_min_us(:,:,f,m1), nylat*nxlon,1)];
            Nro7_10yr_min_model(:,:,f,c) = hist3(ro7_10yr_min_model(:,:,f,c),{xiro7_10yr_min yiro7_10yr_min});
        end
    end
end

%mean
%memory allocation
ro_mean_model  = ones(nylat*nxlon,2,nforce,nchoosek(nmodel,2))*NaN;
Nro_mean_model = ones(length(xiro_mean),length(yiro_mean),nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            ro_mean_model(:,:,f,c)  = [reshape( ro_mean_us(:,:,f,m2), nylat*nxlon,1) reshape( ro_mean_us(:,:,f,m1), nylat*nxlon,1)];
            Nro_mean_model(:,:,f,c) = hist3(ro_mean_model(:,:,f,c),{xiro_mean yiro_mean});
        end
    end
end
%% joint histogram for inter-forcing extreme comparison
%memory allocation
%Qmax20yr
ro_20yr_max_force  = ones(nylat*nxlon,2,nchoosek(nforce,2),nmodel)*NaN;
Nro_20yr_max_force = ones(length(xiro_20yr_max),length(yiro_20yr_max),nchoosek(nforce,2),nmodel)*NaN;
for m = 1:nmodel
    c=0;
    for f1 = 1:nforce-1
        for f2 = f1+1:nforce
            c=c+1;
            ro_20yr_max_force(:,:,c,m)  = [reshape( ro_20yr_max_us(:,:,f2,m), nylat*nxlon,1) reshape( ro_20yr_max_us(:,:,f1,m), nylat*nxlon,1)];
            Nro_20yr_max_force(:,:,c,m) = hist3(ro_20yr_max_force(:,:,c,m),{xiro_20yr_max yiro_20yr_max});
        end
    end
end
ro_20yr_max_intf = ro_20yr_max_force(:,:,[1;2;3;5;6;8],:);
ro_20yr_max_bias = ro_20yr_max_force(:,:,[4;7;9;10],:);
Nro_20yr_max_intf = Nro_20yr_max_force(:,:,[1;2;3;5;6;8],:);
Nro_20yr_max_bias = Nro_20yr_max_force(:,:,[4;7;9;10],:);

%memory allocation
%7Q10
ro7_10yr_min_force  = ones(nylat*nxlon,2,nchoosek(nforce,2),nmodel)*NaN;
Nro7_10yr_min_force = ones(length(xiro7_10yr_min),length(yiro7_10yr_min),nchoosek(nforce,2),nmodel)*NaN;
for m = 1:nmodel
    c=0;
    for f1 = 1:nforce-1
        for f2 = f1+1:nforce
            c=c+1;
            ro7_10yr_min_force(:,:,c,m)  = [reshape( ro7_10yr_min_us(:,:,f2,m), nylat*nxlon,1) reshape( ro7_10yr_min_us(:,:,f1,m), nylat*nxlon,1)];
            Nro7_10yr_min_force(:,:,c,m) = hist3(ro7_10yr_min_force(:,:,c,m),{xiro7_10yr_min yiro7_10yr_min});
        end
    end
end

ro7_10yr_min_intf = ro7_10yr_min_force(:,:,[1;2;3;5;6;8],:);
ro7_10yr_min_bias = ro7_10yr_min_force(:,:,[4;7;9;10],:);
Nro7_10yr_min_intf = Nro7_10yr_min_force(:,:,[1;2;3;5;6;8],:);
Nro7_10yr_min_bias = Nro7_10yr_min_force(:,:,[4;7;9;10],:);
%% -------------------------------------------------------------------------
% Statistics - Spatial bias
% -------------------------------------------------------------------------
ro7_10yr_min_force(ro7_10yr_min_force>2000)=NaN;
ro_20yr_max_force(ro_20yr_max_force>2000 | ro_20yr_max_force<-100)=NaN;
ro7_10yr_min_model(ro7_10yr_min_model>2000)=NaN;
ro_20yr_max_model(ro_20yr_max_model>2000 | ro_20yr_max_model<-100)=NaN;

bias_ro_20yr_max_force=squeeze( nanmean(ro_20yr_max_force(:,1,:,:)-ro_20yr_max_force(:,2,:,:)) );
bias_ro7_10yr_min_force=squeeze( nanmean(ro7_10yr_min_force(:,1,:,:)-ro7_10yr_min_force(:,2,:,:)) );

bias_ro_20yr_max_model=squeeze( nanmean(ro_20yr_max_model(:,1,:,:)-ro_20yr_max_model(:,2,:,:)) );
bias_ro7_10yr_min_model=squeeze( nanmean(ro7_10yr_min_model(:,1,:,:)-ro7_10yr_min_model(:,2,:,:)) );

%% -------------------------------------------------------------------------
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

mycolor3=colormap(jet);
mycolor3(1,:)=[1 1 1];

% line color order
clr4=[0 1 1;1 0 0;0 1 0;0 0 1];
clr5=[0 0 1;0.65 0.16 0.16;1 0 0;0 1 0;0 0 0];

force_name{1,1}='BCCA';
force_name{2,1}='BCSDd';
force_name{3,1}='BCSDm';
force_name{4,1}='AR';
force_name{5,1}='M02';
%% bias maps
figure('Units','inches','position',[1 4 7.25 5.50],'Visible','on','Color',[1 1 1]);
%M02
for i = 1:6
  subplot(5,6,i)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2, 'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
 if i >=1 && i <=3; 
      pcolorm(lat2d, lon2d, ro_20yr_max_us(:,:,5,i));
      caxis([0 120]);
      if i==1; title('CLM','FontSize',9); elseif i==2;title('VIC','FontSize',9); elseif i==3; title('PRMS','FontSize',9);end
  elseif i >=4 && i <=6 ; 
      pcolorm(lat2d, lon2d, ro7_10yr_min_us(:,:,5,i-3));
      caxis([0 0.6]);
      if i==4; title('CLM','FontSize',9); elseif i==5;title('VIC','FontSize',9); elseif i==6; title('PRMS','FontSize',9);end
  end
  plotm([S.Y],[S.X], 'k');
  tightmap; 
  colormap(mycolors);freezeColors;
  framem off
end
%BCCA-M02
for i = 1:6
  subplot(5,6,i+6)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i >=1 && i <=3;    
      pcolorm(lat2d, lon2d, ro_20yr_max_us(:,:,1,i)-ro_20yr_max_us(:,:,5,i));
      caxis([-20 20]);
  elseif i >=4 && i <=6;
      pcolorm(lat2d, lon2d, ro7_10yr_min_us(:,:,1,i-3)-ro7_10yr_min_us(:,:,5,i-3));
      caxis([-0.08 0.08]);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem off
end
%BCSDd-M02
for i = 1:6
  subplot(5,6,i+12)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  framem off
  if i >=1 && i <=3;
      pcolorm(lat2d, lon2d, ro_20yr_max_us(:,:,2,i)-ro_20yr_max_us(:,:,5,i));
      caxis([-20 20]);
  elseif i >=4 && i <=6;
      pcolorm(lat2d, lon2d, ro7_10yr_min_us(:,:,2,i-3)-ro7_10yr_min_us(:,:,5,i-3));
      caxis([-0.08 0.08]);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap;  colormap(mycolors2); freezeColors;
  framem off
end
%BCSDm-M02
for i = 1:6
  subplot(5,6,i+18)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i >=1 && i <=3;
      pcolorm(lat2d, lon2d, ro_20yr_max_us(:,:,3,i)-ro_20yr_max_us(:,:,5,i));
      caxis([-20 20]);
  elseif i >=4 && i <=6;
      pcolorm(lat2d, lon2d, ro7_10yr_min_us(:,:,3,i-3)-ro7_10yr_min_us(:,:,5,i-3));
      caxis([-0.08 0.08]);
  end
  plotm([S.Y],[S.X], 'k');
  tightmap; colormap(mycolors2); freezeColors;
  framem off
end
%AR-MA02
for i = 1:6
  subplot(5,6,i+24)
  ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
  set(ax2, 'Visible','off','layer','top');
  setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
  if i >=1 && i <=3;
      pcolorm(lat2d, lon2d, ro_20yr_max_us(:,:,4,i)-ro_20yr_max_us(:,:,5,i));
      caxis([-20 20]);
  elseif i >=4 && i <=6;
      pcolorm(lat2d, lon2d, ro7_10yr_min_us(:,:,4,i-3)-ro7_10yr_min_us(:,:,5,i-3));
      caxis([-0.08 0.08]);
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
vmove1 = -0.050;
vmove2 = -0.020;
hscale=1.375;
vscale=1.375;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])
set(h(4), 'Position',[ pos(4,1)+hmove4  pos(4,2)+vmove1  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove5  pos(5,2)+vmove1  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove6  pos(6,2)+vmove1  pos(6,3)*hscale pos(6,4)*vscale])

set(h(7), 'Position',[ pos(7,1)+hmove1  pos(7,2)+vmove1  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hmove2  pos(8,2)+vmove1  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove3  pos(9,2)+vmove1  pos(9,3)*hscale pos(9,4)*vscale])
set(h(10),'Position',[pos(10,1)+hmove4 pos(10,2)+vmove1 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)+hmove5 pos(11,2)+vmove1 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)+hmove6 pos(12,2)+vmove1 pos(12,3)*hscale pos(12,4)*vscale])

set(h(13),'Position',[pos(13,1)+hmove1 pos(13,2)+vmove1 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)+hmove2 pos(14,2)+vmove1 pos(14,3)*hscale pos(14,4)*vscale])
set(h(15),'Position',[pos(15,1)+hmove3 pos(15,2)+vmove1 pos(15,3)*hscale pos(15,4)*vscale])
set(h(16),'Position',[pos(16,1)+hmove4 pos(16,2)+vmove1 pos(16,3)*hscale pos(16,4)*vscale])
set(h(17),'Position',[pos(17,1)+hmove5 pos(17,2)+vmove1 pos(17,3)*hscale pos(17,4)*vscale])
set(h(18),'Position',[pos(18,1)+hmove6 pos(18,2)+vmove1 pos(18,3)*hscale pos(18,4)*vscale])

set(h(19),'Position',[pos(19,1)+hmove1 pos(19,2)+vmove1 pos(19,3)*hscale pos(19,4)*vscale])
set(h(20),'Position',[pos(20,1)+hmove2 pos(20,2)+vmove1 pos(20,3)*hscale pos(20,4)*vscale])
set(h(21),'Position',[pos(21,1)+hmove3 pos(21,2)+vmove1 pos(21,3)*hscale pos(21,4)*vscale])
set(h(22),'Position',[pos(22,1)+hmove4 pos(22,2)+vmove1 pos(22,3)*hscale pos(22,4)*vscale])
set(h(23),'Position',[pos(23,1)+hmove5 pos(23,2)+vmove1 pos(23,3)*hscale pos(23,4)*vscale])
set(h(24),'Position',[pos(24,1)+hmove6 pos(24,2)+vmove1 pos(24,3)*hscale pos(24,4)*vscale])

set(h(25),'Position',[pos(25,1)+hmove1 pos(25,2)+vmove2 pos(25,3)*hscale pos(25,4)*vscale])
set(h(26),'Position',[pos(26,1)+hmove2 pos(26,2)+vmove2 pos(26,3)*hscale pos(26,4)*vscale])
set(h(27),'Position',[pos(27,1)+hmove3 pos(27,2)+vmove2 pos(27,3)*hscale pos(27,4)*vscale])
set(h(28),'Position',[pos(28,1)+hmove4 pos(28,2)+vmove2 pos(28,3)*hscale pos(28,4)*vscale])
set(h(29),'Position',[pos(29,1)+hmove5 pos(29,2)+vmove2 pos(29,3)*hscale pos(29,4)*vscale])
set(h(30),'Position',[pos(30,1)+hmove6 pos(30,2)+vmove2 pos(30,3)*hscale pos(30,4)*vscale])

%------ Figure color bar and Annotatio-------------------------------------%
%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.050,0.885,'M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.050,0.675,'BCCA-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.050,0.500,'BCSDd-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.050,0.325,'BCSDm-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.050,0.150,'AR-M02','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
text(0.325,0.975,'RO_{20yr} [mm/day]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.750,0.975,'7RO10 [mm/day]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

%M02 Qmax20yr colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 120]); colormap(mycolors);
colorbar('horizontal','position',[0.175,0.760,0.20,0.0225],...
    'XLim',[0 120],'XTick',0:30:120,'FontSize',8);
cbfreeze;

%M02 7Q10 colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 0.6]); colormap(mycolors);
colorbar('horizontal','position',[0.625,0.760,0.20,0.0225],...
    'XLim',[-0.005 0.605],'XTick',0:0.15:0.6,'FontSize',8);
cbfreeze;

%Qmax20yr difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-20 20]); colormap(mycolors2);
colorbar('horizontal','position',[0.175,0.035,0.20,0.0225],...
   'XLim',[-20.5 20.5],'XTick',-20:10:20,'FontSize',8);
cbfreeze;

%7Q10 difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([-0.08 0.08]); colormap(mycolors2);
colorbar('horizontal','position',[0.625,0.035,0.20,0.0225],...
   'XLim',[-0.081 0.081],'XTick',-0.08:0.04:0.08,'FontSize',8);
cbfreeze;

%save figure
if prnt
    figfile=['./figure/paper/Paper_m7_fig1_extreme_wyr' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng',figfile);
    print('-dpdf','-r300',figfile);
end
%%
figure('Units','inches','position',[1 4 6.75 5.50],'Visible','on','Color',[1 1 1]);
%BCCA-M02
for i = 1:6
  subplot(4,6,i)
  % Qmax20yr
  if i >=1 && i <=3;
      surface(Xro_20yr_max,Yro_20yr_max,log10(Nro_20yr_max_force(:,:,4,i)),'EdgeColor','none');
      hold on
      plot3(xiro_20yr_max,yiro_20yr_max,5*ones(length(xiro_20yr_max),1),'k:','LineWidth',2);
      plot3(ones(length(xiro_20yr_max),1)*50,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*100,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*150,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,5*ones(length(xiro_20yr_max),1),'k:');
      axis equal
      set(gca,'XLim',[0 160]);
      set(gca,'YLim',[0 160]);
      set(gca,'XTick',0:50:150);
      set(gca,'YTick',0:50:150);
      set(gca,'YTickLabel',{'0';'';'';'150'});
      set(gca,'XTickLabel',{'0';'';'';'150'});      
      text(75,-20,force_name{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-20,75,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
      if i==1; title('CLM','FontSize',9); elseif i==2;title('VIC','FontSize',9); elseif i==3; title('PRMS','FontSize',9);end    
  % 7Q10     
  elseif i >=4 && i <=6;
      surface(Xro7_10yr_min,Yro7_10yr_min,log10(Nro7_10yr_min_force(:,:,4,i-3)),'EdgeColor','none');
      hold on
      plot3(xiro7_10yr_min,yiro7_10yr_min,5*ones(length(xiro7_10yr_min),1),'k:','LineWidth',2);
      plot3(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,5*ones(length(xiro7_10yr_min),1),'k:');
      axis equal
      set(gca,'XLim',[0 0.9]);
      set(gca,'YLim',[0 0.9]);
      set(gca,'XTick',0:0.2:0.8);
      set(gca,'YTick',0:0.2:0.8); 
      set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
      set(gca,'XTickLabel',{'0';'';'';'';'0.8'});      
      text(0.4,-0.1,force_name{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-0.1,0.4,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
      if i==4; title('CLM','FontSize',9); elseif i==5;title('VIC','FontSize',9); elseif i==6; title('PRMS','FontSize',9);end
  end
  caxis([0 3]);  
  set(gca,'FontSize',8);
end
%BCSDd-M02
for i = 1:6
  subplot(4,6,i+6)
  if i >=1 && i <=3;
      surface(Xro_20yr_max,Yro_20yr_max,log10(Nro_20yr_max_force(:,:,7,i)),'EdgeColor','none');
      hold on
      plot3(xiro_20yr_max,yiro_20yr_max,5*ones(length(xiro_20yr_max),1),'k:','LineWidth',2);
      plot3(ones(length(xiro_20yr_max),1)*50,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*100,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*150,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,5*ones(length(xiro_20yr_max),1),'k:');
      axis equal
      set(gca,'XLim',[0 160]);
      set(gca,'YLim',[0 160]);
      set(gca,'XTick',0:50:150);
      set(gca,'YTick',0:50:150);
      set(gca,'YTickLabel',{'0';'';'';'150'});
      set(gca,'XTickLabel',{'0';'';'';'150'});      
      text(75,-20,force_name{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-20,75,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);      
  elseif i >=4 && i <=6;
      surface(Xro7_10yr_min,Yro7_10yr_min,log10(Nro7_10yr_min_force(:,:,7,i-3)),'EdgeColor','none');
      hold on
      plot3(xiro7_10yr_min,yiro7_10yr_min,5*ones(length(xiro7_10yr_min),1),'k:','LineWidth',2);
      plot3(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,5*ones(length(xiro7_10yr_min),1),'k:');
      axis equal
      set(gca,'XLim',[0 0.9]);
      set(gca,'YLim',[0 0.9]);
      set(gca,'XTick',0:0.2:0.8);
      set(gca,'YTick',0:0.2:0.8); 
      set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
      set(gca,'XTickLabel',{'0';'';'';'';'0.8'});      
      text(0.4,-0.1,force_name{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-0.1,0.4,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);      
  end
  caxis([0 3]);  
  set(gca,'FontSize',8);
end
%BCSDm-M02
for i = 1:6
  subplot(4,6,i+12)
  if i >=1 && i <=3;
      surface(Xro_20yr_max,Yro_20yr_max,log10(Nro_20yr_max_force(:,:,9,i)),'EdgeColor','none');
      hold on
      plot3(xiro_20yr_max,yiro_20yr_max,5*ones(length(xiro_20yr_max),1),'k:','LineWidth',2);
      plot3(ones(length(xiro_20yr_max),1)*50,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*100,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*150,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,5*ones(length(xiro_20yr_max),1),'k:'); 
      axis equal
      set(gca,'XLim',[0 160]);
      set(gca,'YLim',[0 160]);
      set(gca,'XTick',0:50:150);
      set(gca,'YTick',0:50:150);
      set(gca,'YTickLabel',{'0';'';'';'150'});
      set(gca,'XTickLabel',{'0';'';'';'150'});      
      text(75,-20,force_name{3},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-20,75,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
  elseif i >=4 && i <=6;
      surface(Xro7_10yr_min,Yro7_10yr_min,log10(Nro7_10yr_min_force(:,:,9,i-3)),'EdgeColor','none');
      hold on
      plot3(xiro7_10yr_min,yiro7_10yr_min,5*ones(length(xiro7_10yr_min),1),'k:','LineWidth',2);
      plot3(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,5*ones(length(xiro7_10yr_min),1),'k:');
      axis equal
      set(gca,'XLim',[0 0.9]);
      set(gca,'YLim',[0 0.9]);
      set(gca,'XTick',0:0.2:0.8);
      set(gca,'YTick',0:0.2:0.8);  
      set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
      set(gca,'XTickLabel',{'0';'';'';'';'0.8'});      
      text(0.4,-0.1,force_name{3},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-0.1,0.4,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);      
  end
  caxis([0 3]);  
  set(gca,'FontSize',8);
end
%AR-MA02
for i = 1:6
  subplot(4,6,i+18)
  if i >=1 && i <=3;
      surface(Xro_20yr_max,Yro_20yr_max,log10(Nro_20yr_max_force(:,:,10,i)),'EdgeColor','none');
      hold on
      plot3(xiro_20yr_max,yiro_20yr_max,5*ones(length(xiro_20yr_max),1),'k:','LineWidth',2);
      plot3(ones(length(xiro_20yr_max),1)*50,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*100,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(ones(length(xiro_20yr_max),1)*150,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,5*ones(length(xiro_20yr_max),1),'k:');
      plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,5*ones(length(xiro_20yr_max),1),'k:');
      axis equal
      set(gca,'XLim',[0 160]);
      set(gca,'YLim',[0 160]);
      set(gca,'XTick',0:50:150);
      set(gca,'YTick',0:50:150);
      set(gca,'YTickLabel',{'0';'';'';'150'});
      set(gca,'XTickLabel',{'0';'';'';'150'});
      text(75,-20,force_name{4},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-20,75,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);      
  elseif i >=4 && i <=6;
      surface(Xro7_10yr_min,Yro7_10yr_min,log10(Nro7_10yr_min_force(:,:,10,i-3)),'EdgeColor','none');
      hold on
      plot3(xiro7_10yr_min,yiro7_10yr_min,5*ones(length(xiro7_10yr_min),1),'k:','LineWidth',2);
      plot3(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,5*ones(length(xiro7_10yr_min),1),'k:');
      plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,5*ones(length(xiro7_10yr_min),1),'k:');
      axis equal
      set(gca,'XLim',[0 0.9]);
      set(gca,'YLim',[0 0.9]);
      set(gca,'XTick',0:0.2:0.8);
      set(gca,'YTick',0:0.2:0.8);
      set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
      set(gca,'XTickLabel',{'0';'';'';'';'0.8'});
      text(0.4,-0.1,force_name{4},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-0.1,0.4,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
  end
  caxis([0 3]);
  set(gca,'FontSize',8);
end
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 =  0.0375;
hmove2 =  0.0250;
hmove3 =  0.0125;
hmove4 = -0.0350;
hmove5 = -0.0475;
hmove6 = -0.0600;
vmove1 = +0.010;
vmove2 = -0.000;
vmove3 = -0.015;
vmove4 = -0.030;
hscale=1.165;
vscale=1.165;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])
set(h(4), 'Position',[ pos(4,1)+hmove4  pos(4,2)+vmove1  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove5  pos(5,2)+vmove1  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove6  pos(6,2)+vmove1  pos(6,3)*hscale pos(6,4)*vscale])

set(h(7), 'Position',[ pos(7,1)+hmove1  pos(7,2)+vmove2  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hmove2  pos(8,2)+vmove2  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove3  pos(9,2)+vmove2  pos(9,3)*hscale pos(9,4)*vscale])
set(h(10),'Position',[pos(10,1)+hmove4 pos(10,2)+vmove2 pos(10,3)*hscale pos(10,4)*vscale])
set(h(11),'Position',[pos(11,1)+hmove5 pos(11,2)+vmove2 pos(11,3)*hscale pos(11,4)*vscale])
set(h(12),'Position',[pos(12,1)+hmove6 pos(12,2)+vmove2 pos(12,3)*hscale pos(12,4)*vscale])

set(h(13),'Position',[pos(13,1)+hmove1 pos(13,2)+vmove3 pos(13,3)*hscale pos(13,4)*vscale])
set(h(14),'Position',[pos(14,1)+hmove2 pos(14,2)+vmove3 pos(14,3)*hscale pos(14,4)*vscale])
set(h(15),'Position',[pos(15,1)+hmove3 pos(15,2)+vmove3 pos(15,3)*hscale pos(15,4)*vscale])
set(h(16),'Position',[pos(16,1)+hmove4 pos(16,2)+vmove3 pos(16,3)*hscale pos(16,4)*vscale])
set(h(17),'Position',[pos(17,1)+hmove5 pos(17,2)+vmove3 pos(17,3)*hscale pos(17,4)*vscale])
set(h(18),'Position',[pos(18,1)+hmove6 pos(18,2)+vmove3 pos(18,3)*hscale pos(18,4)*vscale])

set(h(19),'Position',[pos(19,1)+hmove1 pos(19,2)+vmove4 pos(19,3)*hscale pos(19,4)*vscale])
set(h(20),'Position',[pos(20,1)+hmove2 pos(20,2)+vmove4 pos(20,3)*hscale pos(20,4)*vscale])
set(h(21),'Position',[pos(21,1)+hmove3 pos(21,2)+vmove4 pos(21,3)*hscale pos(21,4)*vscale])
set(h(22),'Position',[pos(22,1)+hmove4 pos(22,2)+vmove4 pos(22,3)*hscale pos(22,4)*vscale])
set(h(23),'Position',[pos(23,1)+hmove5 pos(23,2)+vmove4 pos(23,3)*hscale pos(23,4)*vscale])
set(h(24),'Position',[pos(24,1)+hmove6 pos(24,2)+vmove4 pos(24,3)*hscale pos(24,4)*vscale])

% Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.750,0.955,'7RO10 [mm/day]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.325,0.955,'RO_{20yr} [mm/day]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
%colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 3]);
colorbar('horizontal','position',[0.775,0.050,0.100,0.0195],...
    'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
cbfreeze;
%save figure
if prnt
    figfile=['./figure/paper/Paper_m7_fig2_extreme_comparison_force_M02_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
    set(gcf,'PaperPositionMode','auto')
%     print('-dtiff','-r300',figfile);
    print('-dpdf','-r200',figfile);
end
%% Inter-model
figure('Units','inches','position',[1 4 6.00 4.25],'Visible','on','color',[1 1 1]);
% Qmax20yr
for i = 1:3
    subplot(2,3,i)
%     surface(Xro_20yr_max,Yro_20yr_max,log10(Nro_20yr_max_model(:,:,5,i)),'EdgeColor','none');
    imagesc(xiro_20yr_max+1,yiro_20yr_max+1,log10(Nro_20yr_max_model(:,:,5,i)) );
    caxis([0 3]);
    colormap(mycolor3)
    hold on
%     plot(xiro_20yr_max,yiro_20yr_max,'k:','LineWidth',2);
%     plot(ones(length(xiro_20yr_max),1)*50,linspace(0,200,length(xiro_20yr_max)),'k:');
%     plot(ones(length(xiro_20yr_max),1)*100,linspace(0,200,length(xiro_20yr_max)),'k:');
%     plot(ones(length(xiro_20yr_max),1)*150,linspace(0,200,length(xiro_20yr_max)),'k:');
%     plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,'k:');
%     plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,'k:');
%     plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,'k:');
    plot(xiro_20yr_max,yiro_20yr_max,'k:','LineWidth',2);
    plot(ones(length(xiro_20yr_max),1)*50,linspace(0,200,length(xiro_20yr_max)),'k:');
    plot(ones(length(xiro_20yr_max),1)*100,linspace(0,200,length(xiro_20yr_max)),'k:');
    plot(ones(length(xiro_20yr_max),1)*150,linspace(0,200,length(xiro_20yr_max)),'k:');
    plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,'k:');
    plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,'k:');
    plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,'k:');    
    axis equal
    axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0 180]);
    set(gca,'XTick',0:50:150);
    set(gca,'YTick',0:50:150);
    if i == 1
        text(75,-20,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-20,75,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(75,-20,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-20,75,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(75,-20,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-20,75,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end  
    set(gca,'YTickLabel',{'0';'';'';'150'});
    set(gca,'XTickLabel',{'0';'';'';'150'});
end
% 7Q10
for i = 1:3
    subplot(2,3,i+3)
%     surface(Xro7_10yr_min,Yro7_10yr_min,log10(Nro7_10yr_min_model(:,:,5,i)),'EdgeColor','none');
    imagesc(xiro7_10yr_min+0.005,yiro7_10yr_min+0.005,log10(Nro7_10yr_min_model(:,:,5,i)) );
    caxis([0 3]);
    colormap(mycolor3)
    hold on
%     plot3(xiro7_10yr_min,yiro7_10yr_min,5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,5*ones(length(xiro7_10yr_min),1),'k:');
%     plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,5*ones(length(xiro7_10yr_min),1),'k:');
    plot(xiro7_10yr_min,yiro7_10yr_min,'k:','LineWidth',2);
    plot(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,'k:');    
    axis equal
    axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 0.9]);
    set(gca,'YLim',[0 0.9]);
    set(gca,'XTick',0:0.2:0.8);
    set(gca,'YTick',0:0.2:0.8);
    if i == 1
        text(0.4,-0.1,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-0.1,0.4,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(0.4,-0.1,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-0.1,0.4,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(0.4,-0.1,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-0.1,0.4,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end
    set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
    set(gca,'XTickLabel',{'0';'';'';'';'0.8'});

end 
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
%from left to right
hmove1=+0.010;
hmove2=-0.000;
hmove3=-0.010;
vmove1=+0.001;
vmove2=-0.005;
hscale=1.075;
vscale=1.075;
%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])

set(h(4), 'Position',[ pos(4,1)+hmove1  pos(4,2)+vmove2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove2  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove3  pos(6,2)+vmove2  pos(6,3)*hscale pos(6,4)*vscale])
    
% Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.525,0.485,'7Q10 [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.525,0.950,'Q_{20yr} [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
%colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 3]);
colorbar('horizontal','position',[0.775,0.050,0.100,0.0195],...
    'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
cbfreeze;

%save figure
if prnt
    figfile=['./figure/paper/Paper_m7_fig3_extreme_comparison_model_M02_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
%     set(gcf,'PaperPositionMode','auto')
%     print('-dtiff','-r300',figfile);
    print('-dpdf','-r300',figfile);
end

%% Inter-model (including mean)
figure('Units','inches','position',[1 4 6.00 5.25],'Visible','on','color',[1 1 1]);
% Qmax20yr
for i = 1:3
    subplot(3,3,i)
    imagesc(xiro_20yr_max+1,yiro_20yr_max+1,log10(Nro_20yr_max_model(:,:,5,i)) );
    caxis([0 3]);
    colormap(mycolor3)
    hold on
    plot(xiro_20yr_max,yiro_20yr_max,'k:','LineWidth',2);
    plot(ones(length(xiro_20yr_max),1)*50,linspace(0,200,length(xiro_20yr_max)),'k:');
    plot(ones(length(xiro_20yr_max),1)*100,linspace(0,200,length(xiro_20yr_max)),'k:');
    plot(ones(length(xiro_20yr_max),1)*150,linspace(0,200,length(xiro_20yr_max)),'k:');
    plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,'k:');
    plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,'k:');
    plot(linspace(0,200,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,'k:');    
    axis equal
    axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0 180]);
    set(gca,'XTick',0:50:150);
    set(gca,'YTick',0:50:150);
    if i == 1
        text(75,-20,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-20,75,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(75,-20,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-20,75,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(75,-20,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-20,75,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end  
    set(gca,'YTickLabel',{'0';'';'';'150'});
    set(gca,'XTickLabel',{'0';'';'';'150'});
end
% 7Q10
for i = 1:3
    subplot(3,3,i+3)
    imagesc(xiro7_10yr_min+0.005,yiro7_10yr_min+0.005,log10(Nro7_10yr_min_model(:,:,5,i)) );
    caxis([0 3]);
    colormap(mycolor3)
    hold on
    plot(xiro7_10yr_min,yiro7_10yr_min,'k:','LineWidth',2);
    plot(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,'k:');
    plot(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,'k:');    
    axis equal
    axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 0.9]);
    set(gca,'YLim',[0 0.9]);
    set(gca,'XTick',0:0.2:0.8);
    set(gca,'YTick',0:0.2:0.8);
    if i == 1
        text(0.4,-0.1,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-0.1,0.4,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(0.4,-0.1,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-0.1,0.4,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(0.4,-0.1,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-0.1,0.4,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end
    set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
    set(gca,'XTickLabel',{'0';'';'';'';'0.8'});
end
for i=1:3
    % mean
    subplot(3,3,i+6)
    imagesc(xiro_mean,yiro_mean,log10(Nro_mean_model(:,:,5,i)) );
    caxis([0 3]);
    colormap(mycolor3)
    hold on
    plot(xiro_mean,yiro_mean,'k:','LineWidth',2);
    plot(ones(length(xiro_mean),1)*2,linspace(0,10.0,length(xiro_mean)),'k:');
    plot(ones(length(xiro_mean),1)*4,linspace(0,10.0,length(xiro_mean)),'k:');
    plot(ones(length(xiro_mean),1)*6,linspace(0,10.0,length(xiro_mean)),'k:');
    plot(ones(length(xiro_mean),1)*8,linspace(0,10.0,length(xiro_mean)),'k:');
    plot(linspace(0,10.0,length(xiro_mean)),ones(length(xiro_mean),1)*2,'k:');
    plot(linspace(0,10.0,length(xiro_mean)),ones(length(xiro_mean),1)*4,'k:');
    plot(linspace(0,10.0,length(xiro_mean)),ones(length(xiro_mean),1)*6,'k:');    
    plot(linspace(0,10.0,length(xiro_mean)),ones(length(xiro_mean),1)*8,'k:');    
    axis equal
    axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 9]);
    set(gca,'YLim',[0 9]);
    set(gca,'XTick',0:2:8);
    set(gca,'YTick',0:2:8);
    if i == 1
        text(4,-1,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-1,4,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(4,-1,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-1,4,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(4,-1,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-1,4,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end
    set(gca,'YTickLabel',{'0';'';'';'';'8'});
    set(gca,'XTickLabel',{'0';'';'';'';'8'});
end 
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
%from left to right
hmove1=+0.010;
hmove2=-0.000;
hmove3=-0.010;
vmove1=-0.005;
vmove2=-0.000;
vmove3=+0.005;
hscale=1.055;
vscale=1.055;
%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hmove1  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove2  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove1  pos(3,3)*hscale pos(3,4)*vscale])

set(h(4), 'Position',[ pos(4,1)+hmove1  pos(4,2)+vmove2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hmove2  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hmove3  pos(6,2)+vmove2  pos(6,3)*hscale pos(6,4)*vscale])
    
set(h(7), 'Position',[ pos(7,1)+hmove1  pos(7,2)+vmove3  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hmove2  pos(8,2)+vmove3  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hmove3  pos(9,2)+vmove3  pos(9,3)*hscale pos(9,4)*vscale])
% Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.525,0.350,'Mean RO [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.525,0.655,'7RO10 [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.525,0.965,'RO_{20yr} [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
%colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 3]);
colorbar('horizontal','position',[0.775,0.035,0.100,0.0195],...
    'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
cbfreeze;

%save figure
if prnt
    figfile=['./figure/paper/Paper_m7_fig3_extreme_mean_comparison_model_M02_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
    set(gcf,'PaperPositionMode','auto')
%     print('-dtiff','-r300',figfile);
    print('-dpdf','-r300',figfile);
end

%% Scatter plot of inter forcing comparison
for m=1:nmodel
    figure('Units','inches','position',[1 4 6.75 6.75],'Visible','on','Color',[1 1 1]);
    %  - upper right half
    c=0;
    for row=1:nsd-1
        c1=0;
        p1=nsd*(row-1)+row+1;
        for pan=p1:row*nsd
            c1=c1+1;
            c=c+1;
            subplot(4,4,pan)
            surface(Xro_20yr_max,Yro_20yr_max,log10(Nro_20yr_max_intf(:,:,c,m)),'EdgeColor','none');
            caxis([0 3]);
            hold on
            plot3(xiro_20yr_max,yiro_20yr_max,5*ones(length(xiro_20yr_max),1),'k:');
            plot3(ones(length(xiro_20yr_max),1)*50,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
            plot3(ones(length(xiro_20yr_max),1)*100,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
            plot3(ones(length(xiro_20yr_max),1)*150,linspace(0,160,length(xiro_20yr_max)),5*ones(length(xiro_20yr_max),1),'k:');
            plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*50,5*ones(length(xiro_20yr_max),1),'k:');
            plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*100,5*ones(length(xiro_20yr_max),1),'k:');
            plot3(linspace(0,160,length(xiro_20yr_max)),ones(length(xiro_20yr_max),1)*150,5*ones(length(xiro_20yr_max),1),'k:');
            axis equal
            set(gca,'FontSize',8);
            set(gca,'XLim',[0 160]);
            set(gca,'YLim',[0 160]);
            set(gca,'XTick',0:50:150);
            set(gca,'YTick',0:50:150);
            text(75,-20,force_name{row},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
            text(-20,75,force_name{c1+row},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
            if pan~=p1;
                set(gca,'YTickLabel',{''});
                set(gca,'XTickLabel',{''});
            else
                set(gca,'YTickLabel',{'0';'';'';'150'});
                set(gca,'XTickLabel',{'0';'';'';'150'});
            end
        end
    end
    %  - lower left half
    rowoffset=[5;10;15;20];
    c=0;
    for col=1:nsd-1
        p1 = nsd*col+col;
        c1=0;
        for pan=p1:nsd:(nsd-1)*nsd+col
            c=c+1;
            c1=c1+1;
            subplot(4,4,pan)
            surface(Xro7_10yr_min,Yro7_10yr_min,log10(Nro7_10yr_min_intf(:,:,c,m)),'EdgeColor','none');
            caxis([0 3]);
            hold on
            plot3(xiro7_10yr_min,yiro7_10yr_min,5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(ones(length(xiro7_10yr_min),1)*0.2,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(ones(length(xiro7_10yr_min),1)*0.4,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(ones(length(xiro7_10yr_min),1)*0.6,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(ones(length(xiro7_10yr_min),1)*0.8,linspace(0,1.0,length(xiro7_10yr_min)),5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.2,5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.4,5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.6,5*ones(length(xiro7_10yr_min),1),'k:');
            plot3(linspace(0,1.0,length(xiro7_10yr_min)),ones(length(xiro7_10yr_min),1)*0.8,5*ones(length(xiro7_10yr_min),1),'k:');
            axis equal
            set(gca,'FontSize',8);
            set(gca,'XLim',[0 0.9]);
            set(gca,'YLim',[0 0.9]);
            set(gca,'XTick',0:0.2:0.8);
            set(gca,'YTick',0:0.2:0.8);
            text(0.4,-0.1,force_name{col},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
            text(-0.1,0.4,force_name{c1+col},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
            if row~=1
                set(gca,'YTickLabel',{''});
            else
                set(gca,'YTickLabel',{'0';'';'';'';'0.8'});
            end
            if pan<20
                set(gca,'XTickLabel',{''});
            else
                set(gca,'XTickLabel',{'0';'';'';'';'0.8'});
            end
        end
    end
   
    %adjust plot size & location
    h   = get(gcf,'Children');
    for i=1:length(h)
        pos(i,:) =get(h(i),'position');
    end
    %from left to right
    hmove1=+0.015;
    hmove2=+0.000;
    hmove3=-0.015;
    hmove4=-0.030;
    hmove5=-0.045;
    %from bottom to top
    vmove1=-0.015;
    vmove2=-0.010;
    vmove3=-0.005;
    vmove4=+0.000;
    vmove5=+0.005;
    
    hscale=1.1275;
    vscale=1.1275;
    
    %bottom row fromright to left
    set(h(1), 'Position',[ pos(1,1)+hmove2  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
    
    set(h(2), 'Position',[ pos(2,1)+hmove3  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
    set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove2  pos(3,3)*hscale pos(3,4)*vscale])
    
    set(h(4), 'Position',[ pos(4,1)+hmove4  pos(4,2)+vmove1  pos(4,3)*hscale pos(4,4)*vscale])
    set(h(5), 'Position', [ pos(5,1)+hmove4  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
    set(h(6), 'Position', [ pos(6,1)+hmove4  pos(6,2)+vmove3  pos(6,3)*hscale pos(6,4)*vscale])
    
    set(h(7), 'Position', [ pos(7,1)+hmove1  pos(7,2)+vmove2  pos(7,3)*hscale pos(7,4)*vscale])
    
    set(h(8), 'Position', [ pos(8,1)+hmove1  pos(8,2)+vmove3  pos(8,3)*hscale pos(8,4)*vscale])
    set(h(9),  'Position',[ pos(9,1)+hmove2  pos(9,2)+vmove3  pos(9,3)*hscale  pos(9,4)*vscale])
    
    set(h(10), 'Position',[ pos(10,1)+hmove1 pos(10,2)+vmove4 pos(10,3)*hscale pos(10,4)*vscale])
    set(h(11), 'Position',[ pos(11,1)+hmove2 pos(11,2)+vmove4 pos(11,3)*hscale pos(11,4)*vscale]) 
    set(h(12), 'Position',[ pos(12,1)+hmove3 pos(12,2)+vmove4 pos(12,3)*hscale pos(12,4)*vscale])
    
    % Text
    axes('Units','Normalized','position',[0 0 1 1],'visible','off');
    text(0.500,0.045,'Q_{20yr} [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
    text(0.960,0.550,'7Q10 [mm/day]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);   
    text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
    %colorbar, text on the plots
    axes('Units','Normalized','position',[0 0 1 1],'visible','off');
    caxis([0 3]);
    colorbar('horizontal','position',[0.775,0.040,0.100,0.0195],...
        'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
    cbfreeze;
    
    %save figure
    if prnt
        figfile=['./figure/paper/Paper_m7_fig4_extreme_comparison_force_'  model{m} '_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
        %         set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-r300',figfile);
    end
end

toc