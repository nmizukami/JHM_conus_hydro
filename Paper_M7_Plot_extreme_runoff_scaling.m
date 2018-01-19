%{
Plot westen US wide HUC analysis for runoff stat
-monthly mean
-monthly variance
-annual mean
-annual variance
-annual runoff centroid
-
%}
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
% -------------------------------------------------------------------------
% CLM output to be read
% -------------------------------------------------------------------------
roclmvar{1,1}='ro_mean';   % mean annual runoff 
roclmvar{2,1}='ro_var';    
roclmvar{3,1}='ro_02yr_max';
roclmvar{4,1}='ro_05yr_max';
roclmvar{5,1}='ro_10yr_max';
roclmvar{6,1}='ro_20yr_max';
roclmvar{7,1}='ro_50yr_max';
roclmvar{8,1}='ro5_02yr_max';
roclmvar{9,1}='ro5_05yr_max';
roclmvar{10,1}='ro5_10yr_max';
roclmvar{11,1}='ro5_20yr_max';
roclmvar{12,1}='ro5_50yr_max';
roclmvar{13,1}='ro_10yr_min';
roclmvar{14,1}='ro7_10yr_min';
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
nxlon = length(lon1d);
nylat = length(lat1d);
clear us_domain_grid

%Size of dimensions
nforce  = length(force);
nmodel  = length(model);
nregion = length(region);
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
% -------------------------------------------------------------------------
% Some Plot properties
% -------------------------------------------------------------------------
%Import color map for surface map
mycolors  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/rainbow.gp',3,10);
mycolors  = (1/255)*mycolors;
mycolors1  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/BlRe.rgb',3,2);
mycolors1  = (1/255)*mycolors1;
mycolors2  = [255 0 255
              0 0 255
              0 255 0
              0 100 0
              255 102 0
              173 216 230
              0 0 0
              255 0 0];
mycolors2  = (1/255)*mycolors2;         
% line color order
clr4=[0 1 1;1 0 0;0 1 0;0 0 1];
clr5=[0 0 1;0.65 0.16 0.16;1 0 0;0 1 0;0 0 0];
%% Get data - CLM I/O and perform runoff process
% memory allocation
%HUC2 ro50max_us [huc x x prod x model ]
id_huc2_us           = ones(nregion,nforce,nmodel)*NaN;
ele_huc2_us          = ones(nregion,nforce,nmodel)*NaN;
RO_mean_huc2_us      = ones(nregion,nforce,nmodel)*NaN;
RO_var_huc2_us       = ones(nregion,nforce,nmodel)*NaN;
ro_50yr_max_huc2_us  = ones(nregion,nforce,nmodel)*NaN;
ro_20yr_max_huc2_us  = ones(nregion,nforce,nmodel)*NaN;
ro_10yr_max_huc2_us  = ones(nregion,nforce,nmodel)*NaN;
ro_5yr_max_huc2_us   = ones(nregion,nforce,nmodel)*NaN;
ro_2yr_max_huc2_us   = ones(nregion,nforce,nmodel)*NaN;
ro7_10yr_min_huc2_us = ones(nregion,nforce,nmodel)*NaN;
ro_10_min_huc2_us    = ones(nregion,nforce,nmodel)*NaN ;
%HUC4
id_huc4_us           = [];
ele_huc4_us          = [];
RO_mean_huc4_us      = [];
RO_var_huc4_us       = [];
ro_50yr_max_huc4_us  = ones(2000,nforce,nmodel)*NaN;
ro_20yr_max_huc4_us  = [];
ro_10yr_max_huc4_us  = [];
ro_5yr_max_huc4_us   = [];
ro_2yr_max_huc4_us   = [];
ro7_10yr_min_huc4_us = ones(2000,nforce,nmodel)*NaN;
ro_10yr_min_huc4_us  = [];
%HUC8
id_huc8_us           = [];
ele_huc8_us          = [];
RO_mean_huc8_us      = [];
RO_var_huc8_us       = [];
ro_50yr_max_huc8_us  = ones(30000,nforce,nmodel)*NaN;
ro_20yr_max_huc8_us  = [];
ro_10yr_max_huc8_us  = [];
ro_5yr_max_huc8_us   = [];
ro_2yr_max_huc8_us   = [];
ro7_10yr_min_huc8_us = ones(30000,nforce,nmodel)*NaN;
ro_10yr_min_huc8_us  = [];

% CONUS [lat x lon x prod x model]
ro50max_us  = zeros(nylat,nxlon,nforce,nmodel);
ro710min_us = zeros(nylat,nxlon,nforce,nmodel);

%Create simulation area mask
sim_mask = repmat(lon2d,[1,1,nforce,nmodel]);
sim_mask(~isnan(sim_mask)) = 1;

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
    
    %Reading gridded runoff stats    
    % -------------------------------------------------------------------------
    % Read region domain grid
    % -------------------------------------------------------------------------
    region_domain_grid = ['/d3/mizukami/domain_huc/domain_' region{r} '_' num2str(res) 'k.nc'];
    lon_region_1d = netcdf2mat(region_domain_grid,'lon');
    lat_region_1d = netcdf2mat(region_domain_grid,'lat');
    nxlon_region = length(lon_region_1d);
    nylat_region = length(lat_region_1d);
    clear region_domain_grid
    
    % get indices for region domain from western domain array
    i1 = find( lat1d==lat_region_1d(1) );
    i2 = find( lat1d==lat_region_1d(end) );
    j1 = find( lon1d==lon_region_1d(1) );
    j2 = find( lon1d==lon_region_1d(end) );
    % -------------------------------------------------------------------------
    %  Reading NetCDF
    % -------------------------------------------------------------------------
    for p=1:nforce
        for m = 1:nmodel
            if strcmp(model{m},'CLM')
                main_indir ='/d3/mizukami/CLM_OUTPUT';
            elseif strcmp (model{m},'VIC')
                main_indir ='/d3/mizukami/VIC_OUTPUT/netcdf/processed';
            elseif strcmp (model{m},'PRMS')
                main_indir ='/d3/mizukami/PRMS_OUTPUT';                
            end
            ncname = [main_indir '/' force{p} '_'  region{r} '/Extreme_RUNOFF.nc'];
            [var0,fillvalue]=netcdf2mat(ncname,roclmvar{7},'attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro50max_us(i1:i2,j1:j2,p,m) = var0+ro50max_us(i1:i2,j1:j2,p,m);
            clear var0 
            [var0,fillvalue]=netcdf2mat(ncname,roclmvar{14},'attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro710min_us(i1:i2,j1:j2,p,m) = var0+ro710min_us(i1:i2,j1:j2,p,m);        
            clear var0 ncname
        end
    end
    ro50max_us = ro50max_us.*sim_mask;
    ro710min_us = ro710min_us.*sim_mask;
    %Reading huc runoff stats
    % -------------------------------------------------------------------------
    %  Reading mat file
    % -------------------------------------------------------------------------
    %  mat file contains
    %     ID_huc2                nx1
    %     elev_huc2              nx1
    %     ro5_02yr_max_huc2      nx1
    %     ro5_05yr_max_huc2      nx1
    %     ro5_10yr_max_huc2      nx1
    %     ro5_20yr_max_huc2      nx1
    %     ro5_50yr_max_huc2      nx1
    %     ro7_10yr_min_huc2      nx1
    %     ro_02yr_max_huc2       nx1
    %     ro_05yr_max_huc2       nx1
    %     ro_10yr_max_huc2       nx1
    %     ro_10yr_min_huc2       nx1
    %     ro_20yr_max_huc2       nx1
    %     ro_50yr_max_huc2       nx1
    %     ro_mean_huc2           nx1
    %     ro_var_huc2            nx1
    if r==1
        i21 = 1;
        i41 = 1;
        i81 = 1;
    else
        i21=i22+1;
        i41=i42+1;
        i81=i82+1;
    end
    for p=1:nforce
        for m = 1:nmodel
            if strcmp(model{m},'CLM')
                main_indir ='/d3/mizukami/CLM_OUTPUT';
            elseif strcmp (model{m},'VIC')
                main_indir ='/d3/mizukami/VIC_OUTPUT/netcdf/processed';
            elseif strcmp (model{m},'PRMS')
                main_indir ='/d3/mizukami/PRMS_OUTPUT';                
            end
            filename_huc2 = [main_indir '/' force{p} '_'  region{r} '/Extreme_RUNOFF_huc2.mat'];
            load(filename_huc2)           
            filename_huc4 = [main_indir '/' force{p} '_'  region{r} '/Extreme_RUNOFF_huc4.mat'];            
            load(filename_huc4)           
            filename_huc8 = [main_indir '/' force{p} '_'  region{r} '/Extreme_RUNOFF_huc8.mat'];
            load(filename_huc8)
            i22 = i21+length(ro_50yr_max_huc2)-1;
            i42 = i41+length(ro_50yr_max_huc4)-1;
            i82 = i81+length(ro_50yr_max_huc8)-1;
            ro_50yr_max_huc2_us(i21:i22,p,m)  = ro_50yr_max_huc2;
            ro_50yr_max_huc4_us(i41:i42,p,m)  = ro_50yr_max_huc4;
            ro_50yr_max_huc8_us(i81:i82,p,m)  = ro_50yr_max_huc8;

            ro7_10yr_min_huc2_us(i21:i22,p,m)  = ro7_10yr_min_huc2;
            ro7_10yr_min_huc4_us(i41:i42,p,m)  = ro7_10yr_min_huc4;
            ro7_10yr_min_huc8_us(i81:i82,p,m)  = ro7_10yr_min_huc8;            
        end
    end
 
end
ro50max_us(ro50max_us>100)=NaN;
ro710min_us(ro710min_us>100)=NaN;
% -------------------------------------------------------------------------
% Compute Statistics [scale x force x model]
% -------------------------------------------------------------------------
mean_ro_max50  = ones(4,nforce,nmodel)*NaN;
mean_ro7_min10 = ones(4,nforce,nmodel)*NaN;
for p=1:nforce
    for m = 1:nmodel
        mean_ro_max50(:,p,m)  = [squeeze(nanmean(nanmean(ro50max_us(:,:,p,m),1),2))
            squeeze(nanmean(ro_50yr_max_huc8_us(:,p,m),1))
            squeeze(nanmean(ro_50yr_max_huc4_us(:,p,m),1))
            squeeze(nanmean(ro_50yr_max_huc2_us(:,p,m),1))];
        mean_ro7_min10(:,p,m)  = [squeeze(nanmean(nanmean(ro710min_us(:,:,p,m),1),2))
            squeeze(nanmean(ro7_10yr_min_huc8_us(:,p,m),1))
            squeeze(nanmean(ro7_10yr_min_huc4_us(:,p,m),1))
            squeeze(nanmean(ro7_10yr_min_huc2_us(:,p,m),1))];        
    end
end
%% Plot 1 - 50yr daily maximum runoff
for m=1:nmodel
    figure('Name','1 - Scaling of extreme runoff','NumberTitle','off','Units','inches','position',[3 7 5.75 3.5],'color',[1 1 1]);
    set(gcf,'DefaultAxesColorOrder',clr5)
    subplot(1,2,1)
    plot(mean_ro_max50(:,1:4,m),'-','LineWidth',2);hold on
    plot(mean_ro_max50(:,5,m),'-o','LineWidth',2,'color',clr5(5,:))
    title('50yr Daily max. flow')
    set(gca,'XLim',[0.5 4.5]);
    if m<3; set(gca,'YLim',[0 25]); else set(gca,'YLim',[0 1000]); end
    set(gca,'XTick',1:1:4);
    set(gca,'XTickLabel',{'Full Res.';'HUC8';'HUC4';'HUC2'});
    ylabel('mm/day 50yr return');
    set(gca,'FontSize',9)
    
    subplot(1,2,2)
    plot(mean_ro7_min10(:,1:4,m),'-','LineWidth',2);hold on
    plot(mean_ro7_min10(:,5,m),'-o','LineWidth',2,'color',clr5(5,:))
    title('10yr 7-day min. flow')
    set(gca,'XLim',[0.5 4.5]);
    if m<3; set(gca,'YLim',[0 0.08]); else set(gca,'YLim',[0 0.20]); end
    set(gca,'XTick',1:1:4);
    set(gca,'XTickLabel',{'Full Res.';'HUC8';'HUC4';'HUC2'});
    ylabel('mm/day 10yr return');
    set(gca,'FontSize',9)
    
    %adjust plot size & location
    h   = get(gcf,'Children');
    for i=1:length(h)
        pos(i,:) =get(h(i),'position');
    end
    hmove1=+0.025;
    hmove2=-0.025;
    vmove=-0.01;
    hscale=1.10;
    vscale=1.00;
    %bottom row fromright to left
    set(h(1), 'Position',[ pos(1,1)+hmove1 pos(1,2)+vmove pos(1,3)*hscale pos(1,4)*vscale])
    set(h(2), 'Position',[ pos(2,1)+hmove2 pos(2,2)+vmove pos(2,3)*hscale pos(2,4)*vscale])
    
    %legends 2
    axes('unit','Normalized','position',[0 0 1 1],'visible','off');
    dummy1 = plot(1,1,'-','LineWidth',2,'color',clr5(3,:));hold on
    dummy2 = plot(1,1,'-','LineWidth',2,'color',clr5(2,:));
    dummy3 = plot(1,1,'-','LineWidth',2,'color',clr5(1,:));
    dummy4 = plot(1,1,'-','LineWidth',2,'color',clr5(4,:));
    dummy5 = plot(1,1,'-o','LineWidth',2,'color',clr5(5,:));
    set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off');set(dummy5,'Visible','off');
    set(gca,'Visible','off','FontSize',9)
    legend('BCSDm','BCSDd','BCCA','AR','M02',...
        'Orientation','vertical','Location',[0.655 0.765 0.075 0.025]);
    legend('boxon','Color','none');
    
    %save figure
    figfile=['./figure/paper/Paper_m7_fig1-extreme_scaling_WYR' num2str(wyr(1)) '-' num2str(wyr(end)) '_' model{m} '_conus_' num2str(res) 'K'];
    set(gcf,'PaperPositionMode','auto')
    print('-dtiff',figfile);
end
%% Plot 5 - 7Q10 runoff
figure('Name','5 - 7Q10','NumberTitle','off','Units','inches','position',[3 7 6.75 6.75],'Visible','on');
set(gcf,'DefaultAxesColorOrder',clr5)
for r = 1:length(region)
    if r <= 4
        subplot(3,3,r)
    elseif r >= 5
        subplot(3,3,r+1)
    end  
    plot(mean_ro7_min10(:,:,r)','-o','LineWidth',2)
    title(region{r})
    set(gca,'XLim',[0.5 4.5]);
    set(gca,'XTick',1:1:4);
    set(gca,'XTickLabel',{'full res';'HUC8';'HUC4';'HUC2'});
    if r==7
        set(gca,'YLim',[0 5*10^-4]);
    end
    rotateXLabels(gca, 30)
end
subplot(3,3,5)
ax2 = usamap([lat_west1d(1) lat_west1d(end)],[lon_west1d(1)-360 lon_west1d(end)-360]);
set(ax2, 'Visible','off','layer','top');
setm(ax2,'MapProjection','mercator','Grid','off','FontSize',9,'Frame','off');
plotm([S.Y],[S.X], 'k');
tightmap;
mlabel('off'); plabel('off');

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hspace=0.10;
vspace=0.05;
hscale=1.12;
vscale=1.12;
%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)-hspace*0.5  pos(1,2)-vspace*1.4  pos(1,3)*1.5 pos(1,4)*1.5]) 
set(h(2), 'Position',[ pos(2,1)+hspace*0.5  pos(2,2)-vspace*0.9  pos(2,3)*hscale pos(2,4)*vscale]) %RIO
set(h(3), 'Position',[ pos(3,1)-hspace*0.0  pos(3,2)-vspace*0.9  pos(3,3)*hscale pos(3,4)*vscale]) %LCO
set(h(4), 'Position',[ pos(4,1)-hspace*0.5  pos(4,2)-vspace*0.9  pos(4,3)*hscale pos(4,4)*vscale]) %CA
set(h(5), 'Position',[ pos(5,1)+hspace*0.5  pos(5,2)-vspace*0.3  pos(5,3)*hscale pos(5,4)*vscale]) %AR
set(h(6), 'Position',[ pos(6,1)-hspace*0.5  pos(6,2)-vspace*0.3  pos(6,3)*hscale pos(6,4)*vscale]) %GB
set(h(7), 'Position',[ pos(7,1)+hspace*0.5  pos(7,2)+vspace*0.35  pos(7,3)*hscale pos(7,4)*vscale]) %MR
set(h(8), 'Position',[ pos(8,1)-hspace*0.0  pos(8,2)+vspace*0.33  pos(8,3)*hscale pos(8,4)*vscale]) %UCO
set(h(9), 'Position',[ pos(9,1)-hspace*0.5  pos(9,2)+vspace*0.35  pos(9,3)*hscale pos(9,4)*vscale]) %PN

%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-o','LineWidth',2,'color',clr5(1,:));hold on
dummy2 = plot(1,1,'-o','LineWidth',2,'color',clr5(2,:));
dummy3 = plot(1,1,'-o','LineWidth',2,'color',clr5(3,:));
dummy4 = plot(1,1,'-o','LineWidth',2,'color',clr5(4,:));
dummy5 = plot(1,1,'-o','LineWidth',2,'color',clr5(5,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off');set(dummy5,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('BCCA','BCSDd','BCSDm','AR','Maurer',...
    'Orientation','vertical','Location',[0.475 0.225 0.050 0.025]);
legend('boxoff','Color','white');

%save figure
figfile=['./figure/m12_fig6-RO7_MIN10_WYR' num2str(wyr(1)) '-' num2str(wyr(end)) '_west_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dtiff',figfile);  

toc
% text(0.5,0.975,'10yr daily minimum runoff [mm/day]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);
