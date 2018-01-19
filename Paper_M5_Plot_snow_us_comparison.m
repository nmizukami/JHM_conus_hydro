% Comparison of 3month SW, LW, RH climatology between Maurer vs 4 SDs
%          plot_type   temp_res.     climatology  spatial res.   variables
%  -------------------------------------------------------------------------------
%  Fig 1-1 CDF-bias    annual        yes          grid           Peak SWE  

tic
profile on
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
% 
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
elseif strcmp (model,'PRMS')
    main_indir ='/d3/mizukami/PRMS_OUTPUT';
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
us_huc2_shp = '/d3/mizukami/domain_huc/shapefile/HUC02_conus.shp';
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
% CONUS [lat x lon x 365 x prod]
swe_us = zeros(nylat,nxlon,365,nforce,nmodel);

%Create simulation area mask
sim_mask = lon2d+huc2g;
sim_mask(~isnan(sim_mask) ) = 1;
sim_mask = repmat(sim_mask,[1,1,365,nforce,nmodel]);

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
  %showrtwave radiation
  for f=1:nforce
      for m = 1:nmodel 
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
          ncname = [main_indir '/' force{f} '_'  region{r} '/calib/Daily_clim_SWE.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'SWE','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          swe_us(i1:i2,j1:j2,:,f,m) = var0+swe_us(i1:i2,j1:j2,:,f,m);
          clear var0 ncname
      end
  end
end
swe_us = swe_us.*sim_mask;

%Days from October 1st
daynum = linspace(1,365,365);
days = repmat(reshape(daynum,1,1,length(daynum)),[nylat,nxlon,1,length(force),length(model)]);

%Compute stat - peak SWE and center of mass
swe_peak = squeeze( max( swe_us,[],3) );
swe_ct   = squeeze( sum(swe_us.*days,3)./sum(swe_us,3) );
 
swe_peak_mean=mean(mean(swe_peak,4),3);
swe_peak_mean(swe_peak_mean>5000)=NaN;
%Eliminate low swe accumulation area
[rows,cols]=find(swe_peak_mean<50);
for r=1:length(rows)
    swe_peak(rows(r),cols(r),:,:)=NaN;
    swe_ct(rows(r),cols(r),:,:)=NaN;
end
% -------------------------------------------------------------------------
% Make plot 
% -------------------------------------------------------------------------
%% Compute CDF of inter-model and inter-forcing difference
%1. Inter model differece for peak swe and swe ct
% Order
%1 CLM-VIC
%2 CLM-PRMS
%3 VIC-PRMS
%memory allocation
dif_swe_peak_model = ones(nylat*nxlon,nforce,nchoosek(nmodel,2))*NaN;
dif_swe_ct_model = ones(nylat*nxlon,nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            dif_swe_peak_model(:,f,c) = reshape( swe_peak(:,:,f,m1)-swe_peak(:,:,f,m2), nylat*nxlon,1);
            dif_swe_ct_model(:,f,c)   = reshape( swe_ct(:,:,f,m1)-swe_ct(:,:,f,m2), nylat*nxlon,1);
        end
    end
end
clear c m1 m2
%1.2 Inter model differece for peak swe and swe ct
% Order
%1 BCSDd-BCCA
%2 BCSDm-BCCA
%3 AR-BCCA
%4 BCSDm-BCSDd
%5 AR-BCSDd
%6 AR-BCSDm
%memory allocation
dif_swe_peak_force = ones(nylat*nxlon,nchoosek(nforce,2),nmodel)*NaN;
dif_swe_ct_force   = ones(nylat*nxlon,nchoosek(nforce,2),nmodel)*NaN;
for m = 1:nmodel
    c=0;
    for f1 = 1:nforce-1
        for f2 = f1+1:nforce
            c=c+1;
            dif_swe_peak_force(:,c,m) = reshape( swe_peak(:,:,f2,m)-swe_peak(:,:,f1,m), nylat*nxlon,1);
            dif_swe_ct_force(:,c,m)   = reshape( swe_ct(:,:,f2,m)-swe_ct(:,:,f1,m), nylat*nxlon,1);
        end
    end
end
clear c f1 f2
dif_swe_peak_model(dif_swe_peak_model<-500 | dif_swe_peak_model>500)=NaN;
dif_swe_peak_force(dif_swe_peak_force<-500 | dif_swe_peak_force>500)=NaN;
dif_swe_ct_model(dif_swe_ct_model<-200 | dif_swe_ct_model>200)=NaN;
dif_swe_ct_force(dif_swe_ct_force<-200 | dif_swe_ct_force>200)=NaN;

%2   compute cdf for inter model difference in peak SWE
f_peak_model=cell(nforce,nchoosek(nmodel,2));
x_peak_model=cell(nforce,nchoosek(nmodel,2));
for f = 1:nforce
    for m=1:nchoosek(nmodel,2)
        [cdf,x] =ecdf(dif_swe_peak_model(:,f,m));
        f_peak_model{f,m}=cdf;
        x_peak_model{f,m}=x;
    end
end
%3   compute cdf for inter forcing difference in peak SWE
f_peak_force=cell(nchoosek(nforce,2),nmodel);
x_peak_force=cell(nchoosek(nforce,2),nmodel);
for f = 1:nchoosek(nforce,2)
    for m=1:nmodel
        [cdf,x] =ecdf(dif_swe_peak_force(:,f,m));
        f_peak_force{f,m}=cdf;
        x_peak_force{f,m}=x;
    end
end
%4   compute cdf for inter model difference in SWE CT
f_ct_model=cell(nforce,nchoosek(nmodel,2));
x_ct_model=cell(nforce,nchoosek(nmodel,2));
for f = 1:nforce
    for m=1:nchoosek(nmodel,2)
        [cdf,x] =ecdf(dif_swe_ct_model(:,f,m));
        f_ct_model{f,m}=cdf;
        x_ct_model{f,m}=x;
    end
end
%5   compute cdf for inter forcing difference in SWE CT
f_ct_force=cell(nchoosek(nforce,2),nmodel);
x_ct_force=cell(nchoosek(nforce,2),nmodel);
for f = 1:nchoosek(nmodel,2)
    for m=1:nmodel
        [cdf,x] =ecdf(dif_swe_ct_force(:,f,m));
        f_ct_force{f,m}=cdf;
        x_ct_force{f,m}=x;
    end
end

%Bias
bias_swe_peak = ones(nylat*nxlon,4,nmodel)*NaN;
bias_swe_ct = ones(nylat*nxlon,4,nmodel)*NaN;
for m = 1:nmodel
    bias_swe_peak(:,1,m) = reshape( swe_peak(:,:,1,m)-swe_peak(:,:,5,m), nylat*nxlon,1);
    bias_swe_peak(:,2,m) = reshape( swe_peak(:,:,2,m)-swe_peak(:,:,5,m), nylat*nxlon,1);
    bias_swe_peak(:,3,m) = reshape( swe_peak(:,:,3,m)-swe_peak(:,:,5,m), nylat*nxlon,1);
    bias_swe_peak(:,4,m) = reshape( swe_peak(:,:,4,m)-swe_peak(:,:,5,m), nylat*nxlon,1);
    bias_swe_ct(:,1,m) = reshape( swe_ct(:,:,1,m)-swe_ct(:,:,5,m), nylat*nxlon,1);
    bias_swe_ct(:,2,m) = reshape( swe_ct(:,:,2,m)-swe_ct(:,:,5,m), nylat*nxlon,1);
    bias_swe_ct(:,3,m) = reshape( swe_ct(:,:,3,m)-swe_ct(:,:,5,m), nylat*nxlon,1);
    bias_swe_ct(:,4,m) = reshape( swe_ct(:,:,4,m)-swe_ct(:,:,5,m), nylat*nxlon,1);
end
bias_swe_peak(bias_swe_peak<-500 | bias_swe_peak>500)=NaN;
bias_swe_ct(bias_swe_ct<-200 | bias_swe_ct>200)=NaN;

%compute cdf
f_peak_bias=cell(nforce-1,nmodel);
x_peak_bias=cell(nforce-1,nmodel);
for f = 1:nforce-1
    for m=1:nmodel
        [cdf,x] =ecdf(bias_swe_peak(:,f,m));
        f_peak_bias{f,m}=cdf;
        x_peak_bias{f,m}=x;
    end
end
f_ct_bias=cell(nforce-1,nmodel);
x_ct_bias=cell(nforce-1,nmodel);
for f = 1:nforce-1
    for m=1:nmodel
        [cdf,x] =ecdf(bias_swe_ct(:,f,m));
        f_ct_bias{f,m}=cdf;
        x_ct_bias{f,m}=x;
    end
end
%% plot stuf
%Import color map for surface map
%mycolors  = colortable('/glade/u/home/mizukami/mtl/myColortables/WhBlGrYeRe.gp',3,10);
mycolors  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/rainbow.gp',3,10);
mycolors  = (1/255).*mycolors;
mycolors1  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/BlRe.rgb',3,2);
mycolors1  = (1/255).*mycolors1;
mycolors2  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/temp_diff_18lev.rgb',3,6);
mycolors2  = (1/255).*mycolors2;
% mycolors2(11,:)=[];
mycolors3  = colortable('/home/mizukami/hydro_nm/mtl/myColortables/MPL_Blues.rgb',3,2);

% line color order
clr3=[0 0 1;1 0 0;0 1 0];
clr5=[0 0 1;0.65 0.16 0.16;1 0 0;0 1 0;0 0 0];

force_name{1,1}='BCCA';
force_name{2,1}='BCSDd';
force_name{3,1}='BCSDm';
force_name{4,1}='AR';
force_name{5,1}='M02';
%% plot 1 - CDF of inter-model SWE comparison
figure('Units','inches','position',[1 4 7.25 4.50],'Visible','on');
%Plots for peak SWE
for f=1:5
    subplot(2,5,f)
    plot(x_peak_model{f,1},f_peak_model{f,1},'LineWidth',1,'Color',clr3(1,:));hold on   
    plot(x_peak_model{f,2},f_peak_model{f,2},'LineWidth',1,'Color',clr3(2,:));
    plot(x_peak_model{f,3},f_peak_model{f,3},'LineWidth',1,'Color',clr3(3,:));
    grid on;
    set(gca,'FontSize',8);
    set(gca,'XLim',[-120 120]);
    set(gca,'YLim',[0 1]);
    set(gca,'XTick',-100:50:100);
    set(gca,'YTick',0:0.25:1);
    if f==1
        ylabel('CDF');
    else
        set(gca,'YTickLabel',{''});
    end
    axis square
    if f == 1;title('BCCA');
    elseif f == 2;title('BCSDd');
    elseif f == 3;title('BCSDm');    
    elseif f == 4;title('AR');
    elseif f == 5;title('M02');    
    end
end
%Plots for peak SWE
for f=1:5
    subplot(2,5,f+5)
    plot(x_ct_model{f,1},f_ct_model{f,1},'LineWidth',1,'Color',clr3(1,:));hold on   
    plot(x_ct_model{f,2},f_ct_model{f,2},'LineWidth',1,'Color',clr3(2,:));
    plot(x_ct_model{f,3},f_ct_model{f,3},'LineWidth',1,'Color',clr3(3,:));  
    grid on;
    set(gca,'FontSize',8);
    set(gca,'XLim',[-35 35]);
    set(gca,'YLim',[0 1]);
    set(gca,'XTick',-30:15:30);
    set(gca,'YTick',0:0.25:1); 
    if f==1
        ylabel('CDF');
    else
        set(gca,'YTickLabel',{''});
    end
    axis square
    if f == 1;title('BCCA');
    elseif f == 2;title('BCSDd');
    elseif f == 3;title('BCSDm');    
    elseif f == 4;title('AR');
    elseif f == 5;title('M02');    
    end
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hspace1=+0.0200;
hspace2=+0.0025;
hspace3=-0.0150;
hspace4=-0.0325;
hspace5=-0.0500;
vspace1=-0.02;
vspace2=-0.05;
hscale=1.350;
vscale=1.350;
% %bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hspace1  pos(1,2)+vspace1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hspace2  pos(2,2)+vspace1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hspace3  pos(3,2)+vspace1  pos(3,3)*hscale pos(3,4)*vscale])
set(h(4), 'Position',[ pos(4,1)+hspace4  pos(4,2)+vspace1  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hspace5  pos(5,2)+vspace1  pos(5,3)*hscale pos(5,4)*vscale])

set(h(6), 'Position',[ pos(6,1)+hspace1  pos(6,2)+vspace2  pos(6,3)*hscale pos(6,4)*vscale])
set(h(7), 'Position',[ pos(7,1)+hspace2  pos(7,2)+vspace2  pos(7,3)*hscale pos(7,4)*vscale])
set(h(8), 'Position',[ pos(8,1)+hspace3  pos(8,2)+vspace2  pos(8,3)*hscale pos(8,4)*vscale])
set(h(9), 'Position',[ pos(9,1)+hspace4  pos(9,2)+vspace2  pos(9,3)*hscale pos(9,4)*vscale])
set(h(10),'Position',[ pos(10,1)+hspace5 pos(10,2)+vspace2 pos(10,3)*hscale pos(10,4)*vscale])

axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.525,0.970,'Inter-model difference in peak SWE [mm]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0)
text(0.525,0.525,'Inter-model difference in SWE CT [day]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0)
%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',1,'color',clr3(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',1,'color',clr3(2,:));
dummy3 = plot(1,1,'-','LineWidth',1,'color',clr3(3,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('CLM-VIC','CLM-PRMS','VIC-PRMS',...
    'Orientation','horizontal','Location',[0.450 0.050 0.050 0.025]);
legend('boxoff','Color','white');

%save figure
figfile=['./figure/paper/Paper_m5_fig1_swe_comparison_model' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);
%% plot 2 - CDF of bias SWE comparison
figure('Units','inches','position',[1 4 6.55 4.50],'Visible','on');
%Plots for peak SWE
for m=1:3
    subplot(2,3,m)
    plot(x_peak_bias{1,m},f_peak_bias{1,m},'LineWidth',1,'Color',clr5(1,:));hold on   
    plot(x_peak_bias{2,m},f_peak_bias{2,m},'LineWidth',1,'Color',clr5(2,:));
    plot(x_peak_bias{3,m},f_peak_bias{3,m},'LineWidth',1,'Color',clr5(3,:));
    plot(x_peak_bias{4,m},f_peak_bias{4,m},'LineWidth',1,'Color',clr5(4,:)); 
    grid on;
    set(gca,'FontSize',8);
    set(gca,'XLim',[-120 120]);
    set(gca,'YLim',[0 1]);
    set(gca,'XTick',-100:50:100);
    set(gca,'YTick',0:0.25:1);
    if m==1
        ylabel('CDF');
    else
        set(gca,'YTickLabel',{''});
    end
    axis square
    if m == 1;title('CLM');
    elseif m == 2;title('VIC');
    elseif m == 3;title('PRMS');       
    end
end
%Plots for peak SWE
for m=1:3
    subplot(2,3,m+3)
    plot(x_ct_bias{1,m},f_ct_bias{1,m},'LineWidth',1,'Color',clr5(1,:));hold on   
    plot(x_ct_bias{2,m},f_ct_bias{2,m},'LineWidth',1,'Color',clr5(2,:));
    plot(x_ct_bias{3,m},f_ct_bias{3,m},'LineWidth',1,'Color',clr5(3,:));
    plot(x_ct_bias{4,m},f_ct_bias{4,m},'LineWidth',1,'Color',clr5(4,:));  
    grid on;
    set(gca,'FontSize',8);
    set(gca,'XLim',[-35 35]);
    set(gca,'YLim',[0 1]);
    set(gca,'XTick',-30:15:30);
    set(gca,'YTick',0:0.25:1); 
    if m==1
        ylabel('CDF');
    else
        set(gca,'YTickLabel',{''});
    end
    axis square
    if m == 1;title('CLM');
    elseif m == 2;title('VIC');
    elseif m == 3;title('PRMS');    
    end
end

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hspace1=+0.0200;
hspace2=-0.0150;
hspace3=-0.0500;
vspace1=-0.00;
vspace2=+0.00;
hscale=1.005;
vscale=1.005;
% %bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hspace1  pos(1,2)+vspace1  pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hspace2  pos(2,2)+vspace1  pos(2,3)*hscale pos(2,4)*vscale])
set(h(3), 'Position',[ pos(3,1)+hspace3  pos(3,2)+vspace1  pos(3,3)*hscale pos(3,4)*vscale])

set(h(4), 'Position',[ pos(4,1)+hspace1  pos(4,2)+vspace2  pos(4,3)*hscale pos(4,4)*vscale])
set(h(5), 'Position',[ pos(5,1)+hspace2  pos(5,2)+vspace2  pos(5,3)*hscale pos(5,4)*vscale])
set(h(6), 'Position',[ pos(6,1)+hspace3  pos(6,2)+vspace2  pos(6,3)*hscale pos(6,4)*vscale])

axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.525,0.970,'Bias in peak SWE [mm]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0)
text(0.525,0.500,'Bias in SWE CT [day]','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0)
%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',1,'color',clr5(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',1,'color',clr5(2,:));
dummy3 = plot(1,1,'-','LineWidth',1,'color',clr5(3,:));
dummy4 = plot(1,1,'-','LineWidth',1,'color',clr5(4,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('BCCA-M02','BCSDd-M02','BCSDm-M02','AR-M02',...
    'Orientation','horizontal','Location',[0.450 0.035 0.050 0.025]);
legend('boxoff','Color','white');

%save figure
figfile=['./figure/paper/Paper_m5_fig2_swe_comparison_bias' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);

toc
