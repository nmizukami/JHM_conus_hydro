% Inter-model comparison of annual water balance climatology (RR, ET and
% RO) over CONUS.  Climatology is period beteween WY1981 through WY2008.
% 
% Required netCDF
% Model outputs
% - Annual_climate_RUNOFF.nc
% - Annual_climate_ET.nc
% - Annual_climate_RR.nc
% Ancillary files
% - domain_<region>_<res>k.nc
%
%          plot_type     temp_res.  time span   climatology  spatial res.    variables
%  ----------------------------------------------------------------------------------------
%  Fig1    scatter       daily      annual      yes          grid box        ET, RO, RR     CLM vs. VIC
%  Fig2    scatter       daily      annual      yes          grid box        ET, RO, RR     CLM vs. PRMS

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
region{11,1} = 'GUL';
region{12,1} = 'RIO';
region{13,1} = 'UCO';
region{14,1} = 'LCO';
region{15,1} = 'PN';
region{16,1} = 'CA';
region{17,1} = 'GB';
region{18,1} = 'AR';
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
% eyr=2008; emon=9; eday=30;
eyr=1999; emon=9; eday=30;
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
% water year array 
% -------------------------------------------------------------------------
wyr = byr+1:1:eyr;
%% Getting daily clm datm variable from various product files
% -------------------------------------------------------------------------
% Memory allocation
% -------------------------------------------------------------------------
% West [lat x lon x region x prod x model]
pr_us = zeros(nylat,nxlon,nforce,nmodel);
et_us = zeros(nylat,nxlon,nforce,nmodel);
ro_us = zeros(nylat,nxlon,nforce,nmodel);
rr_us = zeros(nylat,nxlon,nforce,nmodel);

%Create simulation area mask
sim_mask = lon2d+huc2g;
sim_mask(~isnan(sim_mask) ) = 1;
sim_mask = repmat(sim_mask,[1,1,nforce,nmodel]);

for r = 1:16%nregion
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
    % Annual climatology value
    % [lat x lon x force x model]
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
            
            % Annual PRCP
            ncname = [main_indir '/' force{f} '_' region{r} '/calib/Annual_climate_PRCP.nc'];
            [var0,fillvalue]=netcdf2mat(ncname,'PRCP','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            pr_us(i1:i2,j1:j2,f,m) = var0*60*60*24*365+pr_us(i1:i2,j1:j2,f,m);
            clear var0 ncname
            
            % Annual RO
            ncname = [main_indir '/' force{f} '_' region{r} '/calib/Annual_climate_RUNOFF.nc'];
            [var0,fillvalue]=netcdf2mat(ncname,'RUNOFF','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            ro_us(i1:i2,j1:j2,f,m) = var0*60*60*24*365+ro_us(i1:i2,j1:j2,f,m);
            clear var0 ncname
            
            % Annual ET
            ncname = [main_indir '/' force{f} '_' region{r} '/calib/Annual_climate_ET.nc'];
            [var0,fillvalue]=netcdf2mat(ncname,'ET','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            et_us(i1:i2,j1:j2,f,m) = var0*60*60*24*365+et_us(i1:i2,j1:j2,f,m);
            clear var0 ncname
            
            % Annual Runoff ratio
            ncname = [main_indir '/' force{f} '_' region{r} '/calib/Annual_climate_RR.nc'];
            [var0,fillvalue]=netcdf2mat(ncname,'RR','attname1','_FillValue');
            var0(var0==fillvalue)=0;
            rr_us(i1:i2,j1:j2,f,m) = var0+rr_us(i1:i2,j1:j2,f,m);
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
pr_M02_max = squeeze( round(nanmax( nanmax( pr_us(:,:,5,:),[],1),[],2)) );
pr_M02_min = squeeze( round(nanmin( nanmin( pr_us(:,:,5,:),[],1),[],2)) );
et_M02_max = squeeze( round(nanmax( nanmax( et_us(:,:,5,:),[],1),[],2)) );
et_M02_min = squeeze( round(nanmin( nanmin( et_us(:,:,5,:),[],1),[],2)) );
ro_M02_max = squeeze( round(nanmax( nanmax( ro_us(:,:,5,:),[],1),[],2)) );
ro_M02_min = squeeze( round(nanmin( nanmin( ro_us(:,:,5,:),[],1),[],2)) );
rr_M02_max = squeeze( round(nanmax( nanmax( rr_us(:,:,5,:),[],1),[],2)) );
rr_M02_min = squeeze( round(nanmin( nanmin( rr_us(:,:,5,:),[],1),[],2)) );

%% joint histogram for inter-model water balance comparison
% hist3 will bin the data
%rr
xirr = 0:0.02:1.0;
yirr = xirr;
[Xrr,Yrr] = meshgrid(xirr,yirr);
%et
xiet = 0:50:2000;
yiet = xiet;
[Xet,Yet] = meshgrid(xiet,yiet);
%ro
xiro = 0:50:2000;
yiro = xiro;
[Xro,Yro] = meshgrid(xiro,yiro);

%rr
%memory allocation
rr_model = ones(nylat*nxlon,2,nforce,nchoosek(nmodel,2))*NaN;
Nrr_model   = ones(length(xirr),length(yirr),nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            rr_model(:,:,f,c)  = [reshape( rr_us(:,:,f,m2), nylat*nxlon,1) reshape( rr_us(:,:,f,m1), nylat*nxlon,1)];
            Nrr_model(:,:,f,c) = hist3(rr_model(:,:,f,c),{xirr yirr});
        end
    end
end
clear c m1 m2
% normalize the histogram data
dxrr = xirr(2)-xirr(1);
dyrr = yirr(2)-yirr(1);
Prr_model = Nrr_model./repmat(sum(sum(Nrr_model,1),2),[length(xirr),length(yirr),1]);

%et
%memory allocation
et_model = ones(nylat*nxlon,2,nforce,nchoosek(nmodel,2))*NaN;
Net_model = ones(length(xiet),length(yiet),nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            et_model(:,:,f,c) = [reshape( et_us(:,:,f,m2), nylat*nxlon,1) reshape( et_us(:,:,f,m1), nylat*nxlon,1)];
            Net_model(:,:,f,c)   = hist3(et_model(:,:,f,c),{xiet yiet});
        end
    end
end
clear c m1 m2
%ro
%memory allocation
ro_model = ones(nylat*nxlon,2,nforce,nchoosek(nmodel,2))*NaN;
Nro_model = ones(length(xiro),length(yiro),nforce,nchoosek(nmodel,2))*NaN;
for f = 1:nforce
    c=0;
    for m1 = 1:nmodel-1
        for m2 = m1+1:nmodel
            c=c+1;
            ro_model(:,:,f,c) = [reshape( ro_us(:,:,f,m2), nylat*nxlon,1) reshape( ro_us(:,:,f,m1), nylat*nxlon,1)];
            Nro_model(:,:,f,c)   = hist3(ro_model(:,:,f,c),{xiro yiro});
        end
    end
end
clear c m1 m2
%% joint histogram for inter-forcing water balance comparison
%memory allocation
rr_force  = ones(nylat*nxlon,2,nchoosek(nforce,2),nmodel)*NaN;
Nrr_force = ones(length(xirr),length(yirr),nchoosek(nforce,2),nmodel)*NaN;
for m = 1:nmodel
    c=0;
    for f1 = 1:nforce-1
        for f2 = f1+1:nforce
            c=c+1;
            rr_force(:,:,c,m)  = [reshape( rr_us(:,:,f2,m), nylat*nxlon,1) reshape( rr_us(:,:,f1,m), nylat*nxlon,1)];
            Nrr_force(:,:,c,m) = hist3(rr_force(:,:,c,m),{xirr yirr});
        end
    end
end
clear c f1 f2

%memory allocation
et_force  = ones(nylat*nxlon,2,nchoosek(nforce,2),nmodel)*NaN;
Net_force = ones(length(xiet),length(yiet),nchoosek(nforce,2),nmodel)*NaN;
for m = 1:nmodel
    c=0;
    for f1 = 1:nforce-1
        for f2 = f1+1:nforce
            c=c+1;
            et_force(:,:,c,m)  = [reshape( et_us(:,:,f2,m), nylat*nxlon,1) reshape( et_us(:,:,f1,m), nylat*nxlon,1)];
            Net_force(:,:,c,m) = hist3(et_force(:,:,c,m),{xiet yiet});
        end
    end
end
clear c f1 f2

et_intf = et_force(:,:,[1;2;3;5;6;8],:);
et_bias = et_force(:,:,[4;7;9;10],:);
Net_intf = Net_force(:,:,[1;2;3;5;6;8],:);
Net_bias = Net_force(:,:,[4;7;9;10],:);

%memory allocation
ro_force  = ones(nylat*nxlon,2,nchoosek(nforce,2),nmodel)*NaN;
Nro_force = ones(length(xiet),length(yiet),nchoosek(nforce,2),nmodel)*NaN;
for m = 1:nmodel
    c=0;
    for f1 = 1:nforce-1
        for f2 = f1+1:nforce
            c=c+1;
            ro_force(:,:,c,m)  = [reshape( ro_us(:,:,f2,m), nylat*nxlon,1) reshape( ro_us(:,:,f1,m), nylat*nxlon,1)];
            Nro_force(:,:,c,m) = hist3(ro_force(:,:,c,m),{xiet yiet});
        end
    end
end
clear c f1 f2

ro_intf = ro_force(:,:,[1;2;3;5;6;8],:);
ro_bias = ro_force(:,:,[4;7;9;10],:);
Nro_intf = Nro_force(:,:,[1;2;3;5;6;8],:);
Nro_bias = Nro_force(:,:,[4;7;9;10],:);
%% -------------------------------------------------------------------------
% Statistics - Spatial bias
% -------------------------------------------------------------------------
ro_force(ro_force>3000)=NaN;
et_force(et_force>3000 | et_force<-100)=NaN;
ro_model(ro_model>3000)=NaN;
et_model(et_model>3000 | et_model<-100)=NaN;

bias_et_force=squeeze( nanmean(et_force(:,1,:,:)-et_force(:,2,:,:)) );
bias_ro_force=squeeze( nanmean(ro_force(:,1,:,:)-ro_force(:,2,:,:)) );

bias_et_model=squeeze( nanmean(et_model(:,1,:,:)-et_model(:,2,:,:)) );
bias_ro_model=squeeze( nanmean(ro_model(:,1,:,:)-ro_model(:,2,:,:)) );

%% -------------------------------------------------------------------------
% Plot setting 
% -------------------------------------------------------------------------
%Import color map for surface map
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
% clr=[0.85 0.6 0.05;0.95 0.8 0.55;0 0.75 0;0.95 0.95 0;0.6 0.6 0.6];
clr=[0.90 0.5 0.05;0.95 0.8 0.55;0 0.75 0;0.15 0.15 0.90;0.5 0.5 0.5];

force_name{1,1}='BCCA';
force_name{2,1}='BCSDd';
force_name{3,1}='BCSDm';
force_name{4,1}='AR';
force_name{5,1}='M02';
%% Scatter-plot between SD and M02 - ET and RO
figure('Units','inches','position',[1 4 6.75 5.50],'Visible','on','Color',[1 1 1]);
% BCCA-M02
for i = 1:6
  subplot(4,6,i)
  if i >=1 && i <=3;  
      imagesc(xiet+25,yiet+25,log10(Net_bias(:,:,1,i)) );
      colormap(mycolor3)
      hold on  
      plot(xiet,yiet,'k:','LineWidth',2);
      plot(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 1; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{''});
      text(1400,-200,force_name{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
      if i==1; title('CLM','FontSize',9); elseif i==2;title('VIC','FontSize',9); elseif i==3; title('PRMS','FontSize',9);end    
  elseif i >=4 && i <=6;
%       surface(Xro,Yro,log10(Nro_bias(:,:,1,i-3)),'EdgeColor','none');
      imagesc(xiro+25,yiro+25,log10(Nro_bias(:,:,1,i-3)) );
      colormap(mycolor3)
      caxis([0 3]);
      hold on
      plot(xiro,yiro,'k:','LineWidth',2);
      plot(ones(length(xiro),1)*700,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*1400,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*2100,linspace(0,2100,length(xiro)),'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*700,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*1400,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*2100,'k:');
      axis equal
      axis xy      
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 4; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{''});
      text(1200,-200,force_name{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90); 
      if i==4; title('CLM','FontSize',9); elseif i==5;title('VIC','FontSize',9); elseif i==6; title('PRMS','FontSize',9);end
  end
  caxis([0 3]);
  set(gca,'FontSize',8);
end
%BCSDd-M02
for i = 1:6
  subplot(4,6,i+6)
  if i >=1 && i <=3;  
%       surface(Xet,Yet,log10(Net_bias(:,:,2,i)),'EdgeColor','none');
      imagesc(xiet+25,yiet+25,log10(Net_bias(:,:,2,i)) );
      colormap(mycolor3);
      caxis([0 3]);
      hold on
      plot(xiet,yiet,'k:','LineWidth',2);
      plot(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 1; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{''});
      text(750,-200,force_name{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,750,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
  elseif i >=4 && i <=6;
%       surface(Xro,Yro,log10(Nro_bias(:,:,2,i-3)),'EdgeColor','none');
      imagesc(xiro+25,yiro+25,log10(Nro_bias(:,:,2,i-3)) );
      caxis([0 3]);
      colormap(mycolor3)
      hold on
      plot(xiro,yiro,'k:','LineWidth',2);
      plot(ones(length(xiro),1)*700,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*1400,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*2100,linspace(0,2100,length(xiro)),'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*700,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*1400,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 4; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{''});
      text(1200,-200,force_name{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);         
  end
  caxis([0 3]);
  set(gca,'FontSize',8);
end
%BCSDm-M02
for i = 1:6
  subplot(4,6,i+12)
  if i >=1 && i <=3;  
%       surface(Xet,Yet,log10(Net_bias(:,:,3,i)),'EdgeColor','none');
      imagesc(xiet+25,yiet+25,log10(Net_bias(:,:,3,i)) );
      colormap(mycolor3);
      caxis([0 3]);
      hold on
      plot(xiet,yiet,'k:','LineWidth',2);
      plot(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 1; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{''});
      text(1200,-200,force_name{3},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
  elseif i >=4 && i <=6;
%       surface(Xro,Yro,log10(Nro_bias(:,:,3,i-3)),'EdgeColor','none');
      imagesc(xiro+25,yiro+25,log10(Nro_bias(:,:,3,i-3)) );
      colormap(mycolor3);
      caxis([0 3]);
      hold on
      plot(xiro,yiro,'k:','LineWidth',2);
      plot(ones(length(xiro),1)*700,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*1400,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*2100,linspace(0,2100,length(xiro)),'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*700,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*1400,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 4; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{''});
      text(1200,-200,force_name{3},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);         
  end
  caxis([0 3]);
  set(gca,'FontSize',8);
end
%AR-M02
for i = 1:6
  subplot(4,6,i+18)
  if i >=1 && i <=3;  
%       surface(Xet,Yet,log10(Net_bias(:,:,4,i)),'EdgeColor','none');
      imagesc(xiet+25,yiet+25,log10(Net_bias(:,:,4,i)) );
      colormap(mycolor3);
      caxis([0 3]);
      hold on
      plot(xiet,yiet,'k:','LineWidth',2);
      plot(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),'k:');
      plot(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,'k:');
      plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 1; set(gca,'YTickLabel',{'0';'';'';'2100'}); 
      else set(gca,'YTickLabel',{''});
      end
      set(gca,'XTickLabel',{'0';'';'';'2100'});
      text(1200,-200,force_name{4},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
  elseif i >=4 && i <=6;
%       surface(Xro,Yro,log10(Nro_bias(:,:,4,i-3)),'EdgeColor','none');
      imagesc(xiro+25,yiro+25,log10(Nro_bias(:,:,4,i-3)) );
      caxis([0 3]);
      colormap(mycolor3);
      hold on
      plot(xiro,yiro,'k:','LineWidth',2);
      plot(ones(length(xiro),1)*700,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*1400,linspace(0,2100,length(xiro)),'k:');
      plot(ones(length(xiro),1)*2100,linspace(0,2100,length(xiro)),'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*700,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*1400,'k:');
      plot(linspace(0,2100,length(xiro)),ones(length(xiro),1)*2100,'k:');
      axis equal
      axis xy
      set(gca,'XLim',[0 2200]);
      set(gca,'YLim',[0 2200]);
      set(gca,'XTick',0:700:2100);
      set(gca,'YTick',0:700:2100);
      if i == 4; set(gca,'YTickLabel',{'0';'';'';'2100'});
      else set(gca,'YTickLabel',{''});
      end          
      set(gca,'XTickLabel',{'0';'';'';'2100'});
      text(1200,-200,force_name{4},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
      text(-200,1200,force_name{5},'HorizontalAlignment','center','FontSize',9,'Rotation',90);         
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
vmove1 = +0.020;
vmove2 = -0.000;
vmove3 = -0.020;
vmove4 = -0.040;
hscale=1.170;
vscale=1.170;
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

axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.750,0.960,'Q [mm/yr]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.325,0.960,'ET [mm/yr]','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
%colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 3]);
colorbar('horizontal','position',[0.775,0.065,0.100,0.0195],...
    'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m4_fig1_wb_comparison_force_M02_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
% set(gcf,'PaperPositionMode','auto')
print('-dtiff','-r300',figfile);

%% Inter-model difference in ET and RO
figure('Units','inches','position',[1 4 6.00 4.25],'Visible','on','color',[1 1 1]);
% ET
for i = 1:3
    subplot(2,3,i)   
    surface(xiet,Yet,log10(Net_model(:,:,5,i)),'EdgeColor','none');
%     imagesc(xiet+25,yiet+25,log10(Net_model(:,:,5,i)) );
    caxis([0 3]);
%     colormap(mycolor3);
    hold on
    plot3(xiet,yiet,5*ones(length(xiet),1),'k:','LineWidth',2);
    plot3(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),5*ones(length(xiet),1),'k:');
    plot3(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),5*ones(length(xiet),1),'k:');
    plot3(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),5*ones(length(xiet),1),'k:');
    plot3(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,5*ones(length(xiet),1),'k:');
    plot3(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,5*ones(length(xiet),1),'k:');
    plot3(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,5*ones(length(xiet),1),'k:');
    axis equal
%     axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 2200]);
    set(gca,'YLim',[0 2200]);
    set(gca,'XTick',0:700:2100);
    set(gca,'YTick',0:700:2100);
    box on
    if i == 1
        text(1200,-200,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-200,1200,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(1200,-200,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-200,1200,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(1200,-200,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-200,1200,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end  
    set(gca,'YTickLabel',{'0';'';'';'2100'});
    set(gca,'XTickLabel',{'0';'';'';'2100'});
end
% RO
for i = 1:3
    subplot(2,3,i+3)
    surface(Xro,Yro,log10(Nro_model(:,:,5,i)),'EdgeColor','none');
%     imagesc(xiro+25,yiro+25,log10(Nro_model(:,:,5,i)) );
    caxis([0 3]);
%     colormap(mycolor3);
    hold on
    plot3(xiro,yiro,5*ones(length(xiro),1),'k:','LineWidth',2);
    plot3(ones(length(xiro),1)*700,linspace(0,2100,length(xiro)),5*ones(length(xiro),1),'k:');
    plot3(ones(length(xiro),1)*1400,linspace(0,2100,length(xiro)),5*ones(length(xiro),1),'k:');
    plot3(ones(length(xiro),1)*2100,linspace(0,2100,length(xiro)),5*ones(length(xiro),1),'k:');
    plot3(linspace(0,2100,length(xiro)),ones(length(xiro),1)*700,5*ones(length(xiro),1),'k:');
    plot3(linspace(0,2100,length(xiro)),ones(length(xiro),1)*1400,5*ones(length(xiro),1),'k:');
    plot3(linspace(0,2100,length(xiro)),ones(length(xiro),1)*2100,5*ones(length(xiro),1),'k:');
    axis equal
%     axis xy
    set(gca,'FontSize',8);
    set(gca,'XLim',[0 2200]);
    set(gca,'YLim',[0 2200]);
    set(gca,'XTick',0:700:2100);
    set(gca,'YTick',0:700:2100);
    box on
    if i == 1
        text(1200,-200,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-200,1200,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==2
        text(1200,-200,model{1},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-200,1200,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    elseif i ==3
        text(1200,-200,model{2},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
        text(-200,1200,model{3},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
    end
    set(gca,'YTickLabel',{'0';'';'';'2100'});
    set(gca,'XTickLabel',{'0';'';'';'2100'});

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
text(0.525,0.485,'Q [mm/yr]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.525,0.950,'ET [mm/yr]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);

%colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 3]);
colorbar('horizontal','position',[0.775,0.050,0.100,0.0195],...
    'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m4_fig2_wb_comparison_model_M02_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
% set(gcf,'PaperPositionMode','auto')
print('-dtiff','-r300',figfile);
%% Scatter plot of inter forcing comparison
% for m=1:nmodel
%     figure('Units','inches','position',[1 4 6.75 6.75],'Visible','on','Color',[1 1 1]);
%     % ET - upper right half
%     c=0;
%     for row=1:nsd-1
%         c1=0;
%         p1=nsd*(row-1)+row+1;
%         for pan=p1:row*nsd
%             c1=c1+1;
%             c=c+1;
%             subplot(4,4,pan)
% %             surface(Xet,Yet,log10(Net_intf(:,:,c,m)),'EdgeColor','none');
%             imagesc(xiet+25,yiet+25,log10(Net_intf(:,:,c,m)) );            
%             caxis([0 3]);
%             colormap(mycolor3);
%             hold on
%             plot3(xiet,yiet,5*ones(length(xiet),1),'k:','LineWidth',2);
%             plot3(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),5*ones(length(xiet),1),'k:');
%             plot3(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),5*ones(length(xiet),1),'k:');
%             plot3(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),5*ones(length(xiet),1),'k:');
%             plot3(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,5*ones(length(xiet),1),'k:');
%             plot3(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,5*ones(length(xiet),1),'k:');
%             plot3(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,5*ones(length(xiet),1),'k:');
%             axis equal
%             axis xy
%             set(gca,'FontSize',8);
%             set(gca,'XLim',[0 2200]);
%             set(gca,'YLim',[0 2200]);
%             set(gca,'XTick',0:700:2100);
%             set(gca,'YTick',0:700:2100);
%             text(1400,200,force_name{row},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
%             text(-200,1200,force_name{c1+row},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
%             if pan~=p1;
%                 set(gca,'YTickLabel',{''});
%                 set(gca,'XTickLabel',{''});
%             else
%                 set(gca,'YTickLabel',{'0';'';'';'2100'});
%                 set(gca,'XTickLabel',{'0';'';'';'2100'});
%             end
%         end
%     end
% 
%     % RO - lower left half
%     rowoffset=[5;10;15;20];
%     c=0;
%     for col=1:nsd-1
%         p1 = nsd*col+col;
%         c1=0;
%         for pan=p1:nsd:(nsd-1)*nsd+col
%             c=c+1;
%             c1=c1+1;
%             subplot(4,4,pan)
% %             surface(Xro,Yro,log10(Nro_intf(:,:,c,m)),'EdgeColor','none');
%             imagesc(xiro+25,yiro+25,log10(Nro_intf(:,:,c,m)) );
%             caxis([0 3]);
%             colormap(mycolor3);
%             hold on
%             plot(xiet,yiet,'k:','LineWidth',2);
%             plot(ones(length(xiet),1)*700,linspace(0,2100,length(xiet)),'k:');
%             plot(ones(length(xiet),1)*1400,linspace(0,2100,length(xiet)),'k:');
%             plot(ones(length(xiet),1)*2100,linspace(0,2100,length(xiet)),'k:');
%             plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*700,'k:');
%             plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*1400,'k:');
%             plot(linspace(0,2100,length(xiet)),ones(length(xiet),1)*2100,'k:');
%             axis equal
%             axis xy
%             set(gca,'FontSize',8);
%             set(gca,'XLim',[0 2200]);
%             set(gca,'YLim',[0 2200]);
%             set(gca,'XTick',0:700:2100);
%             set(gca,'YTick',0:700:2100);
%             text(1200,-200,force_name{col},'HorizontalAlignment','center','FontSize',9,'Rotation',0);
%             text(-200,1200,force_name{c1+col},'HorizontalAlignment','center','FontSize',9,'Rotation',90);
%             if row~=1
%                 set(gca,'YTickLabel',{''});
%             else
%                 set(gca,'YTickLabel',{'0';'';'';'2100'});
%             end
%             if pan<20
%                 set(gca,'XTickLabel',{''});
%             else
%                 set(gca,'XTickLabel',{'0';'';'';'2100'});
%             end
%         end
%     end
%    
%     %adjust plot size & location
%     h   = get(gcf,'Children');
%     for i=1:length(h)
%         pos(i,:) =get(h(i),'position');
%     end
%     %from left to right
%     hmove1=+0.015;
%     hmove2=+0.000;
%     hmove3=-0.015;
%     hmove4=-0.030;
%     hmove5=-0.045;
%     %from bottom to top
%     vmove1=-0.015;
%     vmove2=-0.010;
%     vmove3=-0.005;
%     vmove4=+0.000;
%     vmove5=+0.005;
%     
%     hscale=1.1275;
%     vscale=1.1275;
%     
%     %bottom row fromright to left
%     set(h(1), 'Position',[ pos(1,1)+hmove2  pos(1,2)+vmove1  pos(1,3)*hscale pos(1,4)*vscale])
%     
%     set(h(2), 'Position',[ pos(2,1)+hmove3  pos(2,2)+vmove1  pos(2,3)*hscale pos(2,4)*vscale])
%     set(h(3), 'Position',[ pos(3,1)+hmove3  pos(3,2)+vmove2  pos(3,3)*hscale pos(3,4)*vscale])
%     
%     set(h(4), 'Position',[ pos(4,1)+hmove4  pos(4,2)+vmove1  pos(4,3)*hscale pos(4,4)*vscale])
%     set(h(5), 'Position', [ pos(5,1)+hmove4  pos(5,2)+vmove2  pos(5,3)*hscale pos(5,4)*vscale])
%     set(h(6), 'Position', [ pos(6,1)+hmove4  pos(6,2)+vmove3  pos(6,3)*hscale pos(6,4)*vscale])
%     
%     set(h(7), 'Position', [ pos(7,1)+hmove1  pos(7,2)+vmove2  pos(7,3)*hscale pos(7,4)*vscale])
%     
%     set(h(8), 'Position', [ pos(8,1)+hmove1  pos(8,2)+vmove3  pos(8,3)*hscale pos(8,4)*vscale])
%     set(h(9),  'Position',[ pos(9,1)+hmove2  pos(9,2)+vmove3  pos(9,3)*hscale  pos(9,4)*vscale])
%     
%     set(h(10), 'Position',[ pos(10,1)+hmove1 pos(10,2)+vmove4 pos(10,3)*hscale pos(10,4)*vscale])
%     set(h(11), 'Position',[ pos(11,1)+hmove2 pos(11,2)+vmove4 pos(11,3)*hscale pos(11,4)*vscale]) 
%     set(h(12), 'Position',[ pos(12,1)+hmove3 pos(12,2)+vmove4 pos(12,3)*hscale pos(12,4)*vscale])
%     
%     % Text
%     axes('Units','Normalized','position',[0 0 1 1],'visible','off');
%     text(0.500,0.045,'RO [mm]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
%     text(0.960,0.550,'ET [mm]', 'HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',90);
%     
%     text(0.925,0.025,'log_{10}(N)','HorizontalAlignment','center','FontSize',8,'FontName','Helvetica','Rotation',0);
%     %colorbar, text on the plots
%     axes('Units','Normalized','position',[0 0 1 1],'visible','off');
%     caxis([0 3]);
%     colorbar('horizontal','position',[0.775,0.040,0.100,0.0195],...
%         'XLim',[0 3],'XTick',0:1.5:3,'FontSize',8);
%     cbfreeze;
%     
%     %save figure
%     figfile=['./figure/paper/Paper_m4_fig1_wb_comparison_force_' model{m} '_' num2str(wyr(1)) '-' num2str(wyr(end)) '_us_' num2str(res) 'K'];
% %     set(gcf,'PaperPositionMode','auto')
%     print('-dtiff','-r300',figfile);
% end