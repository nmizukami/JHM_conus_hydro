% Comparison of monthly anomaly
%          plot_type     temp_res.  time span   climatology  spatial res.   variables
%  ---------------------------------------------------------------------------------
%  Fig 1-1
%  Fig 1-2 
%  Fig 1-3 
%  Fig 1-4 

% profile on

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
% byr=2001; bmon=10; bday=1;
%eyr=2008; emon=9; eday=30;
eyr=1999; emon=9; eday=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -Variables that cannot be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% water year array 
% -------------------------------------------------------------------------
wyr = byr+1:1:eyr;
% number of dimension
nwyr    = length(wyr);
nregion = length(region);
nmodel  = length(model);
nprdct  = length(prdct);
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
huc2_shp = '/d3/mizukami/domain_huc/shapefile/HUC/HUC02_conus.shp';
S = shaperead(huc2_shp);
clear huc2_shp

clr3=[0 0 1;1 0 0;0 1 0];
clr5=[0 0 1;0.65 0.16 0.16;1 0 0;0 1 0;0 0 0];
%%
% Memory allocation
% [lat x lon x days x wyrs x prods] 
%pr = zeros(ylat_num,xlon_num,12*nwyr,nmodel,nprdct);
ro = zeros(ylat_num,xlon_num,12*nwyr,nmodel,nprdct);
et = zeros(ylat_num,xlon_num,12*nwyr,nmodel,nprdct);

%prc = zeros(ylat_num,xlon_num,12,nmodel,nprdct);
roc = zeros(ylat_num,xlon_num,12,nmodel,nprdct);
etc = zeros(ylat_num,xlon_num,12,nmodel,nprdct);

%Create simulation area mask
sim_mask = lon2d+huc2g;
sim_mask(~isnan(sim_mask) ) = 1;
sim_mask_yr = repmat(sim_mask,[1,1,12*nwyr,nmodel,nprdct]);

%Create simulation area mask
sim_mask_c = repmat(sim_mask,[1,1,12,nmodel,nprdct]);

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
  for p=1:length(prdct)
      for m=1:length(model)
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
          % Reading monthly grid PRCP
%          ncname = [main_indir '/' prdct{p} '_' region{r} '/calib/Monthly_PRCP.nc'];
%          [var0,fillvalue]=netcdf2mat(ncname,'PRCP','attname1','_FillValue');
%          var0(var0==fillvalue)=0;
%          pr(i1:i2,j1:j2,:,m,p) = var0+squeeze(pr(i1:i2,j1:j2,:,m,p));
%          clear var0 ncname
%          
%          % Reading monthly climatelogy grid PRCP
%          ncname = [main_indir '/' prdct{p} '_' region{r} '/calib/Monthly_clim_PRCP.nc'];
%          [var0,fillvalue]=netcdf2mat(ncname,'PRCP','attname1','_FillValue');
%          var0(var0==fillvalue)=0;
%          prc(i1:i2,j1:j2,:,m,p) = var0+squeeze(prc(i1:i2,j1:j2,:,m,p));
%          clear var0 ncname
          
          % Reading monthly grid ET
          ncname = [main_indir '/' prdct{p} '_' region{r} '/calib/Monthly_ET.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'ET','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          et(i1:i2,j1:j2,:,m,p) = var0+squeeze(et(i1:i2,j1:j2,:,m,p));
          clear var0 ncname
          
          % Reading monthly climatelogy grid ET
          ncname = [main_indir '/' prdct{p} '_' region{r} '/calib/Monthly_clim_ET.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'ET','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          etc(i1:i2,j1:j2,:,m,p) = var0+squeeze(etc(i1:i2,j1:j2,:,m,p));
          clear var0 ncname
          
          % Reading monthly grid RUNOFF
          ncname = [main_indir '/' prdct{p} '_' region{r} '/calib/Monthly_RUNOFF.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'RUNOFF','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          ro(i1:i2,j1:j2,:,m,p) = var0+ro(i1:i2,j1:j2,:,m,p);
          clear var0 ncname
          
          % Reading monthly climatelogy grid runoff
          ncname = [main_indir '/' prdct{p} '_' region{r} '/calib/Monthly_clim_RUNOFF.nc'];
          [var0,fillvalue]=netcdf2mat(ncname,'RUNOFF','attname1','_FillValue');
          var0(var0==fillvalue)=0;
          roc(i1:i2,j1:j2,:,m,p) = var0+squeeze(roc(i1:i2,j1:j2,:,m,p));
          clear var0 ncname          
      end
  end
end
%Mask out the outside the simulation domain
%pr = pr.*sim_mask_yr;
et = et.*sim_mask_yr;
ro = ro.*sim_mask_yr;

%prc = prc.*sim_mask_c;
etc = etc.*sim_mask_c;
roc = roc.*sim_mask_c;

clear sim_mask sim_mask_yr sim_mask_c
%% compute monhtly anomaly series grid
tic
mon = repmat(1:1:12,1,nwyr);
%pra=pr;
eta=et;
roa=ro;
for i = 1:12*nwyr
    for m = 1:nmodel
        for p = 1:nprdct
%            pra(:,:,i,m,p)=pr(:,:,i,m,p) - prc(:,:,mon(i),m,p);
            eta(:,:,i,m,p)=et(:,:,i,m,p) - etc(:,:,mon(i),m,p);
            roa(:,:,i,m,p)=ro(:,:,i,m,p) - roc(:,:,mon(i),m,p);
        end
    end
end
toc
%% inter model difference for each forcing
% nReps=100;
% R_eta_model_ens = ones(ylat_num,xlon_num,nReps,length(prdct))*NaN;
% R_roa_model_ens = ones(ylat_num,xlon_num,nReps,length(prdct))*NaN;
% for j=1:nReps
%     tic
%     eta_samp_mean = ones(ylat_num,xlon_num,12*length(wyr),length(prdct))*NaN;
%     roa_samp_mean = ones(ylat_num,xlon_num,12*length(wyr),length(prdct))*NaN;
%     eta_samp_std  = ones(ylat_num,xlon_num,12*length(wyr),length(prdct))*NaN;
%     roa_samp_std  = ones(ylat_num,xlon_num,12*length(wyr),length(prdct))*NaN;
%     for i = 1:12*length(wyr)
%         sampI = ceil(rand(length(model),1)*3);
%         eta_samp = eta(:,:,i,sampI,:);
%         roa_samp = roa(:,:,i,sampI,:);
%         
%         eta_samp_mean(:,:,i,:)=squeeze( mean(eta_samp,4) );
%         roa_samp_mean(:,:,i,:)=squeeze( mean(roa_samp,4) );
%         
%         eta_samp_std(:,:,i,:)=squeeze( std(eta_samp,1,4) );
%         roa_samp_std(:,:,i,:)=squeeze( std(roa_samp,1,4) );
%     end
%     R_eta_model_ens(:,:,j,:)=squeeze( mean(eta_samp_std,3)./std(eta_samp_mean,0,3) );
%     R_roa_model_ens(:,:,j,:)=squeeze( mean(roa_samp_std,3)./std(roa_samp_mean,0,3) );
%     toc
% end
% 
% %% inter forcing difference for each forcing
% R_eta_force_ens = ones(ylat_num,xlon_num,nReps,nmodel)*NaN;
% R_roa_force_ens = ones(ylat_num,xlon_num,nReps,nmodel)*NaN;
% for j=1:nReps
%     tic
%     eta_samp_mean=ones(ylat_num,xlon_num,12*nwyr,nmodel)*NaN;
%     roa_samp_mean=ones(ylat_num,xlon_num,12*nwyr,nmodel)*NaN;
%     eta_samp_std=ones(ylat_num,xlon_num,12*nwyr,nmodel)*NaN;
%     roa_samp_std=ones(ylat_num,xlon_num,12*nwyr,nmodel)*NaN;
%     for i = 1:12*length(wyr)
%         sampI = ceil(rand(length(prdct)-1,1)*4);
%         eta_samp = eta(:,:,i,:,sampI);
%         roa_samp = roa(:,:,i,:,sampI);
%         
%         eta_samp_mean(:,:,i,:) = squeeze( mean(eta_samp,5));
%         roa_samp_mean(:,:,i,:) = squeeze( mean(roa_samp,5));
%         
%         eta_samp_std (:,:,i,:) = squeeze( std(eta_samp,1,5));
%         roa_samp_std (:,:,i,:) = squeeze( std(roa_samp,1,5));
%     end
%     R_eta_force_ens(:,:,j,:)=squeeze( mean(eta_samp_std,3) ./ std(eta_samp_mean,0,3) );
%     R_roa_force_ens(:,:,j,:)=squeeze( mean(roa_samp_std,3) ./ std(roa_samp_mean,0,3) );
%     toc
% end

% %% inter model difference for each forcing
% nReps=100;
% tic
% R_eta_model=ones(ylat_num,xlon_num,nReps,length(prdct))*NaN;
% R_roa_model=ones(ylat_num,xlon_num,nReps,length(prdct))*NaN;
% for j=1:nReps
%     sampI = ceil(rand(length(model),1)*3);
%     eta_samp = eta(:,:,:,sampI,:);
%     roa_samp = roa(:,:,:,sampI,:);
%     
%     eta_samp_mean=squeeze( mean(eta_samp,4) );
%     roa_samp_mean=squeeze( mean(roa_samp,4) );
%     
%     eta_samp_std=squeeze( std(eta_samp,1,4) );
%     roa_samp_std=squeeze( std(roa_samp,1,4) );
% 
%     R_eta_model(:,:,j,:)=squeeze( mean(eta_samp_std,3)./std(eta_samp_mean,0,3) );
%     R_roa_model(:,:,j,:)=squeeze( mean(roa_samp_std,3)./std(roa_samp_mean,0,3) );
% end
% toc
% %% inter forcing difference for each forcing
% tic

% R_eta_force=ones(ylat_num,xlon_num,nReps,length(model))*NaN;
% R_roa_force=ones(ylat_num,xlon_num,nReps,length(model))*NaN;
% for j=1:nReps
%     sampI = ceil(rand(length(prdct)-1,1)*4);
%     eta_samp = eta(:,:,:,:,sampI);
%     roa_samp = roa(:,:,:,:,sampI);
%     
%     eta_samp_mean=squeeze( mean(eta_samp,5) );
%     roa_samp_mean=squeeze( mean(roa_samp,5) );
%     
%     eta_samp_std=squeeze( std(eta_samp,1,5) );
%     roa_samp_std=squeeze( std(roa_samp,1,5) );
%     
%     R_eta_force(:,:,j,:)=squeeze( mean(eta_samp_std,3) ./std(eta_samp_mean,0,3) );
%     R_roa_force(:,:,j,:)=squeeze( mean(roa_samp_std,3) ./std(roa_samp_mean,0,3) );
% end
% toc
% %%
% R_roa_model_mean=squeeze( mean(R_roa_model_ens,3));
% R_eta_model_mean=squeeze( mean(R_eta_model_ens,3));
% 
% R_roa_force_mean=squeeze( mean(R_roa_force_ens,3));
% R_eta_force_mean=squeeze( mean(R_eta_force_ens,3));

%% inter-model spread for each forcing and inter-forcing spread for each model
tic
R_model = @(x) mean(std(x,1,4),3) ./ std(mean(x,4),0,3);
R_eta_model = squeeze( R_model(eta) );
R_roa_model = squeeze( R_model(roa) );
toc
tic
R_force = @(x) mean(std(x,1,5),3) ./ std(mean(x,5),0,3);
R_eta_force = squeeze( R_force(eta(:,:,:,:,1:4)) );
R_roa_force = squeeze( R_force(roa(:,:,:,:,1:4)) );

R_eta_force_no_bcca = squeeze( R_force(eta(:,:,:,:,2:4)) );
R_roa_force_no_bcca = squeeze( R_force(roa(:,:,:,:,2:4)) );

%R_pra_force = squeeze( R_force(pra(:,:,:,:,2:4)) );
toc
% %% bootstrapped inter model difference for each forcing
% nReps=100;
% tic
% R_eta_model_ens = ones(ylat_num,xlon_num,nReps,nprdct)*NaN;
% R_roa_model_ens = ones(ylat_num,xlon_num,nReps,nprdct)*NaN;
% for j=1:nReps
%     sampI = ceil(rand(12*nwyr,1)*12*nwyr);
%     R_eta_model_ens(:,:,j,:)=squeeze( R_model(eta(:,:,sampI,:,:)) );
%     R_roa_model_ens(:,:,j,:)=squeeze( R_model(roa(:,:,sampI,:,:)) );
% end
% toc
% %% bootstrapped inter forcing difference for each forcing
% tic
% R_eta_force_ens=ones(ylat_num,xlon_num,nReps,nmodel)*NaN;
% R_roa_force_ens=ones(ylat_num,xlon_num,nReps,nmodel)*NaN;
% for j=1:nReps
%     sampI = ceil(rand(12*nwyr,1)*12*nwyr);
%     R_eta_force_ens(:,:,j,:) = squeeze( R_force(eta(:,:,sampI,:,1:4)) );
%     R_roa_force_ens(:,:,j,:) = squeeze( R_force(roa(:,:,sampI,:,1:4)) );
% end
% toc
%% sort array
R_roa_model_1d=reshape(R_roa_model,ylat_num*xlon_num,5);
R_roa_force_1d=reshape(R_roa_force,ylat_num*xlon_num,3);
R_roa_force_no_bcca_1d=reshape(R_roa_force_no_bcca,ylat_num*xlon_num,3);

R_roa_model_sort=sort(R_roa_model_1d,1,'ascend');
R_roa_force_sort=sort(R_roa_force_1d,1,'ascend');
R_roa_force_no_bcca_sort=sort(R_roa_force_no_bcca_1d,1,'ascend');

R_roa_model_sort(any(isnan(R_roa_model_sort),2),:)=[];
R_roa_force_sort(any(isnan(R_roa_force_sort),2),:)=[];
R_roa_force_no_bcca_sort(any(isnan(R_roa_force_no_bcca_sort),2),:)=[];

%% CDF plot
mlen=size(R_roa_model_sort,1);
xm=1:mlen;
flen=size(R_roa_force_sort,1);
xf=1:flen;

figure('Name','1-1 inter-model vs inter-SD - RO','NumberTitle','off','Units','inches','position',[1 7 5.5 5.5],'Visible','on','Color','w');
subplot(2,1,1)
set(gca,'ColorOrder',clr5);
for i=1:size(R_roa_model_sort,2)
    plot(R_roa_model_sort(:,i),xm/mlen,'LineWidth',2,'Color',clr5(i,:));hold on
end
set(gca,'XLim',[0 1.2]);set(gca,'YLim',[0 1.00]);
xlabel('\itR^{model}_{std}')
ylabel('CDF')
title('a) Inter-model spread')
set(gca,'FontSize',9)

subplot(2,1,2)
for i=1:size(R_roa_force_sort,2)
    plot(R_roa_force_sort(:,i),xf/flen,'LineWidth',2,'Color',clr3(i,:));hold on
end
set(gca,'XLim',[0 1.2]);set(gca,'YLim',[0 1.00]);
xlabel('\itR^{forc}_{std}')
ylabel('CDF')
title('b) Inter-SD spread')
set(gca,'FontSize',9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove = +0.0;
vmove1 = -0.005;
vmove2 = +0.005;
hscale=1.01;
vscale=1.05;

%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove pos(1,2)+vmove1 pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove pos(2,2)+vmove2 pos(2,3)*hscale pos(2,4)*vscale])

%legends 1
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr5(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr5(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr5(3,:));
dummy4 = plot(1,1,'-','LineWidth',2,'color',clr5(4,:));
dummy5 = plot(1,1,'-','LineWidth',2,'color',clr5(5,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off'); set(dummy5,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('BCCA','BCSDd','BCSDm','AR','M02',...
    'Orientation','vertical','Location',[0.815 0.825 0.0075 0.030]);
legend('boxoff','Color','w');

%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr3(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr3(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr3(3,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); 
set(gca,'Visible','off','FontSize',9)
legend('CLM','VIC','PRMS',...
    'Orientation','vertical','Location',[0.815 0.375 0.0075 0.030]);
legend('boxoff','Color','none');
%% histogram
figure('Name','1-1 inter-model vs inter-SD - RO','NumberTitle','off','Units','inches','position',[1 7 5.5 6.0],'Visible','on','Color','w');
subplot(2,1,2)
set(gca,'ColorOrder',clr5);
% for i=1:100
%     [n,xout] = hist(reshape(R_roa_model_ens(:,:,i,:),ylat_num*xlon_num,5),100);
%     plot(xout,n,'LineWidth',0.5);hold on
% end
[n,xout] = hist(R_roa_model_1d,80);
for i=1:size(n,2)
    plot(xout,n(:,i),'LineWidth',2,'Color',clr5(i,:));hold on
end
set(gca,'XLim',[0 1.00]);set(gca,'YLim',[0 5500]);
ylabel('count')
title('b) Inter-model spread index')
set(gca,'FontSize',9)

subplot(2,1,1)
% for i=1:100
%     [n,xout] = hist(reshape(R_roa_force_ens(:,:,i,:),ylat_num*xlon_num,3),50);
%     plot(xout,n,'LineWidth',0.5);hold on
% end
[n,xout] = hist(R_roa_force_1d,80);
[n_no_bcca,xout_no_bcca] = hist(R_roa_force_no_bcca_1d,80);
for i=1:size(n,2)
    plot(xout_no_bcca,n_no_bcca(:,i),'LineWidth',2,'Color',clr3(i,:));hold on
    plot(xout,n(:,i),'--','LineWidth',1,'Color',clr3(i,:));
end
set(gca,'XLim',[0 1.00]);set(gca,'YLim',[0 5500]);
ylabel('count')
title('a) Inter-forcing spread index')
set(gca,'FontSize',9)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove = +0.0;
vmove1 = -0.005;
vmove2 = +0.005;
hscale=1.01;
vscale=1.05;

%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove pos(1,2)+vmove1 pos(1,3)*hscale pos(1,4)*vscale])
set(h(2), 'Position',[ pos(2,1)+hmove pos(2,2)+vmove2 pos(2,3)*hscale pos(2,4)*vscale])

%Text
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
text(0.50,0.525,'\itR^{ forcing}_{std}','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);
text(0.50,0.045,'\itR^{ model}_{std}','HorizontalAlignment','center','FontSize',9,'FontName','Helvetica','Rotation',0);

%legends 1
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr5(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr5(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr5(3,:));
dummy4 = plot(1,1,'-','LineWidth',2,'color',clr5(4,:));
dummy5 = plot(1,1,'-','LineWidth',2,'color',clr5(5,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off'); set(dummy5,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('BCCA','BCSDd','BCSDm','AR','M02',...
    'Orientation','vertical','Location',[0.800 0.375 0.0075 0.030]);
legend('boxoff','Color','w');

%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr3(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr3(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr3(3,:));
dummy4 = plot(1,1,'--','LineWidth',1,'color','k');
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off');
set(gca,'Visible','off','FontSize',9)
legend('CLM','VIC','PRMS','w/ BCCA',...
    'Orientation','vertical','Location',[0.800 0.825 0.0075 0.030]);
legend('boxoff','Color','none');

figfile=['./figure/paper/Paper_m10_fig1_histo_spread_RO_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-r300',figfile);
%%
figure('Name','1-2 inter-model-RO','NumberTitle','off','Units','inches','position',[1 7 6 6],'Visible','on','Color','w');
set(0,'DefaultAxesColorOrder',clr5);

subplot(3,2,1)
ax1 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax1, 'Visible','off','layer','top');
setm(ax1,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_model(:,:,1));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('a) BCCA')

subplot(3,2,2)
ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax2, 'Visible','off','layer','top');
setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_model(:,:,2));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('b) BCSDd')

subplot(3,2,3)
ax3 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax3, 'Visible','off','layer','top');
setm(ax3,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_model(:,:,3));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('c) BCSDm')

subplot(3,2,4)
ax4 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax4, 'Visible','off','layer','top');
setm(ax4,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_model(:,:,4));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('d) AR')

subplot(3,2,5)
ax5 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax5, 'Visible','off','layer','top');
setm(ax5,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_model(:,:,5));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('f) M02')

subplot(3,2,6)
% for i=1:100
%     [n,xout] = hist(reshape(R_roa_model_ens(:,:,i,:),ylat_num*xlon_num,5),100);
%     plot(xout,n,'LineWidth',0.5);hold on
% end
[n,xout] = hist(reshape(R_roa_model,ylat_num*xlon_num,5),100);
plot(xout,n,'LineWidth',2);hold on
set(gca,'XLim',[0 1.00]);
xlabel('\itR^{model}_{std}')
ylabel('count')
title('e) histogram')
set(gca,'FontSize',9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.015;
hmove2 = -0.025;
vmove  = -0.005;
hscale1=1.115;
vscale1=1.115;
hscale2=1.135;
vscale2=1.135;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1 pos(1,2)+vmove pos(1,3)*hscale1 pos(1,4)*vscale1])
set(h(2), 'Position',[ pos(2,1)+hmove2 pos(2,2)+vmove pos(2,3)*hscale2 pos(2,4)*vscale2])

set(h(3), 'Position',[ pos(3,1)+hmove1 pos(3,2)+vmove pos(3,3)*hscale2 pos(3,4)*vscale2])
set(h(4), 'Position',[ pos(4,1)+hmove2 pos(4,2)+vmove pos(4,3)*hscale2 pos(4,4)*vscale2])

set(h(5), 'Position',[ pos(5,1)+hmove1 pos(5,2)+vmove pos(5,3)*hscale2 pos(5,4)*vscale2])
set(h(6), 'Position',[ pos(6,1)+hmove2 pos(6,2)+vmove pos(6,3)*hscale2 pos(6,4)*vscale2])

%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.45,0.05,'\itR^{model}_{std}','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);

%legends 1
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr5(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr5(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr5(3,:));
dummy4 = plot(1,1,'-','LineWidth',2,'color',clr5(4,:));
dummy5 = plot(1,1,'-','LineWidth',2,'color',clr5(5,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off'); set(dummy5,'Visible','off');
set(gca,'Visible','off','FontSize',7)
legend('BCCA','BCSDd','BCSDm','AR','M02',...
    'Orientation','vertical','Location',[0.675 0.325 0.0075 0.030]);
legend('boxon','Color','w');

%CLM  difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 1.2]); 
colorbar('horizontal','position',[0.15,0.050,0.25,0.035],...
   'XLim',[0 1.2],'XTick',0:0.3:1.2,'FontSize',10);
cbfreeze;

figfile=['./figure/paper/Paper_m10_fig2_intermodel_RO_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);
%%
figure('Name','1-1 inter-model-ET','NumberTitle','off','Units','inches','position',[1 7 6 6],'Visible','on','Color','w');
set(0,'DefaultAxesColorOrder',clr5);

subplot(3,2,1)
ax1 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax1, 'Visible','off','layer','top');
setm(ax1,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_model(:,:,1));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('a) BCCA')

subplot(3,2,2)
ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax2, 'Visible','off','layer','top');
setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_model(:,:,2));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('b) BCSDd')

subplot(3,2,3)
ax3 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax3, 'Visible','off','layer','top');
setm(ax3,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_model(:,:,3));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('c) BCSDm')

subplot(3,2,4)
ax4 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax4, 'Visible','off','layer','top');
setm(ax4,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_model(:,:,4));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('d) AR')

subplot(3,2,5)
ax5 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax5, 'Visible','off','layer','top');
setm(ax5,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_model(:,:,5));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('e) M02')

subplot(3,2,6)
% for i=1:80
%     [n,xout] = hist(reshape(R_eta_model_ens(:,:,i,:),ylat_num*xlon_num,5),50);
%     plot(xout,n,'LineWidth',0.5);hold on
% end
[n,xout] = hist(reshape(R_eta_model,ylat_num*xlon_num,5),50);
plot(xout,n,'LineWidth',2);hold on
set(gca,'XLim',[0 1.25]);
xlabel('\itR^{model}_{std}')
ylabel('count')
title('f) histogram')
set(gca,'FontSize',9)

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.015;
hmove2 = -0.025;
vmove  = -0.005;
hscale1=1.115;
vscale1=1.115;
hscale2=1.135;
vscale2=1.135;
%bottom row from right to left
set(h(1), 'Position',[ pos(1,1)+hmove1 pos(1,2)+vmove pos(1,3)*hscale1 pos(1,4)*vscale1])
set(h(2), 'Position',[ pos(2,1)+hmove2 pos(2,2)+vmove pos(2,3)*hscale2 pos(2,4)*vscale2])

set(h(3), 'Position',[ pos(3,1)+hmove1 pos(3,2)+vmove pos(3,3)*hscale2 pos(3,4)*vscale2])
set(h(4), 'Position',[ pos(4,1)+hmove2 pos(4,2)+vmove pos(4,3)*hscale2 pos(4,4)*vscale2])

set(h(5), 'Position',[ pos(5,1)+hmove1 pos(5,2)+vmove pos(5,3)*hscale2 pos(5,4)*vscale2])
set(h(6), 'Position',[ pos(6,1)+hmove2 pos(6,2)+vmove pos(6,3)*hscale2 pos(6,4)*vscale2])

%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.45,0.05,'\itR^{model}_{std}','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);

%legends 1
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr5(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr5(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr5(3,:));
dummy4 = plot(1,1,'-','LineWidth',2,'color',clr5(4,:));
dummy5 = plot(1,1,'-','LineWidth',2,'color',clr5(5,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); set(dummy4,'Visible','off'); set(dummy5,'Visible','off');
set(gca,'Visible','off','FontSize',7)
legend('BCCA','BCSDd','BCSDm','BCSAR','M02',...
    'Orientation','vertical','Location',[0.895 0.32 0.01 0.025]);
legend('boxon','Color','w');

%CLM  difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 1.2]); 
colorbar('horizontal','position',[0.15,0.050,0.25,0.035],...
   'XLim',[0 1.2],'XTick',0:0.3:1.2,'FontSize',10);
cbfreeze;

figfile=['./figure/paper/Paper_m10_fig3_intermodel_ET_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);
%%
figure('Name','1-3 inter-forc-RO','NumberTitle','off','Units','inches','position',[1 7 6 5],'Color','w');
set(0,'DefaultAxesColorOrder',clr3);

subplot(2,2,1)
ax1 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax1, 'Visible','off','layer','top');
setm(ax1,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_force(:,:,1));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('a) CLM')

subplot(2,2,2)
ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax2, 'Visible','off','layer','top');
setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_force(:,:,2));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('b) VIC')

subplot(2,2,3)
ax3 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax3, 'Visible','off','layer','top');
setm(ax3,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_roa_force(:,:,3));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('c) PRMS')

subplot(2,2,4)
% for i=1:100
%     [n,xout] = hist(reshape(R_roa_force_ens(:,:,i,:),ylat_num*xlon_num,3),50);
%     plot(xout,n,'LineWidth',0.5);hold on
% end
[n,xout] = hist(reshape(R_roa_force,ylat_num*xlon_num,3),50);
plot(xout,n,'LineWidth',2);hold on
set(gca,'XLim',[0 1.00]);
xlabel('\itR^{forc}_{std}')
ylabel('count')
title('d) histogram')
set(gca,'FontSize',9)

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.010;
hmove2 = -0.075;
vmove1 = +0.005;
vmove2 = -0.005;
hscale1=1.05;
vscale1=1.05;
hscale2=1.20;
vscale2=1.20;
%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hmove1 pos(1,2)+vmove1 pos(1,3)*hscale1 pos(1,4)*vscale1])
set(h(2), 'Position',[ pos(2,1)+hmove2 pos(2,2)+vmove1 pos(2,3)*hscale2 pos(2,4)*vscale2])
set(h(3), 'Position',[ pos(3,1)+hmove1 pos(3,2)+vmove2 pos(3,3)*hscale2 pos(3,4)*vscale2])
set(h(4), 'Position',[ pos(4,1)+hmove2 pos(4,2)+vmove2 pos(4,3)*hscale2 pos(3,4)*vscale2])

%Text
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
text(0.400,0.05,'\itR^{forc}_{std}','HorizontalAlignment','center','FontSize',10,'FontName','Helvetica','Rotation',0);

%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr3(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr3(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr3(3,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); 
set(gca,'Visible','off','FontSize',9)
legend('CLM','VIC','PRMS',...
    'Orientation','vertical','Location',[0.825 0.375 0.05 0.075]);
legend('boxoff','Color','none');

%CLM  difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 1.2]); 
colorbar('horizontal','position',[0.100,0.075,0.25,0.035],...
   'XLim',[0 1.2],'XTick',0:0.3:1.2,'FontSize',10);
cbfreeze;

figfile=['./figure/paper/Paper_m10_fig4_interforce_RO_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);

%%
figure('Name','1-3 inter-forc-ET','NumberTitle','off','Units','inches','position',[1 7 6 5],'Color','w');
set(0,'DefaultAxesColorOrder',clr3);

subplot(2,2,1)
ax1 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax1, 'Visible','off','layer','top');
setm(ax1,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_force(:,:,1));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('CLM')

subplot(2,2,2)
ax2 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax2, 'Visible','off','layer','top');
setm(ax2,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_force(:,:,2));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('VIC')

subplot(2,2,3)
ax3 = usamap([lat1d(1) lat1d(end)],[lon1d(1)-360 lon1d(end)-360]);
set(ax3, 'Visible','off','layer','top');
setm(ax3,'MapProjection','mercator','FontSize',9,'grid','off','MeridianLabel','off','ParallelLabel','off');
pcolorm(lat2d, lon2d, R_eta_force(:,:,3));
plotm([S.Y],[S.X], 'k');
tightmap;
caxis([0 1.2]);freezeColors;
framem('FLineWidth',0.1,'FEdgeColor',[1 1 1])
title('PRMS')

subplot(2,2,4)
% for i=1:100
%     [n,xout] = hist(reshape(R_eta_force_ens(:,:,i,:),ylat_num*xlon_num,3),50);
%     plot(xout,n,'LineWidth',0.5);hold on
% end
[n,xout] = hist(reshape(R_eta_force,ylat_num*xlon_num,3),50);
plot(xout,n,'LineWidth',2);hold on
set(gca,'XLim',[0 1.25]);
xlabel('\itR^{forc}_{std}')
ylabel('count')
set(gca,'FontSize',9)

%adjust plot size & location
h   = get(gcf,'Children');
for i=1:length(h)
    pos(i,:) =get(h(i),'position');
end
hmove1 = +0.010;
hmove2 = -0.075;
vmove1 = +0.005;
vmove2 = -0.005;
hscale1=1.05;
vscale1=1.05;
hscale2=1.20;
vscale2=1.20;
%bottom row fromright to left
set(h(1), 'Position',[ pos(1,1)+hmove1 pos(1,2)+vmove1 pos(1,3)*hscale1 pos(1,4)*vscale1])
set(h(2), 'Position',[ pos(2,1)+hmove2 pos(2,2)+vmove1 pos(2,3)*hscale2 pos(2,4)*vscale2])
set(h(3), 'Position',[ pos(3,1)+hmove1 pos(3,2)+vmove2 pos(3,3)*hscale2 pos(3,4)*vscale2])
set(h(4), 'Position',[ pos(4,1)+hmove2 pos(4,2)+vmove2 pos(4,3)*hscale2 pos(3,4)*vscale2])

%legends 2
axes('unit','Normalized','position',[0 0 1 1],'visible','off');
dummy1 = plot(1,1,'-','LineWidth',2,'color',clr3(1,:));hold on
dummy2 = plot(1,1,'-','LineWidth',2,'color',clr3(2,:));
dummy3 = plot(1,1,'-','LineWidth',2,'color',clr3(3,:));
set(dummy1,'Visible','off'); set(dummy2,'Visible','off'); set(dummy3,'Visible','off'); 
set(gca,'Visible','off','FontSize',9)
legend('CLM','VIC','PRMS',...
    'Orientation','vertical','Location',[0.825 0.375 0.05 0.075]);
legend('boxoff','Color','none');

%CLM  difference colorbar, text on the plots
axes('Units','Normalized','position',[0 0 1 1],'visible','off');
caxis([0 1.2]); 
colorbar('horizontal','position',[0.100,0.075,0.25,0.035],...
   'XLim',[0 1.2],'XTick',0:0.3:1.2,'FontSize',10);
cbfreeze;

%save figure
figfile=['./figure/paper/Paper_m10_fig5_interforce_ET_us_' num2str(res) 'K'];
set(gcf,'PaperPositionMode','auto')
print('-dpng',figfile);

% profile viewer
