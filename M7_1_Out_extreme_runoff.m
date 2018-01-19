%Compute the following runoff metrics for each basin for each huc scale
%-annual mean
%-annual variance
%-annual runoff centroid
%- 2yr daily maximum runoff
%- 5yr daily maximum runoff
%-10yr daily maximum runoff
%-20yr daily maximum runoff
%-50yr daily maximum runoff
%- 2yr 5-day maximum runoff
%- 5yr 5-day maximum runoff
%-10yr 5-day maximum runoff
%-20yr 5-day maximum runoff
%-50yr 5-day maximum runoff
%-10yr daily minimum runoff
%-10yr 7-day minimum runoff

tic

close all
clear
clc
%%  Variable that can be changed
% -------------------------------------------------------------------------
% Forcing product list 
% -------------------------------------------------------------------------
force{1,1}='BCCA12K';
force{2,1}='BCSD12K';
force{3,1}='BCSDdisag12K';
force{4,1}='BCSAR12K';
force{5,1}='MAURER12K';
% -------------------------------------------------------------------------
% region list 
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
%byr=2001; bmon=10; bday=1;
eyr=2008; emon=9; eday=30;
% eyr=1999; emon=9; eday=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------Variables that cannot be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% water year array
% -------------------------------------------------------------------------
wyr = byr+1:1:eyr;
% -------------------------------------------------------------------------
% year starting in CLM output and ending (water year)
% -------------------------------------------------------------------------
swyr_clm = 1981;
% ewyr_clm = 1999;
ewyr_clm = 2008;
swyr_clm_wrf = 2002;
ewyr_clm_wrf = 2008;
% -------------------------------------------------------------------------
% date number for CLM output, and analysis
% -------------------------------------------------------------------------
% -- date number for monthly CLM output (h2, h3 and h4 file) 10/2000-9/2008
tnum_leap_clm = datenum(swyr_clm-1,10,1,12,0,0):1:datenum(ewyr_clm,9,30,12,0,0);
tvec_leap_clm = datevec(tnum_leap_clm);
%remove 2/29 in leap year
ii             = tvec_leap_clm(:,2)==2&tvec_leap_clm(:,3)==29;
tnum_clm       = tnum_leap_clm;
tvec_clm       = tvec_leap_clm;
tnum_clm(ii)   = [];
tvec_clm(ii,:) = [];
clear ii tvec_leap_clm tnum_leap_clm

% -- date number for monthly CLM output (h0 & h1 file) 10/2000-9/2008
ii = (tvec_clm(:,3)==15 & tvec_clm(:,4)==12);
tvec_clm_mon = tvec_clm(ii,:);
clear ii
% -------------------------------------------------------------------------
% date number for CLM output - ONLY WRF
% -------------------------------------------------------------------------
% -- date number for monthly CLM output (h2, h3 and h4 file) 10/2001-9/2008
tnum_leap_clm_wrf = datenum(swyr_clm_wrf-1,10,1,12,0,0):1:datenum(ewyr_clm_wrf,9,30,12,0,0);
tvec_leap_clm_wrf = datevec(tnum_leap_clm_wrf);
%remove 2/29 in leap year
ii             = tvec_leap_clm_wrf(:,2)==2&tvec_leap_clm_wrf(:,3)==29;
tnum_clm_wrf       = tnum_leap_clm_wrf;
tvec_clm_wrf       = tvec_leap_clm_wrf;
tnum_clm_wrf(ii)   = [];
tvec_clm_wrf(ii,:) = [];
clear ii tnum_leap_clm_wrf tvec_leap_clm_wrf

%{
% -------------------------------------------------------------------------
% date number for anaysis period - e.g. 10/1981-9/2008 or  10/2001-9/2008
% -------------------------------------------------------------------------
tnum_leap = datenum(byr,10,1,12,0,0):1:datenum(eyr,9,30,12,0,0);
tvec_leap = datevec(tnum_leap);
%remove 2/29 in leap year
ii         = tvec_leap(:,2)==2&tvec_leap(:,3)==29;
tnum       = tnum_leap;
tvec       = tvec_leap;
tnum(ii)   = [];
tvec(ii,:) = [];
clear ii tnum_leap tvec_leap

% number of days per month
numday= [31;30;31;31;28;31;30;31;30;31;31;30];
month = [10;11;12;1;2;3;4;5;6;7;8;9];
%}
%% READ data - CLM I/O and perform runoff process
% -------------------------------------------------------------------------
%  Process
% -------------------------------------------------------------------------
for r=8:8%1:length(region) %loop for region
    % -------------------------------------------------------------------------
    % HUC ID
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
    lon2d = netcdf2mat(region_domain_grid,'x'); lon2d(lon2d == -999) = NaN;
    lat2d = netcdf2mat(region_domain_grid,'y');
    lon1d = netcdf2mat(region_domain_grid,'lon');
    lat1d = netcdf2mat(region_domain_grid,'lat');
    ele   = netcdf2mat(region_domain_grid,'ele');
    huc2g = netcdf2mat(region_domain_grid,'huc2'); huc2g(huc2g~=huc)=NaN;
    huc4g = netcdf2mat(region_domain_grid,'huc4'); huc4g(huc2g~=huc)=NaN;
    huc8g = netcdf2mat(region_domain_grid,'huc8'); huc8g(huc2g~=huc)=NaN;
    xlon_num = length(lon1d);
    ylat_num = length(lat1d);
    clear region_domain_grid
    
    % HUC2, 4, and 8
    huc2ID = unique(huc2g);  huc2ID(isnan(huc2ID))=[];
    huc4ID = unique(huc4g);  huc4ID(isnan(huc4ID))=[];
    huc8ID = unique(huc8g);  huc8ID(isnan(huc8ID))=[];
    % -------------------------------------------------------------------------
    %  Memory allocation
    % -------------------------------------------------------------------------
    % Annual value for grid box
    % [lat x lon x wyr] -> [n x m x t x 2]
    ro_mean_yr = ones(ylat_num,xlon_num,length(wyr));
    ro_var_yr  = ones(ylat_num,xlon_num,length(wyr));
    ro_ct      = ones(ylat_num,xlon_num,length(wyr));
    RO_MAX_YR  = ones(ylat_num,xlon_num,length(wyr));
    RO5_MAX_YR = ones(ylat_num,xlon_num,length(wyr));
    RO_MIN_YR  = ones(ylat_num,xlon_num,length(wyr));
    RO7_MIN_YR = ones(ylat_num,xlon_num,length(wyr));
    
    % Annual value for huc2
    % [huc2 x wyr x prod] -> [n2 x t]
    ro_mean_yr_huc2   = ones(length(huc2ID),length(wyr));
    RO_VARI_YR_HUC2   = ones(length(huc2ID),length(wyr));
    RO_CT_HUC2        = ones(length(huc2ID),length(wyr));
    RO_MAX_YR_HUC2    = ones(length(huc2ID),length(wyr));
    RO5_MAX_YR_HUC2   = ones(length(huc2ID),length(wyr));
    RO_MIN_YR_HUC2    = ones(length(huc2ID),length(wyr));
    RO7_MIN_YR_HUC2   = ones(length(huc2ID),length(wyr));
    
    % Annual value for huc4
    % [huc4 x wyr x prod] -> [n4 x t]
    RO_MEAN_YR_HUC4   = ones(length(huc4ID),length(wyr));
    RO_VARI_YR_HUC4   = ones(length(huc4ID),length(wyr));
    RO_CT_HUC4        = ones(length(huc4ID),length(wyr));
    RO_MAX_YR_HUC4    = ones(length(huc4ID),length(wyr));
    RO5_MAX_YR_HUC4   = ones(length(huc4ID),length(wyr));
    RO_MIN_YR_HUC4    = ones(length(huc4ID),length(wyr));
    RO7_MIN_YR_HUC4   = ones(length(huc4ID),length(wyr));
    
    % Annual value for huc8
    % [huc8 x wyr x prod] -> [n8 x t]
    RO_MEAN_YR_HUC8   = ones(length(huc8ID),length(wyr));
    RO_VARI_YR_HUC8   = ones(length(huc8ID),length(wyr));
    RO_CT_HUC8        = ones(length(huc8ID),length(wyr));
    RO_MAX_YR_HUC8    = ones(length(huc8ID),length(wyr));
    RO5_MAX_YR_HUC8   = ones(length(huc8ID),length(wyr));
    RO_MIN_YR_HUC8    = ones(length(huc8ID),length(wyr));
    RO7_MIN_YR_HUC8   = ones(length(huc8ID),length(wyr));
    
    for f=1:length(force)    %loop for forcing data
        for m = 3:3%1:length(model)
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
            % -------------------------------------------------------------------------
            % CLM output to be read - from CLM history file
            % -------------------------------------------------------------------------
            if strcmp(model{m},'CLM')
                hclmvar{1,1}='QDRAI';   % sub-surface drainage (mm/s)
                hclmvar{2,1}='QOVER';   % surface runoff (mm/s)
            elseif strcmp (model{m},'VIC')
                hclmvar{1,1}='baseflow';   % sub-surface drainage (mm/s)
                hclmvar{2,1}='runoff';   % surface runoff (mm/s)
            elseif strcmp (model{m},'PRMS')
                hclmvar{1,1}='RUNOFF';   % total runoff (mm/s)
            end
            % --CLM daily history file
            ncname = [main_indir '/' force{f} '_' region{r} '/' force{f} '_' region{r} '.RUNOFF.all.nc'];
            for j=1:size(hclmvar,1)
                [var0,fillvalue]=netcdf2mat(ncname,hclmvar{j},'attname1','_FillValue');
                var0(var0==fillvalue)=NaN;
                hvar{j,1}=var0*60*60*24;
            end
            clear j var0 ncname hclmvar
            
            % The entire time series [lat x lon x 365*wyr]
            if strcmp(model{m},'CLM') || strcmp(model{m},'VIC')
                RO = hvar{1,1}(:,:,:)+hvar{2,1}(:,:,:);
            else
                RO = hvar{1,1};
            end
            
            RO_mask=mean(RO,3); RO_mask(~isnan(RO_mask))=1; %mask based on simulation results [lat x lon]
            RO=RO.*repmat(RO_mask,[1,1,size(RO,3)]);
            clear hvar
            
            % compute 5 days moving average - RO5
            ROa = permute(RO,[3,1,2]); % Make time dimension first one
            ROb = reshape(ROa,size(ROa,1),ylat_num*xlon_num); % Reshape to [365*wyr x lat*lon)
            RO5a = movingave(ROb,5); clear ROa ROb
            RO5b = reshape(RO5a,size(RO5a,1),ylat_num,xlon_num);
            ro5 = permute(RO5b,[2,3,1]);clear RO5a RO5b
            
            % compute 7 days moving average - RO7
            ROa = permute(RO,[3,1,2]); % Make time dimension first one
            ROb = reshape(ROa,size(ROa,1),ylat_num*xlon_num); % Reshape to [365*wyr x lat*lon)
            RO7a = movingave(ROb,7); clear ROa ROb
            RO7b = reshape(RO7a,size(RO7a,1),ylat_num,xlon_num);
            ro7 = permute(RO7b,[2,3,1]);clear RO7a RO7b
            
            for k=1:length(wyr)
                %---- Find row indice for current WY
                if ~strcmp(force{f},'WRF') && ~strcmp(force{f},'WRFag')
                    ii=find((tvec_clm(:,1)==wyr(k)-1 & tvec_clm(:,2)>=10) | (tvec_clm(:,1)==wyr(k) & tvec_clm(:,2)<=9));
                else
                    ii=find((tvec_clm_wrf(:,1)==wyr(k)-1 & tvec_clm_wrf(:,2)>=10) | (tvec_clm_wrf(:,1)==wyr(k) & tvec_clm_wrf(:,2)<=9));
                end
                % Daily runoff [mm/day] per year
                % ro_wyr [lat x lon x day_per_yr]
                ro_wyr    = RO(:,:,ii);
                ro7_wyr   = ro7(:,:,ii);
                ro5_wyr   = ro5(:,:,ii);
                clear ii
                
                % HUC2 average
                % Memory allocation
                % Annual  value - [HUC x day_per_yr]
                ro_wyr_huc2  = ones(length(huc2ID),365)*NaN;
                ro5_wyr_huc2 = ones(length(huc2ID),365)*NaN;
                ro7_wyr_huc2 = ones(length(huc2ID),365)*NaN;
                for l = 1:length(huc2ID)
                    %make mask for current huc id
                    huc2_mask = huc2g;  huc2_mask(~isnan(huc2g))=1;
                    huc2_mask(huc2g~=huc2ID(l) | isnan(lon2d) | isnan(RO_mask))=NaN;
                    
                    % 3D mask (nlat x nlon x 365]
                    huc2_mask_wyr = repmat(huc2_mask,[1,1,365]);
                    
                    ro_wyr_huc2a  = reshape(huc2_mask_wyr.*ro_wyr,ylat_num*xlon_num,365);
                    ro5_wyr_huc2a = reshape(huc2_mask_wyr.*ro5_wyr,ylat_num*xlon_num,365);
                    ro7_wyr_huc2a = reshape(huc2_mask_wyr.*ro7_wyr,ylat_num*xlon_num,365);
                    
                    ro_wyr_huc2(l,:)  = mean(ro_wyr_huc2a (~isnan(ro_wyr_huc2a (:,100)), :),1);
                    ro5_wyr_huc2(l,:) = mean(ro5_wyr_huc2a(~isnan(ro5_wyr_huc2a(:,100)), :),1);
                    ro7_wyr_huc2(l,:) = mean(ro7_wyr_huc2a(~isnan(ro7_wyr_huc2a(:,100)), :),1);
                    
                    clear huc2_mask huc2_mask_wyr ro_wyr_huc2a ro5_wyr_huc2a ro7_wyr_huc2a
                end
                clear l
                %HUC4 average
                % Memory allocation
                % Annual  value - [HUC x day_per_yr]
                ro_wyr_huc4  = ones(length(huc4ID),365)*NaN;
                ro5_wyr_huc4 = ones(length(huc4ID),365)*NaN;
                ro7_wyr_huc4 = ones(length(huc4ID),365)*NaN;
                for l = 1:length(huc4ID)
                    %make mask for current huc id
                    huc4_mask = huc4g;  huc4_mask(~isnan(huc4g))=1;
                    huc4_mask(huc4g~=huc4ID(l) | isnan(lon2d) | isnan(RO_mask))=NaN;
                    
                    % 3D mask (nlat x nlon x 365]
                    huc4_mask_wyr = repmat(huc4_mask,[1,1,365]);
                    
                    ro_wyr_huc4a  = reshape(huc4_mask_wyr.*ro_wyr,ylat_num*xlon_num,365);
                    ro5_wyr_huc4a = reshape(huc4_mask_wyr.*ro5_wyr,ylat_num*xlon_num,365);
                    ro7_wyr_huc4a = reshape(huc4_mask_wyr.*ro7_wyr,ylat_num*xlon_num,365);
                    
                    ro_wyr_huc4(l,:)  = mean(ro_wyr_huc4a (~isnan(ro_wyr_huc4a (:,100)), :),1);
                    ro5_wyr_huc4(l,:) = mean(ro5_wyr_huc4a(~isnan(ro5_wyr_huc4a(:,100)), :),1);
                    ro7_wyr_huc4(l,:) = mean(ro7_wyr_huc4a(~isnan(ro7_wyr_huc4a(:,100)), :),1);
                    clear huc4_mask huc4_mask_wyr ro_wyr_huc4a ro5_wyr_huc4a ro7_wyr_huc4a
                end
                clear ii jj l
                %HUC8 average
                % Memory allocation
                % Annual  value - [HUC x day_per_yr]
                ro_wyr_huc8  = ones(length(huc8ID),365)*NaN;
                ro5_wyr_huc8 = ones(length(huc8ID),365)*NaN;
                ro7_wyr_huc8 = ones(length(huc8ID),365)*NaN;
                for l = 1:length(huc8ID)
                    %make mask for current huc id
                    huc8_mask = huc8g;  huc8_mask(~isnan(huc8g))=1;
                    huc8_mask(huc8g~=huc8ID(l) | isnan(lon2d) | isnan(RO_mask))=NaN;
                    
                    % 3D mask (nlat x nlon x 365]
                    huc8_mask_wyr = repmat(huc8_mask,[1,1,365]);
                    
                    ro_wyr_huc8a  = reshape(huc8_mask_wyr.*ro_wyr,ylat_num*xlon_num,365);
                    ro5_wyr_huc8a = reshape(huc8_mask_wyr.*ro5_wyr,ylat_num*xlon_num,365);
                    ro7_wyr_huc8a = reshape(huc8_mask_wyr.*ro7_wyr,ylat_num*xlon_num,365);
                    
                    ro_wyr_huc8(l,:)  = mean(ro_wyr_huc8a (~isnan(ro_wyr_huc8a (:,100)), :),1);
                    ro5_wyr_huc8(l,:) = mean(ro5_wyr_huc8a(~isnan(ro5_wyr_huc8a(:,100)), :),1);
                    ro7_wyr_huc8(l,:) = mean(ro7_wyr_huc8a(~isnan(ro7_wyr_huc8a(:,100)), :),1);
                    clear huc8_mask huc8_mask_wyr ro_wyr_huc8a ro5_wyr_huc8a ro7_wyr_huc8a
                end
                clear ii jj l
                
                %Days from October 1st
                daynum = linspace(1,365,365);
                days = repmat(reshape(daynum,1,1,length(daynum)),[ylat_num,xlon_num,1]);
                
                % -- Compute annual metric
                %Compute stat per grid box
                ro_ct(:,:,k)      = sum(ro_wyr.*days,3)./sum(ro_wyr,3);
                ro_mean_yr(:,:,k) = mean(ro_wyr,3);
                ro_var_yr(:,:,k)  = std(ro_wyr,0,3);
                
                RO_MAX_YR(:,:,k)  = max(ro_wyr,[],3);
                RO5_MAX_YR(:,:,k) = nanmax(ro5_wyr,[],3);
                RO_MIN_YR(:,:,k)  = min(ro_wyr,[],3);
                RO7_MIN_YR(:,:,k) = nanmin(ro7_wyr,[],3);
                
                clear ro_wyr ro5_wyr ro7_wyr days
                
                %Compute stat per HUC2
                days = repmat(1:1:365,[length(huc2ID),1]);
                RO_CT_HUC2(:,k)      = sum(ro_wyr_huc2.*days,2)./sum(ro_wyr_huc2,2);
                ro_mean_yr_huc2(:,k) = mean(ro_wyr_huc2,2);
                RO_VARI_YR_HUC2(:,k) = std(ro_wyr_huc2,0,2);
                
                RO_MAX_YR_HUC2(:,k)  = max(ro_wyr_huc2,[],2);
                RO5_MAX_YR_HUC2(:,k) = nanmax(ro5_wyr_huc2,[],2);
                RO_MIN_YR_HUC2(:,k)  = min(ro_wyr_huc2,[],2);
                RO7_MIN_YR_HUC2(:,k) = nanmin(ro7_wyr_huc2,[],2);
                clear days
                
                %Compute stat per HUC4
                days = repmat(1:1:365,[length(huc4ID),1]);
                RO_CT_HUC4(:,k)      = sum(ro_wyr_huc4.*days,2)./sum(ro_wyr_huc4,2);
                RO_MEAN_YR_HUC4(:,k) = mean(ro_wyr_huc4,2);
                RO_VARI_YR_HUC4(:,k) = std(ro_wyr_huc4,0,2);
                
                RO_MAX_YR_HUC4(:,k)  = max(ro_wyr_huc4,[],2);
                RO5_MAX_YR_HUC4(:,k) = nanmax(ro5_wyr_huc4,[],2);
                RO_MIN_YR_HUC4(:,k)  = min(ro_wyr_huc4,[],2);
                RO7_MIN_YR_HUC4(:,k) = nanmin(ro7_wyr_huc4,[],2);
                clear days
                
                %Compute stat per HUC8
                days = repmat(1:1:365,[length(huc8ID),1]);
                RO_CT_HUC8(:,k)      = sum(ro_wyr_huc8.*days,2)./sum(ro_wyr_huc8,2);
                RO_MEAN_YR_HUC8(:,k) = mean(ro_wyr_huc8,2);
                RO_VARI_YR_HUC8(:,k) = std(ro_wyr_huc8,0,2);
                
                RO_MAX_YR_HUC8(:,k)  = max(ro_wyr_huc8,[],2);
                RO5_MAX_YR_HUC8(:,k) = nanmax(ro5_wyr_huc8,[],2);
                RO_MIN_YR_HUC8(:,k)  = min(ro_wyr_huc8,[],2);
                RO7_MIN_YR_HUC8(:,k) = nanmin(ro7_wyr_huc8,[],2);
                clear days
            end
            clear RO RO7 RO5 i k l
            
            %% Grid box analysis
            % Annual mean and variability
            RO_CT_C      = squeeze( mean(ro_ct,3) );
            RO_MEAN_YR_C = squeeze( mean(ro_mean_yr,3) );
            RO_VARI_YR_C = squeeze( mean(ro_var_yr,3) );
            
            %Frequency analysis
            %Ascending order and compute CDF which is nonexceedance probability
            %--------------------------------
            %Low flow analysis - 7 day event
            %--------------------------------
            % compute 7Q10 - RO7_MIN10
            RO7_MIN10b=prctile(RO7_MIN_YR,10,3);
            RO7_MIN10 = squeeze( RO7_MIN10b );
            clear RO7_MIN10b
            
            %--------------------------------
            %Low flow analysis - 1 day event
            %--------------------------------
            RO_MIN10b=prctile(RO_MIN_YR,10,3);
            RO_MIN10 = squeeze( RO_MIN10b );
            clear RO_MIN10b
            
            %Decending order and compute CDF which is exceedance probability
            %--------------------------------
            %High flow analysis - 5 day event
            %--------------------------------
            RO5_MAX2b=prctile(RO5_MAX_YR,50,3);
            RO5_MAX2 = squeeze( RO5_MAX2b );
            clear RO5_MAX2b
            
            % compute 5 year 5 day events -RO5_MAX5
            RO5_MAX5b=prctile(RO5_MAX_YR,80,3);
            RO5_MAX5 = squeeze( RO5_MAX5b );
            clear RO5_MAX5b
            
            % compute 10 year 5 day events -RO5_MAX10
            RO5_MAX10b=prctile(RO5_MAX_YR,90,3);
            RO5_MAX10 = squeeze( RO5_MAX10b );
            clear RO5_MAX10b
            
            % compute 20 year 5 day events -RO5_MAX20
            RO5_MAX20b=prctile(RO5_MAX_YR,98,3);
            RO5_MAX20 = squeeze( RO5_MAX20b );
            clear RO5_MAX20b
            
            % compute 50 year 5 day events -RO5_MAX50
            RO5_MAX50b=prctile(RO5_MAX_YR,98,3);
            RO5_MAX50 = squeeze( RO5_MAX50b );
            clear RO5_MAX50b
            
            %--------------------------------
            %High flow analysis - 1 day event
            %--------------------------------
            RO_MAX2b=prctile(RO_MAX_YR,50,3);
            RO_MAX2 = squeeze( RO_MAX2b );
            clear RO_MAX2b
            
            % compute 5 year one day events -RO_MAX5
            RO_MAX5b=prctile(RO_MAX_YR,80,3);
            RO_MAX5 = squeeze( RO_MAX5b );
            clear RO_MAX5b
            
            % compute 20 year one day events -RO_MAX10
            RO_MAX10b=prctile(RO_MAX_YR,90,3);
            RO_MAX10 = squeeze( RO_MAX10b );
            clear RO_MAX10b
            
            % compute 20 year one day events -RO_MAX20
            RO_MAX20b=prctile(RO_MAX_YR,95,3);
            RO_MAX20 = squeeze(RO_MAX20b );
            clear RO_MAX20b
            
            % compute 50 year one day events -RO_MAX50
            RO_MAX50b=prctile(RO_MAX_YR,98,3);
            RO_MAX50 = squeeze( RO_MAX50b );
            clear RO_MAX50b
            
            %Change nan to -999
            RO_MEAN_YR_C(isnan(RO_MEAN_YR_C)) = -999;
            RO_VARI_YR_C(isnan(RO_VARI_YR_C)) = -999;
            RO_CT_C(isnan(RO_CT_C)) = -999;
            RO7_MIN10(isnan(RO7_MIN10)) = -999;
            RO_MIN10(isnan(RO_MIN10)) = -999;
            RO5_MAX2(isnan(RO5_MAX2)) = -999;
            RO5_MAX5(isnan(RO5_MAX5)) = -999;
            RO5_MAX10(isnan(RO5_MAX10)) = -999;
            RO5_MAX20(isnan(RO5_MAX20)) = -999;
            RO5_MAX50(isnan(RO5_MAX50)) = -999;
            RO_MAX2(isnan(RO_MAX2)) = -999;
            RO_MAX5(isnan(RO_MAX5)) = -999;
            RO_MAX10(isnan(RO_MAX10)) = -999;
            RO_MAX20(isnan(RO_MAX20)) = -999;
            RO_MAX50(isnan(RO_MAX50)) = -999;
            
            %% Output in netCDF
            %---- Open netCDF file.
%             outnc = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF.nc'];
            outnc = [main_indir '/' force{f} '_' region{r} '/Extreme_RUNOFF.nc'];
            ncid = netcdf.create(outnc,'CLOBBER');
            
            %---- Define the dimensions of the variable
            ylat_dimid = netcdf.defDim( ncid, 'lat'  ,ylat_num);
            xlon_dimid = netcdf.defDim( ncid, 'lon'  ,xlon_num);
            
            %---- Define a new variable in the file.
            ro_mean_varid      = netcdf.defVar ( ncid, 'ro_mean',     'float', [xlon_dimid,ylat_dimid] );
            ro_var_varid       = netcdf.defVar ( ncid, 'ro_var',      'float', [xlon_dimid,ylat_dimid] );
            ro7_10yr_min_varid = netcdf.defVar ( ncid, 'ro7_10yr_min','float', [xlon_dimid,ylat_dimid] );
            ro_10yr_min_varid  = netcdf.defVar ( ncid, 'ro_10yr_min', 'float', [xlon_dimid,ylat_dimid] );
            ro5_02yr_max_varid = netcdf.defVar ( ncid, 'ro5_02yr_max','float', [xlon_dimid,ylat_dimid] );
            ro5_05yr_max_varid = netcdf.defVar ( ncid, 'ro5_05yr_max','float', [xlon_dimid,ylat_dimid] );
            ro5_10yr_max_varid = netcdf.defVar ( ncid, 'ro5_10yr_max','float', [xlon_dimid,ylat_dimid] );
            ro5_20yr_max_varid = netcdf.defVar ( ncid, 'ro5_20yr_max','float', [xlon_dimid,ylat_dimid] );
            ro5_50yr_max_varid = netcdf.defVar ( ncid, 'ro5_50yr_max','float', [xlon_dimid,ylat_dimid] );
            ro_02yr_max_varid  = netcdf.defVar ( ncid, 'ro_02yr_max', 'float', [xlon_dimid,ylat_dimid] );
            ro_05yr_max_varid  = netcdf.defVar ( ncid, 'ro_05yr_max', 'float', [xlon_dimid,ylat_dimid] );
            ro_10yr_max_varid  = netcdf.defVar ( ncid, 'ro_10yr_max', 'float', [xlon_dimid,ylat_dimid] );
            ro_20yr_max_varid  = netcdf.defVar ( ncid, 'ro_20yr_max', 'float', [xlon_dimid,ylat_dimid] );
            ro_50yr_max_varid  = netcdf.defVar ( ncid, 'ro_50yr_max', 'float', [xlon_dimid,ylat_dimid] );
            
            %---- Define the coordinate variables.
            lat_d_varid = netcdf.defVar ( ncid, 'lat', 'float', ylat_dimid );
            lon_d_varid = netcdf.defVar ( ncid, 'lon', 'float', xlon_dimid );
            
            %---- Leave define mode and enter data mode to write data.
            netcdf.endDef (ncid);
            
            %---- Write the coordinate variable data.
            netcdf.putVar ( ncid, lat_d_varid, lat1d );
            netcdf.putVar ( ncid, lon_d_varid, lon1d );
            
            %---- Write data to variable
            netcdf.putVar ( ncid, ro_mean_varid,      permute(RO_MEAN_YR_C,[2,1]));
            netcdf.putVar ( ncid, ro_var_varid,       permute(RO_VARI_YR_C,[2,1]));
            netcdf.putVar ( ncid, ro7_10yr_min_varid, permute(RO7_MIN10,[2,1]));
            netcdf.putVar ( ncid, ro_10yr_min_varid,  permute(RO_MIN10,[2,1]));
            netcdf.putVar ( ncid, ro5_02yr_max_varid, permute(RO5_MAX2,[2,1]));
            netcdf.putVar ( ncid, ro5_05yr_max_varid, permute(RO5_MAX5,[2,1]));
            netcdf.putVar ( ncid, ro5_10yr_max_varid, permute(RO5_MAX10,[2,1]));
            netcdf.putVar ( ncid, ro5_20yr_max_varid, permute(RO5_MAX20,[2,1]));
            netcdf.putVar ( ncid, ro5_50yr_max_varid, permute(RO5_MAX50,[2,1]));
            netcdf.putVar ( ncid, ro_02yr_max_varid,  permute(RO_MAX2,[2,1]));
            netcdf.putVar ( ncid, ro_05yr_max_varid,  permute(RO_MAX5,[2,1]));
            netcdf.putVar ( ncid, ro_10yr_max_varid,  permute(RO_MAX10,[2,1]));
            netcdf.putVar ( ncid, ro_20yr_max_varid,  permute(RO_MAX20,[2,1]));
            netcdf.putVar ( ncid, ro_50yr_max_varid,  permute(RO_MAX50,[2,1]));
            
            %---- Re-enter define mode.
            netcdf.reDef(ncid);
            %
            %---- Create an attribute associated with the variable.
            
            %- Mean
            netcdf.putAtt(ncid, ro_mean_varid, 'units','mm/day');
            netcdf.putAtt(ncid, ro_mean_varid, 'long_name','Annual mean runoff');
            %- Variability
            netcdf.putAtt(ncid, ro_var_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_var_varid ,'long_name','Annual runoff varialce');
            %- low flow frequency
            netcdf.putAtt(ncid, ro7_10yr_min_varid, 'units','mm/day');
            netcdf.putAtt(ncid, ro7_10yr_min_varid, 'long_name','Minimum 7days 10yr runoff');
            %-
            netcdf.putAtt(ncid, ro_10yr_min_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_10yr_min_varid ,'long_name','Minimum daily 10yr runoff');
            %- high flow frequency
            netcdf.putAtt(ncid, ro5_02yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro5_02yr_max_varid ,'long_name','Maximum 5days 2yr runoff');
            %-
            netcdf.putAtt(ncid, ro5_05yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro5_05yr_max_varid ,'long_name','Maximum 5days 5yr runoff');
            %-
            netcdf.putAtt(ncid, ro5_10yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro5_10yr_max_varid ,'long_name','Maximum 5days 10yr runoff');
            %-
            netcdf.putAtt(ncid, ro5_20yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro5_20yr_max_varid ,'long_name','Maximum 5days 20yr runoff');
            %-
            netcdf.putAtt(ncid, ro5_50yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro5_50yr_max_varid ,'long_name','Maximum 5days 50yr runoff');
            %-
            netcdf.putAtt(ncid, ro_02yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_02yr_max_varid ,'long_name','Maximum daily 2yr runoff');
            %-
            netcdf.putAtt(ncid, ro_05yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_05yr_max_varid ,'long_name','Maximum daily 5yr runoff');
            %-
            netcdf.putAtt(ncid, ro_10yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_10yr_max_varid ,'long_name','Maximum daily 10yr runoff');
            %-
            netcdf.putAtt(ncid, ro_20yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_20yr_max_varid ,'long_name','Maximum daily 20yr runoff');
            %-
            netcdf.putAtt(ncid, ro_50yr_max_varid ,'units','mm/day');
            netcdf.putAtt(ncid, ro_50yr_max_varid ,'long_name','Maximum daily 50yr runoff');
            
            %-Global
            gattid=netcdf.getConstant('NC_GLOBAL');
            netcdf.putAtt(ncid,gattid,'projection','lat/lon grid');
            netcdf.putAtt(ncid,gattid,'x-resolution',[num2str(res) ' km']);
            netcdf.putAtt(ncid,gattid,'y-resolution',[num2str(res) ' km']);
            netcdf.putAtt(ncid,gattid,'matlab file',[pwd '/M12_Out_annual_runoff_stat_28yr.m']);
            netcdf.close(ncid);
            %
            %---- Add a fill (missing) value for each grid. use NCO
            %command ncatted
            system(['ncatted -h -a _FillValue,ro_mean,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_var,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro7_10yr_min,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_10yr_min,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro5_02yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro5_05yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro5_10yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro5_20yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro5_50yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_02yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_05yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_10yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_20yr_max,o,f,-999. ' outnc]);
            system(['ncatted -h -a _FillValue,ro_50yr_max,o,f,-999. ' outnc]);
            
            %% HUC analysis
            % HUC2 ID and huc ave. elevation
            % Memory allocation [HUC x 1]
            ID_huc2  = ones(length(huc2ID),1)*NaN;
            elev_huc2 = ones(length(huc2ID),1)*NaN;
            for l = 1:length(huc2ID)
                %make mask for current huc id
                huc2_mask = huc2g;  huc2_mask(~isnan(huc2g))=1;
                huc2_mask(huc2g~=huc2ID(l))=NaN;
                
                elev_huc2a = reshape(huc2_mask.*ele,   ylat_num*xlon_num,1);
                ID_huc2a  = reshape(huc2_mask.*huc2g, ylat_num*xlon_num,1);
                
                elev_huc2(l,:) = mean(elev_huc2a(~isnan(elev_huc2a (:,1)), :),1);
                ID_huc2(l,:)  = mean(ID_huc2a (~isnan(ID_huc2a(:,1)),   :),1);
                
                clear huc2_mask elev_huc2a ID_huc2a
            end
            clear l
            %HUC4 ID and huc ave. elevation
            % Memory allocation [HUC x 1]
            elev_huc4 = ones(length(huc4ID),1)*NaN;
            ID_huc4  = ones(length(huc4ID),1)*NaN;
            for l = 1:length(huc4ID)
                %make mask for current huc id
                huc4_mask = huc4g;  huc4_mask(~isnan(huc4g))=1;
                huc4_mask(huc4g~=huc4ID(l))=NaN;
                
                elev_huc4a  = reshape(huc4_mask.*ele,  ylat_num*xlon_num ,1);
                ID_huc4a  = reshape(huc4_mask.*huc4g, ylat_num*xlon_num,1);
                
                elev_huc4(l,:) = mean(elev_huc4a(~isnan(elev_huc4a (:,1)), :),1);
                ID_huc4(l,:)  = mean(ID_huc4a (~isnan(ID_huc4a(:,1)),   :),1);
                clear huc4_mask elev_huc4a ID_huc4a
            end
            clear l
            %HUC8 ID and huc ave. elevation
            % Memory allocation [HUC x 1]
            elev_huc8 = ones(length(huc8ID),1)*NaN;
            ID_huc8  = ones(length(huc8ID),1)*NaN;
            for l = 1:length(huc8ID)
                %make mask for current huc id
                huc8_mask = huc8g;  huc8_mask(~isnan(huc8g))=1;
                huc8_mask(huc8g~=huc8ID(l))=NaN;
                
                elev_huc8a  = reshape(huc8_mask.*ele,  ylat_num*xlon_num,1);
                ID_huc8a  = reshape(huc8_mask.*huc8g, ylat_num*xlon_num,1);
                
                elev_huc8(l,:) = mean(elev_huc8a(~isnan(elev_huc8a (:,1)), :),1);
                ID_huc8(l,:)  = mean(ID_huc8a (~isnan(ID_huc8a(:,1)),   :),1);
                
                clear huc8_mask elev_huc8a ID_huc8a
            end
            clear l
            
            % Annual mean and variability
            RO_CT_HUC2_C = squeeze( mean(RO_CT_HUC2,2) );
            ro_mean_huc2 = squeeze( mean(ro_mean_yr_huc2,2) );
            ro_var_huc2  = squeeze( mean(RO_VARI_YR_HUC2,2) );
            
            RO_CT_HUC4_C = squeeze( mean(RO_CT_HUC4,2) );
            ro_mean_huc4 = squeeze( mean(RO_MEAN_YR_HUC4,2) );
            ro_var_huc4  = squeeze( mean(RO_VARI_YR_HUC4,2) );
            
            RO_CT_HUC8_C = squeeze( mean(RO_CT_HUC8,2) );
            ro_mean_huc8 = squeeze( mean(RO_MEAN_YR_HUC8,2) );
            ro_var_huc8  = squeeze( mean(RO_VARI_YR_HUC8,2) );
            
            %Frequency analysis
            %Ascending order and compute CDF which is nonexceedance probability
            %--------------------------------
            %Low flow analysis - 7 day event
            %--------------------------------
            % compute 7Q10 - RO7_MIN10
            RO7_MIN10b=prctile(RO7_MIN_YR_HUC2,10,2);
            ro7_10yr_min_huc2 = squeeze( RO7_MIN10b );
            clear RO7_MIN10b
            RO7_MIN10b=prctile(RO7_MIN_YR_HUC4,10,2);
            ro7_10yr_min_huc4 = squeeze( RO7_MIN10b );
            clear RO7_MIN10b
            RO7_MIN10b=prctile(RO7_MIN_YR_HUC8,10,2);
            ro7_10yr_min_huc8 = squeeze( RO7_MIN10b );
            clear RO7_MIN10b
            %--------------------------------
            %Low flow analysis - 1 day event
            %--------------------------------
            RO_MIN10b=prctile(RO_MIN_YR_HUC2,10,2);
            ro_10yr_min_huc2 = squeeze( RO_MIN10b );
            clear RO_MIN10b
            RO_MIN10b=prctile(RO_MIN_YR_HUC4,10,2);
            ro_10yr_min_huc4 = squeeze( RO_MIN10b );
            clear RO_MIN10b
            RO_MIN10b=prctile(RO_MIN_YR_HUC8,10,2);
            ro_10yr_min_huc8 = squeeze( RO_MIN10b );
            clear RO_MIN10b
            
            %Decending order and compute CDF which is exceedance probability
            %--------------------------------
            %High flow analysis - 5 day event
            %--------------------------------
            RO5_MAX2b=prctile(RO5_MAX_YR_HUC2,50,2);
            ro5_02yr_max_huc2 = squeeze( RO5_MAX2b );
            clear RO5_MAX2b
            RO5_MAX2b=prctile(RO5_MAX_YR_HUC4,50,2);
            ro5_02yr_max_huc4 = squeeze( RO5_MAX2b );
            clear RO5_MAX2b
            RO5_MAX2b=prctile(RO5_MAX_YR_HUC8,50,2);
            ro5_02yr_max_huc8 = squeeze( RO5_MAX2b );
            clear RO5_MAX2b
            
            % compute 5 year 5 day events -ro5_05yr_max
            RO5_MAX5b=prctile(RO5_MAX_YR_HUC2,80,2);
            ro5_05yr_max_huc2 = squeeze( RO5_MAX5b );
            clear RO5_MAX5b
            RO5_MAX5b=prctile(RO5_MAX_YR_HUC4,80,2);
            ro5_05yr_max_huc4 = squeeze( RO5_MAX5b );
            clear RO5_MAX5b
            RO5_MAX5b=prctile(RO5_MAX_YR_HUC8,80,2);
            ro5_05yr_max_huc8 = squeeze( RO5_MAX5b );
            clear RO5_MAX5b
            
            % compute 10 year 5 day events -ro5_10yr_max
            RO5_MAX10b=prctile(RO5_MAX_YR_HUC2,90,2);
            ro5_10yr_max_huc2 = squeeze( RO5_MAX10b );
            clear RO5_MAX10b
            RO5_MAX10b=prctile(RO5_MAX_YR_HUC4,90,2);
            ro5_10yr_max_huc4 = squeeze( RO5_MAX10b );
            clear RO5_MAX10b
            RO5_MAX10b=prctile(RO5_MAX_YR_HUC8,90,2);
            ro5_10yr_max_huc8 = squeeze( RO5_MAX10b );
            clear RO5_MAX10b
            
            % compute 20 year 5 day events -ro5_20yr_max
            RO5_MAX20b=prctile(RO5_MAX_YR_HUC2,98,2);
            ro5_20yr_max_huc2 = squeeze( RO5_MAX20b );
            clear RO5_MAX20b
            RO5_MAX20b=prctile(RO5_MAX_YR_HUC4,98,2);
            ro5_20yr_max_huc4 = squeeze( RO5_MAX20b );
            clear RO5_MAX20b
            RO5_MAX20b=prctile(RO5_MAX_YR_HUC8,98,2);
            ro5_20yr_max_huc8 = squeeze( RO5_MAX20b );
            clear RO5_MAX20b
            
            % compute 50 year 5 day events -ro5_50yr_max
            RO5_MAX50b=prctile(RO5_MAX_YR_HUC2,98,2);
            ro5_50yr_max_huc2 = squeeze( RO5_MAX50b );
            clear RO5_MAX50b
            RO5_MAX50b=prctile(RO5_MAX_YR_HUC4,98,2);
            ro5_50yr_max_huc4 = squeeze( RO5_MAX50b );
            clear RO5_MAX50b
            RO5_MAX50b=prctile(RO5_MAX_YR_HUC8,98,2);
            ro5_50yr_max_huc8 = squeeze( RO5_MAX50b );
            clear RO5_MAX50b
            
            %--------------------------------
            %High flow analysis - 1 day event
            %--------------------------------
            RO_MAX2b=prctile(RO_MAX_YR_HUC2,50,2);
            ro_02yr_max_huc2 = squeeze( RO_MAX2b );
            clear RO_MAX2b
            RO_MAX2b=prctile(RO_MAX_YR_HUC4,50,2);
            ro_02yr_max_huc4 = squeeze( RO_MAX2b );
            clear RO_MAX2b
            RO_MAX2b=prctile(RO_MAX_YR_HUC8,50,2);
            ro_02yr_max_huc8 = squeeze( RO_MAX2b );
            clear RO_MAX2b
            
            % compute 5 year one day events -RO_MAX5
            RO_MAX5b=prctile(RO_MAX_YR_HUC2,80,2);
            ro_05yr_max_huc2 = squeeze( RO_MAX5b );
            clear RO_MAX5b
            RO_MAX5b=prctile(RO_MAX_YR_HUC4,80,2);
            ro_05yr_max_huc4 = squeeze( RO_MAX5b );
            clear RO_MAX5b
            RO_MAX5b=prctile(RO_MAX_YR_HUC8,80,2);
            ro_05yr_max_huc8 = squeeze( RO_MAX5b );
            clear RO_MAX5b
            
            % compute 20 year one day events -RO_MAX10
            RO_MAX10b=prctile(RO_MAX_YR_HUC2,90,2);
            ro_10yr_max_huc2 = squeeze( RO_MAX10b );
            clear RO_MAX10b
            RO_MAX10b=prctile(RO_MAX_YR_HUC4,90,2);
            ro_10yr_max_huc4 = squeeze( RO_MAX10b );
            clear RO_MAX10b
            RO_MAX10b=prctile(RO_MAX_YR_HUC8,90,2);
            ro_10yr_max_huc8 = squeeze( RO_MAX10b );
            clear RO_MAX10b
            
            % compute 20 year one day events -RO_MAX20
            RO_MAX20b=prctile(RO_MAX_YR_HUC2,95,2);
            ro_20yr_max_huc2 = squeeze(RO_MAX20b );
            clear RO_MAX20b
            RO_MAX20b=prctile(RO_MAX_YR_HUC4,95,2);
            ro_20yr_max_huc4 = squeeze(RO_MAX20b );
            clear RO_MAX20b
            RO_MAX20b=prctile(RO_MAX_YR_HUC8,95,2);
            ro_20yr_max_huc8 = squeeze(RO_MAX20b );
            clear RO_MAX20b
            
            % compute 50 year one day events -RO_MAX50
            RO_MAX50b=prctile(RO_MAX_YR_HUC2,98,2);
            ro_50yr_max_huc2 = squeeze( RO_MAX50b );
            clear RO_MAX50b
            RO_MAX50b=prctile(RO_MAX_YR_HUC4,98,2);
            ro_50yr_max_huc4 = squeeze( RO_MAX50b );
            clear RO_MAX50b
            RO_MAX50b=prctile(RO_MAX_YR_HUC8,98,2);
            ro_50yr_max_huc8 = squeeze( RO_MAX50b );
            clear RO_MAX50b
            
            %% Matlab binary output for HUC aggregated values
%             huc2_mat = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF_huc2.mat'];
%             huc4_mat = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF_huc4.mat'];
%             huc8_mat = [main_indir '/' force{f} '_' region{r} '/calib/Extreme_RUNOFF_huc8.mat'];
            huc2_mat = [main_indir '/' force{f} '_' region{r} '/Extreme_RUNOFF_huc2.mat'];
            huc4_mat = [main_indir '/' force{f} '_' region{r} '/Extreme_RUNOFF_huc4.mat'];
            huc8_mat = [main_indir '/' force{f} '_' region{r} '/Extreme_RUNOFF_huc8.mat']; 
            
            save(huc2_mat,'ID_huc2','elev_huc2','ro_mean_huc2','ro_var_huc2','ro_10yr_min_huc2','ro7_10yr_min_huc2',...
                'ro_02yr_max_huc2','ro_05yr_max_huc2','ro_10yr_max_huc2','ro_20yr_max_huc2','ro_50yr_max_huc2',...
                'ro5_02yr_max_huc2','ro5_05yr_max_huc2','ro5_10yr_max_huc2','ro5_20yr_max_huc2','ro5_50yr_max_huc2');
            
            save(huc4_mat,'ID_huc4','elev_huc4','ro_mean_huc4','ro_var_huc4','ro_10yr_min_huc4','ro7_10yr_min_huc4',...
                'ro_02yr_max_huc4','ro_05yr_max_huc4','ro_10yr_max_huc4','ro_20yr_max_huc4','ro_50yr_max_huc4',...
                'ro5_02yr_max_huc4','ro5_05yr_max_huc4','ro5_10yr_max_huc4','ro5_20yr_max_huc4','ro5_50yr_max_huc4');
            
            save(huc8_mat,'ID_huc8','elev_huc8','ro_mean_huc8','ro_var_huc8','ro_10yr_min_huc8','ro7_10yr_min_huc8',...
                'ro_02yr_max_huc8','ro_05yr_max_huc8','ro_10yr_max_huc8','ro_20yr_max_huc8','ro_50yr_max_huc8',...
                'ro5_02yr_max_huc8','ro5_05yr_max_huc8','ro5_10yr_max_huc8','ro5_20yr_max_huc8','ro5_50yr_max_huc8');
            
            fprintf(['Finished ' model{m} '_' force{f} '_' region{r} '\n'])
        end % end of model
    end % end of forcing data
end % end of region

toc

% % compute 1Q10 - RO_MIN10
% [ F, X ] = FDC( RO_MIN_YR, 3, 'ascend' );
% Xa = permute(X,[3,1,2,4]); % Make time dimension first one [365*wyr x lat x lon x force]
% Xb = reshape(Xa,size(Xa,1),ylat_num*xlon_num,length(force)); % Reshape to [365*wyr x lat*lon x force]
% %memory allocation
% RO_MIN10a = ones(1,ylat_num*xlon_num,length(force))*NaN;
% for i = 1:length(force)
%     RO_MIN10a(:,:,i) = interp1(F,Xb(:,:,i),0.1,'linear');
% end
% clear Xa Xb i
% RO_MIN10b = reshape(RO_MIN10a,size(RO_MIN10a,1),ylat_num,xlon_num,length(force));
% RO_MIN10 = squeeze( permute(RO_MIN10b,[2,3,1,4]) );
% clear RO_MIN10a RO_MIN10b
% clear X F
