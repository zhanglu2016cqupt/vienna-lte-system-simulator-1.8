classdef miscUtils
    % Implements miscellaneous functions.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2011 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
    end
    
    methods(Static)
        function Y = fillnans(X, n, search_radius, varargin)
            % (c) Josep Colom Ikuno, INTHFT, 2010
            % Substitute NaNs and values smaller than a threshold (if given) by an
            % interpolation. Initial idea by fillnans from Ian Howat
            
            size_X = size(X);
            
            if ~isempty(varargin)
                threshold = varargin{1};
                [NaNs_pos_row NaNs_pos_col] = find(isnan(X) & (X<threshold));
            else
                [NaNs_pos_row NaNs_pos_col] = find(isnan(X));
            end
            
            N_nans = length(NaNs_pos_row);
            Y = X;
            
            for k_=1:N_nans
                % row / col position of NaN
                row_pos            = NaNs_pos_row(k_);
                col_pos            = NaNs_pos_col(k_);
                % rectangular region around NaN with width 2*searchradius+1
                % borders considered with min / max operation
                rows               = max((row_pos-search_radius),1):min((row_pos+search_radius),size_X(1));
                cols               = max((col_pos-search_radius),1):min((col_pos+search_radius),size_X(2));
                neighbors          = X(rows,cols);                  % The neighbors (NaNs included). We will take out the NaNs afterwards
                % Grid of neighbours
                rows_mat           = repmat(rows',[1 length(cols)]); % row position of all of the neighbors
                cols_mat           = repmat(cols, [length(rows) 1]); % col postion of all of the neighbors
                % Distance from current NaN
                D_mat              = sqrt((row_pos-rows_mat).^2 + (col_pos-cols_mat).^2);
                % Generate boolean matrices for masking out NaN neighbors within correct
                % distance
                no_NaN_neighbors   = ~isnan(neighbors);
                correct_distance   = D_mat<=search_radius; % from rectangular to circular distance measure
                % not_yourself     = D_mat>0; % --> not needed, as you know that you are a NaN yourself
                X_i                = neighbors(no_NaN_neighbors & correct_distance);
                D_i                = D_mat(no_NaN_neighbors & correct_distance).^n;
                Y(row_pos,col_pos) = sum(X_i./D_i)/sum(1./D_i); % Substitute the NaNs in the output
            end
        end
        
        function mod_angle = wrapTo359(the_angle)
            % Equivalent to wrapTo360, but it maps to the interval [0 359] (360 is set to zero).
            mod_angle = mod(the_angle, 360);
        end
        
        function hpol = polar2(varargin)
            %POLAR  Polar coordinate plot. Edited to plot antenna gain
            %   patterns.
            %   POLAR(THETA, RHO) makes a plot using polar coordinates of
            %   the angle THETA, in radians, versus the radius RHO.
            %   POLAR(THETA,RHO,R) uses the radial limits specified by the two element
            %   vector R.
            %   POLAR(THETA,RHO,S) uses the linestyle specified in string S.
            %   See PLOT for a description of legal linestyles.
            %   POLAR(THETA,RHO,R,S) uses the linestyle specified in string S and the
            %   radial limits in R, where R is [r_min r_max r_increase].
            %
            %   POLAR(AX,...) plots into AX instead of GCA.
            %
            %   H = POLAR(...) returns a handle to the plotted object in H.
            %
            %   Example:
            %      t = 0:.01:2*pi;
            %      polar(t,sin(2*t).*cos(2*t),'--r')
            %
            %   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
            %
            %   Revised version by Daniel Armyr, 2009. Based on Mathworks original.
            
            %   Copyright 1984-2007 The MathWorks, Inc.
            %   $Revision: 5.22.4.9 $  $Date: 2007/08/27 17:06:52 $
            
            % Parse possible Axes input
            [cax,args,nargs] = axescheck(varargin{:});
            error(nargchk(1,4,nargs,'struct'));
            
            if nargs < 1 || nargs > 4
                error('MATLAB:polar:InvalidInput', 'Requires 2 to 4 data arguments.')
            elseif nargs == 2
                theta = args{1};
                rho = args{2};
                if ischar(rho)
                    line_style = rho;
                    rho = theta;
                    [mr,nr] = size(rho);
                    if mr == 1
                        theta = 1:nr;
                    else
                        th = (1:mr)';
                        theta = th(:,ones(1,nr));
                    end
                else
                    line_style = 'auto';
                end
                radial_limits = [];
            elseif nargs == 1
                theta = args{1};
                line_style = 'auto';
                rho = theta;
                [mr,nr] = size(rho);
                if mr == 1
                    theta = 1:nr;
                else
                    th = (1:mr)';
                    theta = th(:,ones(1,nr));
                end
                radial_limits = [];
            elseif nargs == 3
                if ( ischar(args{3}) )
                    [theta,rho,line_style] = deal(args{1:3});
                    radial_limits = [];
                else
                    [theta,rho,radial_limits] = deal(args{1:3});
                    line_style = 'auto';
                    if ( ~(numel(radial_limits) == 3) )
                        error ( 'R must be a 2 element vector' );
                    end
                    
                    % Data clipping
                    rho(rho<radial_limits(1)) = radial_limits(1);
                    rho(rho>radial_limits(2)) = radial_limits(2);
                end
            else %nargs == 4
                [theta,rho,radial_limits,line_style] = deal(args{1:4});
                if ( ~(numel(radial_limits) == 3) )
                    error ( 'R must be a 2 element vector' );
                end
                
                % Data clipping
                rho(rho<radial_limits(1)) = radial_limits(1);
                rho(rho>radial_limits(2)) = radial_limits(2);
            end
            
            if ischar(theta) || ischar(rho)
                error('MATLAB:polar:InvalidInputType', 'Input arguments must be numeric.');
            end
            if ~isequal(size(theta),size(rho))
                error('MATLAB:polar:InvalidInput', 'THETA and RHO must be the same size.');
            end
            
            % get hold state
            cax = newplot(cax);
            
            next = lower(get(cax,'NextPlot'));
            hold_state = ishold(cax);
            
            % get x-axis text color so grid is in same color
            tc = get(cax,'xcolor');
            ls = get(cax,'gridlinestyle');
            
            % Hold on to current Text defaults, reset them to the
            % Axes' font attributes so tick marks use them.
            fAngle  = get(cax, 'DefaultTextFontAngle');
            fName   = get(cax, 'DefaultTextFontName');
            fSize   = get(cax, 'DefaultTextFontSize');
            fWeight = get(cax, 'DefaultTextFontWeight');
            fUnits  = get(cax, 'DefaultTextUnits');
            set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
                'DefaultTextFontName',   get(cax, 'FontName'), ...
                'DefaultTextFontSize',   get(cax, 'FontSize'), ...
                'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
                'DefaultTextUnits','data')
            
            % only do grids if hold is off
            if ~hold_state
                
                % make a radial grid
                hold(cax,'on');
                set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
                
                % ensure that Inf values don't enter into the limit calculation.
                arho = abs(rho(:));
                if ( isempty(radial_limits) )
                    maxrho = max(arho(arho ~= Inf));
                    minrho = 0;
                    hhh=line([minrho minrho maxrho maxrho],[minrho maxrho maxrho minrho],'parent',cax);
                    v = [get(cax,'xlim') get(cax,'ylim')];
                    ticks = numel(get(cax,'ytick'));
                    delete(hhh);
                    % check radial limits and ticks
                    rmin = v(1); rmax = v(4); rticks = max(ticks-1,2);
                    if rticks > 5   % see if we can reduce the number
                        if rem(rticks,2) == 0
                            rticks = rticks/2;
                        elseif rem(rticks,3) == 0
                            rticks = rticks/3;
                        end
                    end
                    rinc = (rmax-rmin)/rticks;
                    
                else
                    rmax = radial_limits(2);
                    rmin = radial_limits(1);
                    %         order = (10^floor(log10(rmax-rmin)));
                    %         firstDigit = floor((rmax-rmin)/order);
                    %         if ( firstDigit <= 1 )
                    %             step = 0.2*order;
                    %         elseif ( firstDigit <= 3 )
                    %             step = 0.5*order;
                    %         elseif ( firstDigit <= 7 )
                    %             step = order;
                    %         else
                    %             step = 2*order;
                    %         end
                    %         rinc = step;
                    rinc = radial_limits(3);
                end
                
                % define a circle
                th = 0:pi/50:2*pi;
                xunit = cos(th);
                yunit = sin(th);
                % now really force points on x/y axes to lie on them exactly
                inds = 1:(length(th)-1)/4:length(th);
                xunit(inds(2:2:4)) = zeros(2,1);
                yunit(inds(1:2:5)) = zeros(3,1);
                % plot background if necessary
                if ~ischar(get(cax,'color')),
                    patch('xdata',xunit*(rmax-rmin),'ydata',yunit*(rmax-rmin), ...
                        'edgecolor',tc,'facecolor',get(cax,'color'),...
                        'handlevisibility','off','parent',cax);
                end
                
                % draw radial circles
                c82 = cos(82*pi/180);
                s82 = sin(82*pi/180);
                for i=(rmin+rinc):rinc:rmax
                    hhh = line(xunit*(i-rmin),yunit*(i-rmin),'linestyle',ls,'color',tc,'linewidth',1,...
                        'handlevisibility','off','parent',cax);
                    text((i-rmin+rinc/20)*c82,(i-rmin+rinc/20)*s82, ...
                        ['  ' num2str(i)],'verticalalignment','bottom',...
                        'handlevisibility','off','parent',cax)
                end
                set(hhh,'linestyle','-') % Make outer circle solid
                
                % plot spokes
                th = (1:6)*2*pi/12;
                cst = cos(th); snt = sin(th);
                cs = [-cst; cst];
                sn = [-snt; snt];
                line((rmax-rmin)*cs,(rmax-rmin)*sn,'linestyle',ls,'color',tc,'linewidth',1,...
                    'handlevisibility','off','parent',cax)
                
                % annotate spokes in degrees
                rt = 1.1*(rmax-rmin);
                for i = 1:length(th)
                    text(rt*cst(i),rt*snt(i),int2str(i*30),...
                        'horizontalalignment','center',...
                        'handlevisibility','off','parent',cax);
                    if i == length(th)
                        loc = int2str(0);
                    else
                        loc = int2str(180+i*30);
                    end
                    text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
                        'handlevisibility','off','parent',cax)
                end
                
                % set view to 2-D
                view(cax,2);
                % set axis limits
                axis(cax,(rmax-rmin)*[-1 1 -1.15 1.15]);
                
                setappdata( cax, 'rMin', rmin );
                
            else
                %Try to find the inner radius of the current axis.
                if ( isappdata ( cax, 'rMin' ) )
                    rmin = getappdata( cax, 'rMin' );
                else
                    rmin = 0;
                end
            end
            
            % Reset defaults.
            set(cax, 'DefaultTextFontAngle', fAngle , ...
                'DefaultTextFontName',   fName , ...
                'DefaultTextFontSize',   fSize, ...
                'DefaultTextFontWeight', fWeight, ...
                'DefaultTextUnits',fUnits );
            
            % transform data to Cartesian coordinates.
            xx = (rho - rmin).*cos(theta);
            yy = (rho - rmin).*sin(theta);
            
            % plot data on top of grid
            if strcmp(line_style,'auto')
                q = plot(xx,yy,'parent',cax);
            else
                q = plot(xx,yy,line_style,'parent',cax);
            end
            
            if nargout == 1
                hpol = q;
            end
            
            if ~hold_state
                set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
            end
            set(get(cax,'xlabel'),'visible','on')
            set(get(cax,'ylabel'),'visible','on')
            
            if ~isempty(q) && ~isdeployed
                makemcode('RegisterHandle',cax,'IgnoreHandle',q,'FunctionName','polar');
            end
        end
        
        function ecdf_data = ecdf(input_vector)
            % Filter-out NaNs and Infs
            input_vector_finite = input_vector(isfinite(input_vector));
            
            % ECDF plot plus some extras (percentiles)
            [ecdf_data.f ecdf_data.x ecdf_data.flo ecdf_data.fup] = ecdf(input_vector_finite);
            
            % To allow for inverse mappings
            [unique_x m n] = unique(ecdf_data.x);
            unique_y = ecdf_data.f(m);
            
            mean_data = mean(input_vector_finite);
            min_data  = min(input_vector_finite);
            max_data  = max(input_vector_finite);
            
            ecdf_data.min        = min_data;
            ecdf_data.max        = max_data;
            ecdf_data.mean_x     = mean_data;
            ecdf_data.mean_f     = interp1(unique_x,unique_y,mean_data);
            ecdf_data.mean_log   = mean(log(input_vector_finite));
            ecdf_data.fairness   = (sum(input_vector_finite)^2)/(length(input_vector_finite)*sum(input_vector_finite.^2));
            ecdf_data.p05        = ecdf_data.x(find(ecdf_data.f>=0.05,1,'first'));
            ecdf_data.p50        = ecdf_data.x(find(ecdf_data.f>=0.5,1,'first'));
            ecdf_data.p95        = ecdf_data.x(find(ecdf_data.f>=0.95,1,'first'));
            ecdf_data.input_data = input_vector;
            ecdf_data.what       = []; % To be filled outside of this function (this avoids runtime errors)
            ecdf_data.unit       = []; % To be filled outside of this function (this avoids runtime errors)
        end
        
        function [confidence_interval confidence_interval_relative] = confidence_interval(input_data,percentile)
            % Calculates for the input vector consisting of repeated
            % measures of the same parameter the confidence interval
            sorted_input = sort(input_data);
            mean_input   = mean(sorted_input);
            min_input    = min(sorted_input);
            max_input    = max(sorted_input);
            range        = max_input - min_input;
            N_intervals  = 100;
            range_step   = range / N_intervals;
            limits       = [
                mean_input + (1:100)*range_step
                mean_input - (1:100)*range_step
                ];
            target_num = percentile/100*length(input_data);
            N_in_set = zeros(1,N_intervals);
            for i_=1:N_intervals
                N_in_set(i_) = sum(input_data<=limits(1,i_) & input_data>=limits(2,i_));
            end
            conf_int_idx = find(N_in_set>=target_num,1,'first');
            confidence_interval          = conf_int_idx*range_step;
            confidence_interval_relative = confidence_interval / mean_input;
        end
        
        function tidy_up_memory_before_closing(UEs,eNodeBs_sectors,eNodeBs)
            % It was seen that if this is not performed PRIOR to finishing a
            % simulation, the 'clear' command takes extremely long to
            % finish. Probably due to the high number of cross-references
            % and handles in memory.
            for u_=1:length(UEs)
                if ~isstruct(UEs(u_))
                    UEs(u_).clear;
                end
            end
            
            for bs_=1:length(eNodeBs_sectors)
                eNodeBs_sectors(bs_).clear;
            end
            
            for site_=1:length(eNodeBs)
                eNodeBs(site_).clear;
            end
        end
        
        function tx_mode_string = tx_mode_to_string(tx_mode)
            switch tx_mode
                case 1
                    tx_mode_string = 'SIXO';
                case 2
                    tx_mode_string = 'TxD';
                case 3
                    tx_mode_string = 'OLSM';
                case 4
                    tx_mode_string = 'CLSM';
                otherwise
                    tx_mode_string = sprintf('%d',tx_mode);
            end
        end
        
        function tx_mode_string_long = tx_mode_to_string_long(tx_mode,nTX,nRX)
            tx_mode_string_long = sprintf('%gx%g%s',nTX,nRX,utils.miscUtils.tx_mode_to_string(tx_mode));
        end
        
        function [tx_mode nTX nRX] = string_long_to_tx_mode(tx_mode_string_long)
            if strfind(tx_mode_string_long,'SISO')
                tx_mode = 1;
                nTX     = 1;
                nRX     = 1;
                return
            else
                nTX = str2double(tx_mode_string_long(1));
                nRX = str2double(tx_mode_string_long(3));
                if strfind(tx_mode_string_long,'TxD')
                    tx_mode = 2;
                    return
                end
                if strfind(tx_mode_string_long,'OLSM')
                    tx_mode = 3;
                    return
                end
                if strfind(tx_mode_string_long,'CLSM')
                    tx_mode = 4;
                    return
                end
            end
        end
        
        % perform a fitting based on the mean value of a scatterplot
        function [bin_x bin_mean points_in_bin bin_max bin_min bin_max_idx bin_min_idx] = fit_scatterplot_data(x_data,y_data,num_bins)
            the_linspace  = linspace(min(x_data),max(x_data),num_bins+1);
            bin_starts    = the_linspace(1:(end-1));
            bin_ends      = the_linspace(2:end);
            bin_x         = (bin_starts+bin_ends)/2;
            bin_mean         = zeros(size(bin_x));
            bin_max     = zeros(size(bin_x));
            bin_min     = zeros(size(bin_x));
            bin_max_idx = zeros(size(bin_x));
            bin_min_idx = zeros(size(bin_x));
            points_in_bin = zeros(size(bin_x));
            for bin_idx=1:num_bins
                points_in_bin_idxs         = (x_data>=bin_starts(bin_idx))&(x_data<bin_ends(bin_idx));
                points_in_current_bin      = sum(points_in_bin_idxs);
                points_in_bin(bin_idx)     = points_in_current_bin;
                
                if points_in_current_bin>0
                    points_in_bin_idxs_lin     = find(points_in_bin_idxs);
                    points_to_average          = y_data(points_in_bin_idxs);
                    
                    [bin_max(bin_idx) I_max] = max(points_to_average);
                    [bin_min(bin_idx) I_min] = min(points_to_average);
                    bin_max_idx(bin_idx) = points_in_bin_idxs_lin(I_max);
                    bin_min_idx(bin_idx) = points_in_bin_idxs_lin(I_min);
                    
                    points_to_average          = points_to_average(isfinite(points_to_average));
                    bin_mean(bin_idx)             = mean(points_to_average); % Average only the finite part of the vector
                else
                    bin_mean(bin_idx)    = NaN;
                    bin_max(bin_idx)     = NaN;
                    bin_min(bin_idx)     = NaN;
                    bin_max_idx(bin_idx) = NaN;
                    bin_min_idx(bin_idx) = NaN;
                end
            end
        end
        
        function [sitesCopy eNodeBsCopy] = copySitesAndeNodeBs(sites,eNodeBs)
            sitesSectors = cell(1,length(sites));
            for siteIdx=1:length(sites)
                allSectors            = [sites(siteIdx).sectors];
                sitesSectors{siteIdx} = [allSectors.eNodeB_id];
            end
            
            eNodeBsCopy = network_elements.eNodeB_sector; % Initialization
            for enodebIdx = 1:length(eNodeBs)
                eNodeBsCopy(enodebIdx) = network_elements.eNodeB_sector;
                eNodeBsCopy(enodebIdx).id                         = eNodeBs(enodebIdx).id;
                eNodeBsCopy(enodebIdx).eNodeB_id                  = eNodeBs(enodebIdx).eNodeB_id;
                eNodeBsCopy(enodebIdx).azimuth                    = eNodeBs(enodebIdx).azimuth;
                eNodeBsCopy(enodebIdx).antenna                    = eNodeBs(enodebIdx).antenna;
                eNodeBsCopy(enodebIdx).attached_UEs               = eNodeBs(enodebIdx).attached_UEs;
                eNodeBsCopy(enodebIdx).max_power                  = eNodeBs(enodebIdx).max_power;
                eNodeBsCopy(enodebIdx).signaling_power            = eNodeBs(enodebIdx).signaling_power;
                eNodeBsCopy(enodebIdx).nTX                        = eNodeBs(enodebIdx).nTX;
                eNodeBsCopy(enodebIdx).always_on                  = eNodeBs(enodebIdx).always_on;
                eNodeBsCopy(enodebIdx).macroscopic_pathloss_model = eNodeBs(enodebIdx).macroscopic_pathloss_model;
                eNodeBsCopy(enodebIdx).unquantized_CQI_feedback   = eNodeBs(enodebIdx).unquantized_CQI_feedback;
                eNodeBsCopy(enodebIdx).transmitter                = eNodeBs(enodebIdx).transmitter;
                eNodeBsCopy(enodebIdx).frequency_band             = eNodeBs(enodebIdx).frequency_band;
                eNodeBsCopy(enodebIdx).antenna_name               = eNodeBs(enodebIdx).antenna_name;
                eNodeBsCopy(enodebIdx).antenna_type               = eNodeBs(enodebIdx).antenna_type;
                eNodeBsCopy(enodebIdx).tx_height                  = eNodeBs(enodebIdx).tx_height;
                eNodeBsCopy(enodebIdx).electrical_downtilt        = eNodeBs(enodebIdx).electrical_downtilt;
                eNodeBsCopy(enodebIdx).mechanical_downtilt        = eNodeBs(enodebIdx).mechanical_downtilt;
            end
            
            sitesCopy = network_elements.eNodeB; % Initialization
            for siteIdx = 1:length(sites)
                sitesCopy(siteIdx).id        = sites(siteIdx).id;
                sitesCopy(siteIdx).pos       = sites(siteIdx).pos;
                sitesCopy(siteIdx).name      = sites(siteIdx).name;
                sitesCopy(siteIdx).altitude  = sites(siteIdx).altitude;
                sitesCopy(siteIdx).site_name = sites(siteIdx).site_name;
                sitesCopy(siteIdx).site_type = sites(siteIdx).site_type;
                for sectorIdx = 1:length(sitesSectors{siteIdx})
                    if isempty(sitesCopy(siteIdx).sectors)
                        sitesCopy(siteIdx).sectors               = eNodeBsCopy(sitesSectors{siteIdx}(sectorIdx));
                        sitesCopy(siteIdx).sectors.parent_eNodeB = sitesCopy(siteIdx);
                    else
                        sitesCopy(siteIdx).sectors(sectorIdx)               = eNodeBsCopy(sitesSectors{siteIdx}(sectorIdx));
                        sitesCopy(siteIdx).sectors(sectorIdx).parent_eNodeB = sitesCopy(siteIdx);
                    end
                end
            end
        end
    end
end

