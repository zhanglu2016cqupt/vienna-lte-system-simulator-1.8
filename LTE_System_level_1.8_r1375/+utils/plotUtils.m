classdef plotUtils
    % Some functions to make plotting easier
    
    properties
    end
    
    methods(Static)
        function figure_handle = add_eNodeBs_positions(figure_handle,eNodeBs,networkPathlossMap,varargin)
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end
            
            if isempty(varargin)
                MarkerEdgeColor = 'k';
                MarkerFaceColor = 'w';
            else
                MarkerEdgeColor = varargin{1};
                MarkerFaceColor = varargin{2};
            end
            hold all
            set(gca,'YDir','normal');
            
            % Plot all of the eNodeBs' positions
            for b_ = 1:length(eNodeBs)
                scatter(eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceColor',MarkerFaceColor);
                text(eNodeBs(b_).pos(1)+7*networkPathlossMap.data_res,eNodeBs(b_).pos(2),num2str(b_),'Color','k');
            end
        end
        
        function figure_handle = add_eNodeBs_antenna_direction_vectors(figure_handle,eNodeBs,vector_length,LineColor)
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end
            hold all
            set(gca,'YDir','normal');
            
            % Plot a line that tells where the antennas are pointing
            for b_=1:length(eNodeBs)
                origin = eNodeBs(b_).pos;
                for s_=1:length(eNodeBs(b_).sectors)
                    angle = wrapTo360(-eNodeBs(b_).sectors(s_).azimuth+90);
                    vector = vector_length*[ cosd(angle) sind(angle) ];
                    destiny = vector + origin;
                    plot([origin(1) destiny(1)],[origin(2) destiny(2)],LineColor);
                end
            end
        end
        
        function figure_handle = add_UEs_positions(figure_handle,UEs,DotColor,add_UE_ids)
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end
            hold all
            set(gca,'YDir','normal');
            
            % Add the position of each UE
            UE_pos = reshape([UEs.pos],2,length(UEs))';
            scatter(UE_pos(:,1),UE_pos(:,2),'Marker','.','MarkerFaceColor',DotColor,'MarkerEdgeColor',DotColor);
            if add_UE_ids
                for u_=1:length(UEs)
                    text(UEs(u_).pos(1)+5*1,UEs(u_).pos(2),num2str(UEs(u_).id),'FontSize',8);
                end
            end
        end
        
        function figure_handle = add_cell_centers(figure_handle,eNodeBs_sectors,networkPathlossMap)
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end
            hold all
            set(gca,'YDir','normal');
            
            % Mark with a text the center of the cells
            for s_ = 1:length(eNodeBs_sectors)
                text(networkPathlossMap.sector_centers(s_,1),networkPathlossMap.sector_centers(s_,2),num2str(s_),'HorizontalAlignment','center','VerticalAlignment','middle','Color',0.75*[1 1 1]);
            end
        end
        
        function figure_handle = plot_SINR(figure_handle,full_SINR_map,cell_idx,networkPathlossMap,clims,the_title)
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end        
            if isempty(clims)
                imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,full_SINR_map(:,:,cell_idx));
            else
                imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,full_SINR_map(:,:,cell_idx),clims);
            end
            colorbar;
            xlabel('x-pos [m]');
            ylabel('y-pos [m]');
            title(the_title);
        end
        
        function figure_handle = plot_sector_SINR_cdfs(...
                figure_handle,...
                networkPathlossMap,networkShadowFadingMap,...
                LTE_config, eNodeBs)
            % Plots the ecdf of the target cell's (not sector) SINR
            % (c) Josep Colom Ikuno, INTHFT, 2010-2012
            % www.nt.tuwien.ac.at
            %
            % Left this function without deleting it, as it still works,
            % but since it does not make much sense if RRHs are employed,
            % it is not called anymore.
            
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end
            
            if ~isempty(networkPathlossMap.SINR2) && ~isa(networkShadowFadingMap,'channel_gain_wrappers.shadowFadingDummyMap')
                use_shadowing = true;
            else
                use_shadowing = false;
            end
            
            if isfield(LTE_config ,'compute_only_UEs_from_this_eNodeBs') && ~isempty(LTE_config.compute_only_UEs_from_this_eNodeBs)
                consider_all_map = false;
                considered_eNodeBs = LTE_config.compute_only_UEs_from_this_eNodeBs;
            else
                consider_all_map = true;
            end
            
            if use_shadowing
                if consider_all_map
                    [f1, x1] = ecdf(networkPathlossMap.SINR(:));
                    [f2, x2] = ecdf(networkPathlossMap.SINR2(:));
                    
                    figure(figure_handle);
                    hold on;
                    plot(x1,f1,'k','DisplayName','SINR CDF, macro and shadow fading');
                    plot(x2,f2,':k','DisplayName','SINR CDF, macro fading only');
                    legend('show','Location','Best');
                    grid('on');
                    xlabel('SINR (dB)');
                    ylabel('F(x)');
                    title('ROI wideband SINR CDF');
                else
                    shadowing_assignment     = false(size(networkPathlossMap.sector_assignment));
                    non_shadowing_assignment = false(size(networkPathlossMap.sector_assignment2));
                    for i_=1:length(considered_eNodeBs)
                        shadowing_assignment = shadowing_assignment|(networkPathlossMap.sector_assignment==considered_eNodeBs(i_));
                        non_shadowing_assignment = non_shadowing_assignment|(networkPathlossMap.sector_assignment2==considered_eNodeBs(i_));
                    end
                    [f1, x1] = ecdf(networkPathlossMap.SINR(shadowing_assignment));
                    [f2, x2] = ecdf(networkPathlossMap.SINR2(non_shadowing_assignment));
                    
                    % Extra map (pathloss)
                    nSites = size(networkShadowFadingMap.pathloss,3);
                    all_pathloss_dB = 10*log10(networkPathlossMap.pathloss);
                    all_shadow_dB   = zeros([size(all_pathloss_dB,1) size(all_pathloss_dB,2) nSites]);
                    for s_=1:nSites
                        all_shadow_dB(:,:,s_) = imresize(10*log10(networkShadowFadingMap.pathloss(:,:,s_)),[size(all_pathloss_dB,1),size(all_pathloss_dB,2)]);
                    end
                    all_pathloss_plus_shadow_dB = all_pathloss_dB;
                    for s_=1:length(eNodeBs)
                        all_pathloss_plus_shadow_dB(:,:,s_) = all_pathloss_plus_shadow_dB(:,:,s_) + all_shadow_dB(:,:,eNodeBs(s_).parent_eNodeB.id);
                    end
                    
                    % Cut out the unneeded parts and generate a combined pathloss map
                    combined_PL_shadowing    = zeros(size(networkPathlossMap.sector_assignment));
                    combined_PL_no_shadowing = zeros(size(networkPathlossMap.sector_assignment2));
                    for s_=1:length(eNodeBs)
                        current_PL_no_shadowing = all_pathloss_dB(:,:,s_);
                        current_PL_shadowing    = all_pathloss_plus_shadow_dB(:,:,s_);
                        current_PL_no_shadowing(networkPathlossMap.sector_assignment2~=s_) = NaN;
                        current_PL_shadowing(networkPathlossMap.sector_assignment~=s_) = NaN;
                        
                        combined_PL_shadowing(networkPathlossMap.sector_assignment==s_)     = current_PL_shadowing(networkPathlossMap.sector_assignment==s_);
                        combined_PL_no_shadowing(networkPathlossMap.sector_assignment2==s_) = current_PL_no_shadowing(networkPathlossMap.sector_assignment2==s_);
                        all_pathloss_dB(:,:,s_)             = current_PL_no_shadowing;
                        all_pathloss_plus_shadow_dB(:,:,s_) = current_PL_shadowing;
                    end
                    
                    [f3, x3] = ecdf(reshape(all_pathloss_plus_shadow_dB(:,:,considered_eNodeBs),1,[]));
                    [f4, x4] = ecdf(reshape(all_pathloss_dB(:,:,considered_eNodeBs),1,[]));
                    
                    figure(figure_handle);
                    subplot(1,2,1);
                    hold on;
                    plot(x1,f1,'k','DisplayName','Macro and shadow fading');
                    plot(x2,f2,':k','DisplayName','Macro fading only');
                    legend('show','Location','SouthEast');
                    grid('on');
                    xlabel('SINR (dB)');
                    ylabel('F(x)');
                    title('Considered cells wideband SINR CDF');
                    
                    subplot(1,2,2);
                    hold on;
                    plot(x3,f3,'k','DisplayName','With shadow fading');
                    plot(x4,f4,':k','DisplayName','Without shadow fading');
                    legend('show','Location','SouthEast');
                    grid('on');
                    xlabel('Signal loss (dB)');
                    ylabel('F(x)');
                    title('Considered cells pathloss CDF');
                    
                    % Extra plots that were used for the validation scenarios
                    extra_pathloss_plots = false;
                    
                    if extra_pathloss_plots
                        % Plot cell pathloss
                        combined_PL_shadowing_cut = combined_PL_shadowing;
                        combined_PL_shadowing_cut(~shadowing_assignment) = NaN;
                        combined_PL_no_shadowing_cut = combined_PL_no_shadowing;
                        combined_PL_no_shadowing_cut(~non_shadowing_assignment) = NaN;
                        
                        TX_power = 10*log10(LTE_config.eNodeB_tx_power)+30;
                        pLoss    = LTE_config.additional_penetration_loss;
                        RX_power_shadow    = TX_power-pLoss-reshape(all_pathloss_plus_shadow_dB(:,:,considered_eNodeBs),1,[]);
                        RX_power_no_shadow = TX_power-pLoss-reshape(all_pathloss_dB(:,:,considered_eNodeBs),1,[]);
                        RX_power_shadow    = RX_power_shadow(isfinite(RX_power_shadow));
                        RX_power_no_shadow = RX_power_no_shadow(isfinite(RX_power_no_shadow));
                        [f3 x3 flo3 fup3] = ecdf(RX_power_shadow);
                        [f4 x4 flo4 fup4] = ecdf(RX_power_no_shadow);
                        figure;
                        hold on;
                        plot(x3,f3,'k','DisplayName','With shadow fading');
                        plot(x4,f4,':k','DisplayName','Without shadow fading');
                        legend('show','Location','SouthEast');
                        grid('on');
                        xlabel('Signal loss (dB)');
                        ylabel('F(x)');
                        title('Average received power CDF');
                        
                        data_to_save.SINR_CDF_shadow.f    = f1;
                        data_to_save.SINR_CDF_shadow.x    = x1;
                        data_to_save.SINR_CDF_no_shadow.f = f2;
                        data_to_save.SINR_CDF_no_shadow.x = x2;
                        data_to_save.RSRP_CDF_shadow.f    = f3;
                        data_to_save.RSRP_CDF_shadow.x    = x3;
                        data_to_save.RSRP_CDF_no_shadow.f = f4;
                        data_to_save.RSRP_CDF_no_shadow.x = x4;
                        
                        save('CDF_data.mat','-struct','data_to_save');
                        
                        figure;
                        imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,combined_PL_shadowing,[50 145]);
                        xlabel('x pos [m]');
                        ylabel('y pos [m]');
                        title('cell pathloss with shadowing');
                        colorbar;
                        
                        figure;
                        imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,combined_PL_no_shadowing,[50 145]);
                        xlabel('x pos [m]');
                        ylabel('y pos [m]');
                        title('cell pathloss, no shadowing');
                        colorbar;
                        
                        figure;
                        imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,combined_PL_shadowing_cut,[50 145]);
                        xlabel('x pos [m]');
                        ylabel('y pos [m]');
                        title('cell pathloss with shadowing');
                        colorbar;
                        
                        figure;
                        imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,combined_PL_no_shadowing_cut,[50 145]);
                        xlabel('x pos [m]');
                        ylabel('y pos [m]');
                        title('cell pathloss, no shadowing');
                        colorbar;
                    end
                end
            else
                if consider_all_map
                    [f1, x1] = ecdf(networkPathlossMap.SINR(:));
                else
                    assignment  = false(size(networkPathlossMap.sector_assignment));
                    for i_=1:length(considered_eNodeBs)
                        assignment = assignment|(networkPathlossMap.sector_assignment==considered_eNodeBs(i_));
                    end
                    [f1, x1] = ecdf(networkPathlossMap.SINR(assignment));
                end
                
                figure(figure_handle);
                hold on;
                plot(x1,f1,'k','DisplayName','SINR CDF');
                legend('show','Location','Best');
                grid('on');
                xlabel('SINR (dB)');
                ylabel('F(x)');
                title('ROI SINR CDF');
            end
        end
        
        % Plot eNodeB and UE positions
        function figure_handle = plot_eNodeBs_and_UEs(figure_handle,eNodeBs,UEs,networkPathlossMap,add_UE_ids)
            if ~isempty(figure_handle)
                figure(figure_handle);
            else
                figure_handle = figure;
            end
            
            hold all;
            grid on;
            
            utils.plotUtils.add_UEs_positions(figure_handle,UEs,'black',add_UE_ids);
            utils.plotUtils.add_eNodeBs_antenna_direction_vectors(figure_handle,eNodeBs,40,'blue');
            utils.plotUtils.add_eNodeBs_positions(figure_handle,eNodeBs,networkPathlossMap,'black','red');
            
            xlim(networkPathlossMap.roi_x);
            ylim(networkPathlossMap.roi_y);
            
            title('eNodeB and UE positions');
            xlabel('x pos [m]');
            ylabel('y pos [m]');
            axis equal;
        end
    end
    
end

