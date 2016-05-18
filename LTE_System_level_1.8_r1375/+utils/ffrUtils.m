classdef ffrUtils
    % Implements functions for calculation of FFR-related mappings and the such.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2011 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
    end
    
    methods(Static)
        % According to the previously calculated parameters, assigns the UEs to either the PR or FR part. Stores the information in an object to be shared by all FFR schedulers
        function per_UE_FFR_mapping_object = assign_FFR_band_to_UEs(UEs,FR_SINR_switching_dB,networkPathlossMap)
            all_UEs_pos                 = reshape([UEs.pos],2,[])';
            all_UEs_pos_pix             = LTE_common_pos_to_pixel( all_UEs_pos,networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
            all_UEs_pos_pix_lin_idx     = sub2ind(size(networkPathlossMap.SINR), all_UEs_pos_pix(:,2), all_UEs_pos_pix(:,1));
            SINRs_FR_dB                 = networkPathlossMap.SINR(all_UEs_pos_pix_lin_idx);
            FFR_mapping                 = double(SINRs_FR_dB>=FR_SINR_switching_dB); % FR users
            FFR_mapping(FFR_mapping==0) = 2;                                         % PR users
            per_UE_FFR_mapping_object   = utils.ffrUEmapping(FFR_mapping);
        end
        
        % Copy and filter out a RB grid and leaves only the values specified by the indexing vector
        function a_RB_grid = filter_RB_grid(a_RB_grid,indexing_vector)
           a_RB_grid = a_RB_grid.clone();
           a_RB_grid.user_allocation            = a_RB_grid.user_allocation(indexing_vector);
           a_RB_grid.power_allocation           = a_RB_grid.power_allocation(indexing_vector);
           a_RB_grid.power_allocation_signaling = a_RB_grid.power_allocation_signaling(indexing_vector);
           a_RB_grid.n_RB                       = sum(indexing_vector);
           a_RB_grid.PMI                        = a_RB_grid.PMI(indexing_vector);
        end
       
        % Copy and separate a feedback struct according to the indexing vector
        function [FR_feedback PR_feedback] = separate_feedback_FFR(last_received_feedbacks,FR_assignment,PR_assignment,FR_UEs_idx,PR_UEs_idx)
            FR_feedback.UE_id             = last_received_feedbacks.UE_id(FR_UEs_idx);
            FR_feedback.tx_mode           = last_received_feedbacks.tx_mode(FR_UEs_idx);
            FR_feedback.nCodewords        = last_received_feedbacks.nCodewords(FR_UEs_idx);
            FR_feedback.CQI               = last_received_feedbacks.CQI(FR_assignment,:,FR_UEs_idx);
            FR_feedback.RI                = last_received_feedbacks.RI(FR_UEs_idx);
            FR_feedback.PMI               = last_received_feedbacks.PMI(FR_assignment,FR_UEs_idx);
            FR_feedback.feedback_received = last_received_feedbacks.feedback_received(FR_UEs_idx);
            
            PR_feedback.UE_id             = last_received_feedbacks.UE_id(PR_UEs_idx);
            PR_feedback.tx_mode           = last_received_feedbacks.tx_mode(PR_UEs_idx);
            PR_feedback.nCodewords        = last_received_feedbacks.nCodewords(PR_UEs_idx);
            PR_feedback.CQI               = last_received_feedbacks.CQI(PR_assignment,:,PR_UEs_idx);
            PR_feedback.RI                = last_received_feedbacks.RI(PR_UEs_idx);
            PR_feedback.PMI               = last_received_feedbacks.PMI(PR_assignment,PR_UEs_idx);
            PR_feedback.feedback_received = last_received_feedbacks.feedback_received(PR_UEs_idx);
        end
        
        % Merges the RB_grid objects from the FR and PR parts into the parent RB grid
        function merge_RB_grids(parent_RB_grid,FR_RB_grid,PR_RB_grid,FR_RB_mapping,PR_RB_mapping)
            % clear UE allocation
            parent_RB_grid.user_allocation(:) = 0;
            
            % Process FR part
            if ~isempty(FR_RB_grid)
                parent_RB_grid.user_allocation(FR_RB_mapping)            = FR_RB_grid.user_allocation;
                parent_RB_grid.power_allocation(FR_RB_mapping)           = FR_RB_grid.power_allocation;
                parent_RB_grid.power_allocation_signaling(FR_RB_mapping) = FR_RB_grid.power_allocation_signaling;
                parent_RB_grid.PMI(FR_RB_mapping)                        = FR_RB_grid.PMI;
            end
            
            % Process PR part
            if ~isempty(PR_RB_grid)
                parent_RB_grid.user_allocation(PR_RB_mapping)            = PR_RB_grid.user_allocation;
                parent_RB_grid.power_allocation(PR_RB_mapping)           = PR_RB_grid.power_allocation;
                parent_RB_grid.power_allocation_signaling(PR_RB_mapping) = PR_RB_grid.power_allocation_signaling;
                parent_RB_grid.PMI(PR_RB_mapping)                        = PR_RB_grid.PMI;
            end
            
            % Direct sum of the TX bits
            if ~isempty(FR_RB_grid)&&~isempty(PR_RB_grid)&&(length(FR_RB_grid.size_bits)~=length(PR_RB_grid.size_bits))
                parent_RB_grid.size_bits = zeros(1,2);
                parent_RB_grid.size_bits(1:length(FR_RB_grid.size_bits)) = FR_RB_grid.size_bits;
                parent_RB_grid.size_bits(1:length(PR_RB_grid.size_bits)) = parent_RB_grid.size_bits(1:length(PR_RB_grid.size_bits)) + PR_RB_grid.size_bits;
            else
                if isempty(FR_RB_grid)
                    FR_bits = 0;
                else
                    FR_bits = FR_RB_grid.size_bits;
                end
                if isempty(PR_RB_grid)
                    PR_bits = 0;
                else
                    PR_bits = PR_RB_grid.size_bits;
                end
                parent_RB_grid.size_bits = FR_bits + PR_bits;
            end
            
            % Null out the rest of the data allocation: note that the signaling power is not changed
            unused_RBs                                            = ~(FR_RB_mapping|PR_RB_mapping);
            parent_RB_grid.user_allocation(unused_RBs)            = 0;
            parent_RB_grid.power_allocation(unused_RBs)           = 0;
            parent_RB_grid.PMI(unused_RBs)                        = NaN;
        end
        
        % Assign the FFR frequencies (reuse-3) to a hexagonal grid
        function freq_assignments = assign_frequencies_to_hex_grid(target_cell_id,eNodeBs_sectors,cell_centers)
            freq_assignments      = zeros(1,size(cell_centers,1));
            plot_freq_assignments = false;
            plot_final            = false;
            
            if plot_freq_assignments
                utils.ffrUtils.plot_freq_assignment(eNodeBs_sectors,target_cell_id,cell_centers,freq_assignments);
            end
            
            % Set the target cell to use f1 and first process process the first ring
            freq_assignments(target_cell_id) = 1;
            six_nearest      = utils.ffrUtils.get_six_closest_cells(target_cell_id,cell_centers);
            freq_assignments = utils.ffrUtils.assign_frequencies(target_cell_id,six_nearest,freq_assignments);
            
            % Get the rest of the cells
            all_neighbors_by_distance = utils.ffrUtils.get_cells_ordered_by_distance(target_cell_id,cell_centers);
            all_neighbors_by_distance = all_neighbors_by_distance(8:end);
            
            if plot_freq_assignments
                utils.ffrUtils.plot_freq_assignment(eNodeBs_sectors,target_cell_id,cell_centers,freq_assignments);
            end
            
            for neighbor_idx=1:length(all_neighbors_by_distance)
                target_cell_id_tmp                = all_neighbors_by_distance(neighbor_idx);
                six_nearest                       = utils.ffrUtils.get_six_closest_cells(target_cell_id_tmp,cell_centers);
                [freq_assignments done_something] = utils.ffrUtils.assign_frequencies(target_cell_id_tmp,six_nearest,freq_assignments);
                if done_something&&plot_freq_assignments
                    utils.ffrUtils.plot_freq_assignment(eNodeBs_sectors,target_cell_id_tmp,cell_centers,freq_assignments);
                end
                if sum(freq_assignments~=0)==length(freq_assignments)
                    break;
                end
            end
            
            if plot_final
                utils.ffrUtils.plot_freq_assignment(eNodeBs_sectors,target_cell_id,cell_centers,freq_assignments);
            end
        end
        
        % Get the six (by distance) closest cells and return them ordered by angle
        function six_nearest_sorted_per_angle = get_six_closest_cells(target_cell_id,cell_centers)
            N_cells          = size(cell_centers,1);
            cell_pos        = cell_centers(target_cell_id,:);
            % Get the closest cells
            [distances,index]       = sortrows(sqrt(sum((cell_centers-cell_pos(ones(N_cells,1),:)).^2,2)));
            
            % safeguard for the border cases (just 5 neighbors)
            if distances(7)/mean(distances(2:6))>1.3
                six_nearest = index(2:6);
            else
                six_nearest = index(2:7);
            end
            six = length(six_nearest);

            six_nearest_pos = cell_centers(six_nearest,:);
             % Order the cells in angular distance
            [~,index]       = sortrows(angle(sum((six_nearest_pos-cell_pos(ones(six,1),:)).*[ones(six,1) 1i*ones(six,1)],2)));
            six_nearest_sorted_per_angle = six_nearest(index);
        end
        
        % Get all of the cells ordered by distance
        function [ordered_cells ordered_distances] = get_cells_ordered_by_distance(target_cell_id,cell_centers)
             N_cells                          = size(cell_centers,1);
            cell_pos                          = cell_centers(target_cell_id,:);
            [ordered_distances,ordered_cells] = sortrows(sqrt(sum((cell_centers-cell_pos(ones(N_cells,1),:)).^2,2)));
        end
        
        % Assign an appropriate frequency to the six neighboring cells
        function [freq_assignments did_something] = assign_frequencies(target_cell_id,six_nearest_sorted_per_angle,freq_assignments)
            total_freqs  = 1:3;
            
            freq_assignments_six = freq_assignments(six_nearest_sorted_per_angle);
            total_neighbors = length(six_nearest_sorted_per_angle);
            
            % check if there is actually anything to do
            if ~(sum(freq_assignments_six~=0)==total_neighbors)
                did_something = true;
                
                % Choose what frequencies to use
                if sum(freq_assignments_six==0)==total_neighbors %freq_assignments(target_cell_id)
                    % Probably the first cell (no neighbour frequency is assigned yet)
                    freqs_to_use = total_freqs(total_freqs~=freq_assignments(target_cell_id));
                else
                    first_non_zero       = find(freq_assignments_six,1,'first'); % find the first non-null neighbor frequency assignment
                    
                    if ~isempty(first_non_zero)
                        six_nearest_sorted_per_angle = circshift(six_nearest_sorted_per_angle,-(first_non_zero-1)); % make the first non-null the first to freq-allocate
                    end
                    
                    freqs_to_use = nonzeros(unique(freq_assignments_six));
                    if isempty(find(freqs_to_use==1,1))
                        freq_assignments(target_cell_id) = 1;
                    elseif isempty(find(freqs_to_use==2,1))
                        freq_assignments(target_cell_id) = 2;
                    else
                        freq_assignments(target_cell_id) = 3;
                    end
                    
                    if ~isempty(first_non_zero)
                        first_non_zero_freq = freq_assignments(six_nearest_sorted_per_angle(1));
                        freqs_to_use        = circshift(freqs_to_use,-(find(freqs_to_use==first_non_zero_freq,1)-1)); % begin with the frequency of the first non-null
                    end
                end
                
                % Assign frequencies
                for cell_idx=1:length(six_nearest_sorted_per_angle)
                    % To avoid overwriting
                    if freq_assignments(six_nearest_sorted_per_angle(cell_idx))==0
                        freq_assignments(six_nearest_sorted_per_angle(cell_idx)) = freqs_to_use(mod(cell_idx-1,length(freqs_to_use))+1);
                    end
                end
            else
                did_something = false;
            end
        end
        
        function plot_freq_assignment(eNodeBs_sectors,center,cell_centers,freq_assignments,varargin)
            if ~isempty(varargin)
                SizeData = varargin{1};
            else
                SizeData = [];
            end
            figure;
            hold all
            if ~isempty(center)
                scatter(cell_centers(center,1),cell_centers(center,2),'ok');
            end
            for cell_=1:length(eNodeBs_sectors)
                switch freq_assignments(cell_)
                    case 1
                        color = 'r';
                    case 2
                        color = 'g';
                    case 3
                        color = 'b';
                    otherwise
                        color = 0.25*[1 1 1];
                end
                if isempty(SizeData)
                    scatter(cell_centers(cell_,1),cell_centers(cell_,2),'.','MarkerEdgeColor', color);
                else
                    scatter(cell_centers(cell_,1),cell_centers(cell_,2),'.','MarkerEdgeColor', color,'SizeData',1000);
                end
            end
            % grid on
            axis equal
        end
        
        function plot_FFR_sim_results_over_BFR(sim_results)
            % Get curves over BFR
            reduction_factor = 1000;
            BFR_vect  = sim_results.BFR_vect;
            FFR_means = [sim_results.FFR_capacity_density_ECDF.mean_x];
            FFR_p05   = [sim_results.FFR_capacity_density_ECDF.p05];
            FFR_p95   = [sim_results.FFR_capacity_density_ECDF.p95];
            
            R1_mean   = sim_results.reuse1_capacity_density_ECDF.mean_x(ones(1,length(BFR_vect)));
            R1_p05    = sim_results.reuse1_capacity_density_ECDF.p05(ones(1,length(BFR_vect)));
            R1_p95    = sim_results.reuse1_capacity_density_ECDF.p95(ones(1,length(BFR_vect)));
            
            R3_mean   = sim_results.reuse3_capacity_density_ECDF.mean_x(ones(1,length(BFR_vect)));
            R3_p05    = sim_results.reuse3_capacity_density_ECDF.p05(ones(1,length(BFR_vect)));
            R3_p95    = sim_results.reuse3_capacity_density_ECDF.p95(ones(1,length(BFR_vect)));
            
            p95_diffs  = sim_results.maximums.p95_diffs_percent;
            mean_diffs = sim_results.maximums.mean_diffs_percent;
            p05_diffs  = sim_results.maximums.p05_diffs_percent;
            
            line_width = 1.5;
            
            figure('units','normalized','outerposition',[0.5 0.5 0.5 0.5]);
            hold all;
            plot(BFR_vect,FFR_p95/reduction_factor,  'b-','DisplayName','FFR peak (95%) capacity density','LineWidth',line_width);
            plot(BFR_vect,FFR_means/reduction_factor,'r-','DisplayName','FFR mean capacity density','LineWidth',line_width);
            plot(BFR_vect,FFR_p05/reduction_factor,  'g-','DisplayName','FFR edge (5%) capacity density','LineWidth',line_width);
            
            plot(BFR_vect,R1_p95/reduction_factor,  'b--','DisplayName','Reuse-1 peak (95%) capacity density','LineWidth',line_width);
            plot(BFR_vect,R1_mean/reduction_factor, 'r--','DisplayName','Reuse-1 mean capacity density','LineWidth',line_width);
            plot(BFR_vect,R1_p05/reduction_factor,  'g--','DisplayName','Reuse-1 edge (5%) capacity density','LineWidth',line_width);
            
            the_ylim = ylim;
            ylim([0 the_ylim(2)]);
            
            % plot(BFR_vect,R3_p95/reduction_factor,  'b:','DisplayName','Reuse-3 peak (95%) capacity density','LineWidth',line_width);
            % plot(BFR_vect,R3_mean/reduction_factor, 'r:','DisplayName','Reuse-3 mean capacity density','LineWidth',line_width);
            % plot(BFR_vect,R3_p05/reduction_factor,  'g:','DisplayName','Reuse-3 edge (5%) capacity density','LineWidth',line_width);
            
            legend('show','Location','SouthEastOutside');
            
            offset = 0.02;
            
            text(sim_results.maximums.p95_BFR+offset , sim_results.maximums.p95 /reduction_factor,  sprintf('  %3.2f%% of R_{1-95%%}\n  %3.2f%% of R_{1-mean}\n  %3.2f%% of R_{1-5%%}',p95_diffs(1),p95_diffs(2),p95_diffs(3)),'HorizontalAlignment','left','VerticalAlignment','middle','BackgroundColor',[0.75 0.75 0.75]);
            text(sim_results.maximums.mean_BFR+offset, sim_results.maximums.mean/reduction_factor,sprintf('  %3.2f%% of R_{1-95%%}\n  %3.2f%% of R_{1-mean}\n  %3.2f%% of R_{1-5%%}',mean_diffs(1),mean_diffs(2),mean_diffs(3)),'HorizontalAlignment','left','VerticalAlignment','middle','BackgroundColor',[0.75 0.75 0.75]);
            text(sim_results.maximums.p05_BFR+offset , sim_results.maximums.p05 /reduction_factor,  sprintf('  %3.2f%% of R_{1-95%%}\n  %3.2f%% of R_{1-mean}\n  %3.2f%% of R_{1-5%%}',p05_diffs(1),p05_diffs(2),p05_diffs(3)),'HorizontalAlignment','left','VerticalAlignment','middle','BackgroundColor',[0.75 0.75 0.75]);
            scatter(sim_results.maximums.p95_BFR     , sim_results.maximums.p95 /reduction_factor,'.k','SizeData',500);
            scatter(sim_results.maximums.mean_BFR    , sim_results.maximums.mean/reduction_factor,'.k','SizeData',500);
            scatter(sim_results.maximums.p05_BFR     , sim_results.maximums.p05 /reduction_factor,'.k','SizeData',500);
            
            grid on;
            xlabel('B_{FR}');
            switch sim_results.FFR_mode
                case 'shannon'
                    title(sprintf('Capacity density over B_{FR}: Blog_2(1+SINR)'));
                    ylabel('Capacity density [kbit/s/m^2]');
                otherwise
                    title(sprintf('Throughput density over B_{FR}: %s',sim_results.FFR_mode));
                    ylabel('Throughput density [kbit/s/m^2]');
            end
        end
        
        function plot_capacity_density_ECDF(sim_results)
            reduction_factor = 1000;
            BFR_vect         = sim_results.BFR_vect;
            cdf_figure       = figure;
            cdf_axes         = axes('Parent',cdf_figure,'LineWidth',1,'FontSize',12);
            color_scheme     = ['r' 'g' 'b'];
            hold(cdf_axes,'all');
            plot(cdf_axes,sim_results.reuse1_capacity_density_ECDF.x/reduction_factor, sim_results.reuse1_capacity_density_ECDF.f, color_scheme(1), 'LineWidth', 1,'DisplayName','Reuse-1');
            plot(cdf_axes,sim_results.reuse3_capacity_density_ECDF.x/reduction_factor, sim_results.reuse3_capacity_density_ECDF.f, color_scheme(2), 'LineWidth', 1,'DisplayName','Reuse-3');
            for BFR_idx = 1:length(BFR_vect)
                plot(cdf_axes,sim_results.FFR_capacity_density_ECDF(BFR_idx).x/reduction_factor, sim_results.FFR_capacity_density_ECDF(BFR_idx).f, color_scheme(3), 'LineWidth', 1,'DisplayName','FFR');
                if BFR_idx==1
                    legend(cdf_axes, 'Location', 'SouthEast');
                end
                scatter(cdf_axes,sim_results.FFR_capacity_density_ECDF(BFR_idx).mean_x/reduction_factor, sim_results.FFR_capacity_density_ECDF(BFR_idx).mean_f,  ['.' color_scheme(3)]);
                % text(sim_results.FFR_capacity_density_ECDF(BFR_idx).mean(1), sim_results.FFR_capacity_density_ECDF(BFR_idx).mean(2),sprintf('%g',BFR_vect(BFR_idx)),'HorizontalAlignment','center','VerticalAlignment','middle');
            end
            % To avoid overlapping
            plot(cdf_axes,sim_results.reuse1_capacity_density_ECDF.x/reduction_factor, sim_results.reuse1_capacity_density_ECDF.f, color_scheme(1), 'LineWidth', 1,'DisplayName','Reuse-1');
            plot(cdf_axes,sim_results.reuse3_capacity_density_ECDF.x/reduction_factor, sim_results.reuse3_capacity_density_ECDF.f, color_scheme(2), 'LineWidth', 1,'DisplayName','Reuse-3');
            scatter(cdf_axes,sim_results.reuse1_capacity_density_ECDF.mean_x/reduction_factor, sim_results.reuse1_capacity_density_ECDF.mean_f,  ['.' color_scheme(1)]);
            scatter(cdf_axes,sim_results.reuse3_capacity_density_ECDF.mean_x/reduction_factor, sim_results.reuse3_capacity_density_ECDF.mean_f,  ['.' color_scheme(2)]);
            xlabel(cdf_axes,'capacity density c [kbit/s/m^2]');
            ylabel(cdf_axes,'F(x)');
            if length(unique(BFR_vect))==1
                BRF_string = sprintf('%g',BFR_vect(1));
            else
                BRF_string = sprintf('%g-%g',min(BFR_vect),max(BFR_vect));
            end
            title(cdf_axes,sprintf('\\beta_{FR}=%s',BRF_string));
            grid(cdf_axes,'on');
            hold(cdf_axes,'off');
        end
        
        function plot_capacity_maps(sim_results,BFR_idx)
            if sim_results.compact_results
                fprintf('Maps cannot be plotted from compact results. Please generate the files again with sotrage mode set to "full"');
            end
            
            BFR                    = sim_results.BFR_vect(BFR_idx);
            target_position_matrix = sim_results.target_position_matrix;
            FR_area                = sim_results.FR_area(:,:,BFR_idx);
            PR_area                = sim_results.PR_area(:,:,BFR_idx);
            if BFR==1
                the_edges = edge(FR_area);
            else
                the_edges = edge(PR_area);
            end
            target_position_matrix_edges = edge(target_position_matrix);
            
            SINR_caxis             = [-10 25];
            capacity_caxis         = [min([sim_results.FFR_capacity_density_ECDF.min]) max([sim_results.FFR_capacity_density_ECDF.max])]/1000;
            area_caxis             = [0 2];
            nRow_plots = 3;
            nCol_plots = 3;
            
            SINR_FR_log = sim_results.SINR_FR_log;
            SINR_PR_log = sim_results.SINR_PR_log;
            
            [hor ver] = find(target_position_matrix);
            hor_roi = round([min(hor) max(hor)]);
            ver_roi = round([min(ver) max(ver)]);
            
            FR_SINR_image      = SINR_FR_log;
            PR_SINR_image      = SINR_PR_log;
            FFR_SINR_image     = SINR_FR_log.*FR_area + SINR_PR_log.*PR_area;
            FR_capacity_image  = sim_results.capacity_density_FR(:,:,BFR_idx);
            PR_capacity_image  = sim_results.capacity_density_PR(:,:,BFR_idx);
            FFR_capacity_image = FR_capacity_image.*FR_area + PR_capacity_image.*PR_area;
            [PR_area_edge_hor PR_area_edge_ver] = find(the_edges(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2)));
            [area_edge_hor    area_edge_ver]    = find(target_position_matrix_edges(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2)));
            
            FR_SINR_image_plot      = FR_SINR_image(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            PR_SINR_image_plot      = PR_SINR_image(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            FFR_SINR_image_plot     = FFR_SINR_image(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            FR_capacity_image_plot  = FR_capacity_image(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            PR_capacity_image_plot  = PR_capacity_image(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            FFR_capacity_image_plot = FFR_capacity_image(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            FR_area_image_plot      = FR_area(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            PR_area_image_plot      = area_caxis(2)*PR_area(hor_roi(1):hor_roi(2),ver_roi(1):ver_roi(2));
            FFR_area_image_plot     = FR_area_image_plot + area_caxis(2)*PR_area_image_plot;
            
            figure('Name', sprintf('BFR=%g',BFR),'units','normalized','outerposition',[0 0 1 1]);
            subplot(nRow_plots,nCol_plots,1);
            imagesc(FR_SINR_image_plot);
            colorbar;
            caxis(SINR_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(PR_area_edge_ver,PR_area_edge_hor,'.k','SizeData',0.5);
            title('SINR_{FR} [dB]');
            
            subplot(nRow_plots,nCol_plots,2);
            imagesc(PR_SINR_image_plot);
            colorbar;
            caxis(SINR_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(PR_area_edge_ver,PR_area_edge_hor,'.k','SizeData',0.5);
            title('SINR_{PR} [dB]');
            
            subplot(nRow_plots,nCol_plots,3);
            imagesc(FFR_SINR_image_plot);
            colorbar;
            caxis(SINR_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(PR_area_edge_ver,PR_area_edge_hor,'.k','SizeData',0.5);
            title('SINR_{FFR} [dB]');
            
            subplot(nRow_plots,nCol_plots,4);
            imagesc(FR_capacity_image_plot/1000); colorbar;
            caxis(capacity_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(PR_area_edge_ver,PR_area_edge_hor,'.k','SizeData',0.5);
            title('Capacity-density_{FR} [kbit/s/m^2]');
            
            subplot(nRow_plots,nCol_plots,5);
            imagesc(PR_capacity_image_plot/1000); colorbar;
            caxis(capacity_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(PR_area_edge_ver,PR_area_edge_hor,'.k','SizeData',0.5);
            title('Capacity-density_{PR} [kbit/s/m^2]');
            
            subplot(nRow_plots,nCol_plots,6);
            imagesc(FFR_capacity_image_plot/1000); colorbar;
            caxis(capacity_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(PR_area_edge_ver,PR_area_edge_hor,'.k','SizeData',0.5);
            title('Capacity-density_{FFR} [kbit/s/m^2]');
            
            subplot(nRow_plots,nCol_plots,7);
            imagesc(uint8(FR_area_image_plot)); colorbar;
            caxis(area_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(area_edge_ver,area_edge_hor,'.k','SizeData',0.5);
            title('FR area (target sector)');
            
            subplot(nRow_plots,nCol_plots,8);
            imagesc(uint8(PR_area_image_plot));  colorbar;
            caxis(area_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(area_edge_ver,area_edge_hor,'.k','SizeData',0.5);
            title('PR area (target sector)');
            
            subplot(nRow_plots,nCol_plots,9);
            imagesc(uint8(FFR_area_image_plot));  colorbar;
            caxis(area_caxis);
            set(gca,'YDir','normal');
            hold all
            scatter(area_edge_ver,area_edge_hor,'.k','SizeData',0.5);
            title('FFR area (target sector)');
        end
        
        function sim_results = calculate_optimum_BFR(LTE_config,eNodeBs,eNodeBs_sectors,networkPathlossMap,BFR_stepsize,storage_mode,FFR_mode,varargin)
            % Martin Taranetz, 2010. Modified and extended by Josep COlom Ikuno, 2011
            % (c) 2011 by INTHFT
            % www.nt.tuwien.ac.at
            
            %% Initial Parameters
            full_bandwidth          = LTE_config.bandwidth;                              %  [Hz]
            noise_power_density     = 10^(0.1*LTE_config.UE.thermal_noise_density)/1000; %  [W/Hz]: e.g. 3.9811e-021
            pathlossmaps_lin        = 10.^(-networkPathlossMap.pathloss/10);
            frequency_assignment    = LTE_config.scheduler_params.frequency_assignment;
            N_sectors               = length(eNodeBs_sectors);
            eff_bandwidth           = LTE_config.N_RB*12*LTE_config.N_sym*2/LTE_config.TTI_length;
            tx_bandwidth            = LTE_config.N_RB*LTE_config.RB_bandwidth;
            
            switch storage_mode
                case 'full'
                    compact_results = false;
                case 'compact'
                    compact_results = true;
                otherwise
                    error('"mode" is either "full" or "compact"');
            end
            sim_results.compact_results = compact_results;
            
            switch FFR_mode
                case 'shannon'
                    capacity = @(SINR_lin)log2(1+SINR_lin);
                otherwise
                    throughput_mapper                                  = varargin{1};
                    [LTE_config.tx_mode LTE_config.nTX LTE_config.nRX] = utils.miscUtils.string_long_to_tx_mode(FFR_mode);
                    capacity = @(SINR_lin)throughput_mapper.SNR_to_spectral_efficiency(10*log10(SINR_lin),LTE_config.tx_mode,LTE_config.nTX,LTE_config.nRX);
            end
            
            %% BFR_values and Memory
            BFR_vect =  0:BFR_stepsize:1;% 1;[0:0.05:1]; % 0.754
            BPR_vect = (1-BFR_vect)/3;
            
            % Use assignment for FRF = 3
            % Target eNodeB and Sector: considering only center eNodeB
            b_ = LTE_config.target_sector(1);
            s_ = LTE_config.target_sector(2);
            c_ = LTE_config.target_cell;
            
            % Pixel position of target eNodeB
            roi_x            = networkPathlossMap.roi_x;
            roi_y            = networkPathlossMap.roi_y;
            sim_radius       = LTE_config.inter_eNodeB_distance/2;
            sim_frequency    = LTE_config.frequency;
            map_res          = networkPathlossMap.data_res;
            target_pos       = eNodeBs(b_).pos;                     % Target eNodeB Position [m]
            target_pos_pixel = LTE_common_pos_to_pixel( [roi_x(2) roi_y(2)], [target_pos(1) target_pos(2)], map_res);
            
            tx_powers       = [eNodeBs_sectors.max_power];
            numel_tx_powers = length(unique(tx_powers));
            
            if numel_tx_powers==1
                transmit_power_density = tx_powers(1)/tx_bandwidth; % [W/Hz]
                power_density_all_lin  = transmit_power_density*pathlossmaps_lin;
            end
            
            %% Calculate power of signal/(noise+interferers)
            cells_to_calculate = 1:length(eNodeBs_sectors);
            SINR_FR_lin = zeros([size(pathlossmaps_lin,1) size(pathlossmaps_lin,2) length(cells_to_calculate)]);
            SINR_PR_lin = zeros([size(pathlossmaps_lin,1) size(pathlossmaps_lin,2) length(cells_to_calculate)]);
            for current_cell=cells_to_calculate
                target_cell_freq = LTE_config.scheduler_params.frequency_assignment(current_cell);
                interferers_FR   = [1:(current_cell-1) (current_cell+1):N_sectors]; % Take out yourself
                interferers_PR   = interferers_FR(frequency_assignment(interferers_FR)==target_cell_freq);
                
                if numel_tx_powers~=1
                    transmit_power_density = tx_powers(current_cell)/tx_bandwidth; % [W/Hz]
                    power_density_all_lin  = transmit_power_density*pathlossmaps_lin;
                end
                
                signal_power_density            = power_density_all_lin(:,:,current_cell);
                interferer_power_density_FR_lin = power_density_all_lin(:,:,interferers_FR);
                interferer_power_density_PR_lin = power_density_all_lin(:,:,interferers_PR);
                SINR_FR_lin(:,:,current_cell)   = signal_power_density./(noise_power_density + sum(interferer_power_density_FR_lin,3));
                SINR_PR_lin(:,:,current_cell)   = signal_power_density./(noise_power_density + sum(interferer_power_density_PR_lin,3));
            end
            SINR_FR_lin = SINR_FR_lin(:,:,cells_to_calculate);
            SINR_PR_lin = SINR_PR_lin(:,:,cells_to_calculate);
            SINR_FR_dB  = 10.*log10(SINR_FR_lin);
            SINR_PR_dB  = 10.*log10(SINR_PR_lin);
            
            %% SINR map
            [sector_SINR_FR sector_assignment_FR] = max(SINR_FR_dB,[],3);
            [sector_SINR_PR sector_assignment_PR] = max(SINR_PR_dB,[],3);
            
            %% Determine positions of target cell
            target_position_matrix      = (networkPathlossMap.sector_assignment == c_);
            target_cell_area            = sum(target_position_matrix(:));
            target_cell.SINR_FR_lin     = SINR_FR_lin(:,:,c_);
            target_cell.SINR_PR_lin     = SINR_PR_lin(:,:,c_);
            target_cell.SINR_FR_dB      = 10*log10(target_cell.SINR_FR_lin);
            target_cell.SINR_PR_dB      = 10*log10(target_cell.SINR_PR_lin);
            target_cell.SINR_FR_max_lin = max(target_cell.SINR_FR_lin(target_position_matrix));
            target_cell.SINR_FR_min_lin = min(target_cell.SINR_FR_lin(target_position_matrix));
            
            target_cell.capacity_only_FR  = eff_bandwidth / target_cell_area    .* capacity(target_cell.SINR_FR_lin); % The reuse 1 case
            target_cell.capacity_only_PR  = eff_bandwidth / (3*target_cell_area).* capacity(target_cell.SINR_PR_lin); % The reuse 3 case
            
            %% Preallocate matrix for simulation results
            sim_results.BFR_vect               = BFR_vect;
            nBFR                               = length(BFR_vect);
            pathloss_map_size_pixels           = size(networkPathlossMap.pathloss);
            sim_results.FFR_mode               = FFR_mode;
            
            if ~compact_results
                sim_results.target_position_matrix = target_position_matrix;
                sim_results.capacity_density_FFR   = zeros([pathloss_map_size_pixels(1),pathloss_map_size_pixels(2),nBFR]);
                sim_results.capacity_density_FR    = zeros([pathloss_map_size_pixels(1),pathloss_map_size_pixels(2),nBFR]);
                sim_results.capacity_density_PR    = zeros([pathloss_map_size_pixels(1),pathloss_map_size_pixels(2),nBFR]);
                sim_results.FR_area                = false([pathloss_map_size_pixels(1),pathloss_map_size_pixels(2),nBFR]);
                sim_results.PR_area                = false([pathloss_map_size_pixels(1),pathloss_map_size_pixels(2),nBFR]);
                sim_results.SINR_FR_log            = SINR_FR_dB(:,:,c_);
                sim_results.SINR_PR_log            = SINR_PR_dB(:,:,c_);
            end
            
            %% Capacity, Mean and CDFs using Reuse 1 and Reuse 3 for comparison
            % Restrict to target cell and cell edge
            target_cell_capacity_FR_vect = target_cell.capacity_only_FR(target_position_matrix);
            target_cell_capacity_PR_vect = target_cell.capacity_only_PR(target_position_matrix);
            
            % Calculate ECDF-related values
            sim_results.reuse1_capacity_density_ECDF = utils.miscUtils.ecdf(target_cell_capacity_FR_vect(:));
            sim_results.reuse3_capacity_density_ECDF = utils.miscUtils.ecdf(target_cell_capacity_PR_vect(:));
            
            %%  Entering simulation loop :
            
            % For calculating the equivalent area
            number_of_SINR_steps = 200;
            max_SINR_log  = 10*log10(target_cell.SINR_FR_max_lin-0.01);
            min_SINR_log  = 10*log10(target_cell.SINR_FR_min_lin+0.01);
            SINR_stepsize = (max_SINR_log - min_SINR_log)/number_of_SINR_steps;
            
            capacity_FR_scaled      = eff_bandwidth*capacity(target_cell.SINR_FR_lin);
            capacity_PR_scaled      = eff_bandwidth*capacity(target_cell.SINR_PR_lin);
            target_cell_FR_SINR_lin = target_cell.SINR_FR_lin.*target_position_matrix;
            
            SINR_switching_vect_log = max_SINR_log:-SINR_stepsize:min_SINR_log;
            SINR_switching_vect_lin = 10.^(SINR_switching_vect_log/10);
            
            fprintf('Entering simulation loop : %s\n',FFR_mode);
            for BFR_idx = 1:length(BFR_vect)
                
                BFR        = BFR_vect(BFR_idx);
                BPR        = BPR_vect(BFR_idx);
                break_loop = false;
                
                fprintf('BFR : %g', BFR);
                
                %%  Capacity Calculation
                for SINR_idx = 1:length(SINR_switching_vect_lin)
                    % Handle the extreme cases
                    if BFR==0
                        capacity_density_FR_  = zeros(size(capacity_PR_scaled));
                        capacity_density_PR_  = BPR/sum(target_position_matrix(:)) .* capacity_PR_scaled;
                        capacity_density_FFR  = capacity_density_PR_;
                        FR_area_posterior     = false(size(target_position_matrix));
                        PR_area_posterior     = target_position_matrix;
                        SINR_switching_log    = NaN;
                        break_loop            = true;
                    elseif BFR==1
                        capacity_density_FR_  = BFR/sum(target_position_matrix(:)) .* capacity_FR_scaled;
                        capacity_density_PR_  = zeros(size(capacity_PR_scaled));
                        capacity_density_FFR  = capacity_density_FR_;
                        PR_area_posterior     = false(size(target_position_matrix));
                        FR_area_posterior     = target_position_matrix;
                        SINR_switching_log    = NaN;
                        break_loop            = true;
                    else
                        SINR_switching_lin = SINR_switching_vect_lin(SINR_idx);
                        SINR_switching_log = SINR_switching_vect_log(SINR_idx);
                        
                        FR_area_prior         = target_cell_FR_SINR_lin >= SINR_switching_lin;
                        FR_area_prior_sum     = sum(FR_area_prior(:));
                        PR_area_prior_sum     = target_cell_area - FR_area_prior_sum;
                        
                        capacity_density_FR_  = BFR/FR_area_prior_sum .* capacity_FR_scaled;
                        capacity_density_PR_  = BPR/PR_area_prior_sum .* capacity_PR_scaled;
                        
                        FR_area_posterior = (capacity_density_FR_ > capacity_density_PR_)&target_position_matrix;
                        PR_area_posterior = (~FR_area_posterior)&target_position_matrix;
                        
                        % Found the SINR threshold for which we should not use FR anymore: where the computed sizes of the PR/FR zones are equal
                        if (FR_area_prior_sum >= sum(FR_area_posterior(:)))
                            capacity_density_FFR   = (capacity_density_FR_.*FR_area_posterior + capacity_density_PR_.*(~FR_area_posterior));
                            break_loop     = true;
                            
                            %%******************************
                            %%*        DEBUG PLOT          *
                            %%******************************
                            %figure_subplot = figure('Name', [num2str(BFR_vect(BFR_index))]);
                            %subplot(2,1,1);
                            %imagesc(FR_area_prior);
                            %title('FR area prior');
                            %subplot(2,1,2);
                            %imagesc(FR_area_posterior);
                            %title('FR area posterior');
                        end
                    end
                    
                    if break_loop
                        % Fill in results
                        if ~compact_results
                            sim_results.capacity_density_FFR(:,:,BFR_idx)  = capacity_density_FFR;
                            sim_results.capacity_density_FR(:,:,BFR_idx)   = capacity_density_FR_;
                            sim_results.capacity_density_PR(:,:,BFR_idx)   = capacity_density_PR_;
                            sim_results.FR_area(:,:,BFR_idx)               = FR_area_posterior;
                            sim_results.PR_area(:,:,BFR_idx)               = PR_area_posterior;
                        end
                        
                        sim_results.FFR_capacity_density_ECDF(BFR_idx) = utils.miscUtils.ecdf(capacity_density_FFR(target_position_matrix));
                        sim_results.FR_SINR_switching(BFR_idx)         = SINR_switching_log;
                        break;
                    end
                end
                
                fprintf(' (SINR threshold: %g dB)\n', SINR_switching_log);
            end
            
            %% Post-processing of the results
            FFR_means = [sim_results.FFR_capacity_density_ECDF.mean_x];
            FFR_p05   = [sim_results.FFR_capacity_density_ECDF.p05];
            FFR_p95   = [sim_results.FFR_capacity_density_ECDF.p95];
            [max_p95  max_p95_idx ] = max(FFR_p95);
            [max_mean max_mean_idx] = max(FFR_means);
            [max_p05  max_p05_idx ] = max(FFR_p05);
            R1_mean   = sim_results.reuse1_capacity_density_ECDF.mean_x(ones(1,length(BFR_vect)));
            R1_p05    = sim_results.reuse1_capacity_density_ECDF.p05(ones(1,length(BFR_vect)));
            R1_p95    = sim_results.reuse1_capacity_density_ECDF.p95(ones(1,length(BFR_vect)));
            sim_results.maximums.p95      = max_p95;
            sim_results.maximums.p95_BFR  = BFR_vect(max_p95_idx);
            sim_results.maximums.mean     = max_mean;
            sim_results.maximums.mean_BFR = BFR_vect(max_mean_idx);
            sim_results.maximums.p05      = max_p05;
            sim_results.maximums.p05_BFR  = BFR_vect(max_p05_idx);
            sim_results.maximums.p95_diffs_percent  = ([FFR_p95(max_p95_idx)  FFR_means(max_p95_idx)  FFR_p05(max_p95_idx)]-[R1_p95(1) R1_mean(1) R1_p05(1)])./[R1_p95(1) R1_mean(1) R1_p05(1)]*100;
            sim_results.maximums.mean_diffs_percent = ([FFR_p95(max_mean_idx) FFR_means(max_mean_idx) FFR_p05(max_mean_idx)]-[R1_p95(1) R1_mean(1) R1_p05(1)])./[R1_p95(1) R1_mean(1) R1_p05(1)]*100;
            sim_results.maximums.p05_diffs_percent  = ([FFR_p95(max_p05_idx)  FFR_means(max_p05_idx)  FFR_p05(max_p05_idx)]-[R1_p95(1) R1_mean(1) R1_p05(1)])./[R1_p95(1) R1_mean(1) R1_p05(1)]*100;
            sim_results.maximums.FR_SINR_threshold  = [ sim_results.FR_SINR_switching(max_p05_idx) sim_results.FR_SINR_switching(max_mean_idx) sim_results.FR_SINR_switching(max_p95_idx) ];
            
            fprintf('\n');
            
            %% Cool colorful plots
            if LTE_config.show_network > 0
                clims_SINR = [-100 25];
                figure_handle = utils.plotUtils.plot_SINR(LTE_config.plots.FFR_FR_SINR,SINR_FR_dB,c_,networkPathlossMap,clims_SINR,'SINR (reuse-1)');
                utils.plotUtils.add_eNodeBs_positions(figure_handle,eNodeBs,networkPathlossMap);
                utils.plotUtils.add_cell_centers(figure_handle,eNodeBs_sectors,networkPathlossMap);
                figure_handle = utils.plotUtils.plot_SINR(LTE_config.plots.FFR_FR_sector_assignment,sector_assignment_FR,1,networkPathlossMap,[],'Cells (reuse-1)');
                utils.plotUtils.add_eNodeBs_positions(figure_handle,eNodeBs,networkPathlossMap);
                utils.plotUtils.add_cell_centers(figure_handle,eNodeBs_sectors,networkPathlossMap);
                figure_handle = utils.plotUtils.plot_SINR(LTE_config.plots.FFR_PR_SINR,SINR_PR_dB,c_,networkPathlossMap,clims_SINR,'SINR (reuse-3)');
                utils.plotUtils.add_eNodeBs_positions(figure_handle,eNodeBs,networkPathlossMap);
                utils.plotUtils.add_cell_centers(figure_handle,eNodeBs_sectors,networkPathlossMap);
                figure_handle = utils.plotUtils.plot_SINR(LTE_config.plots.FFR_PR_sector_assignment,sector_assignment_PR,1,networkPathlossMap,[],'Cells (reuse-3)');
                utils.plotUtils.add_eNodeBs_positions(figure_handle,eNodeBs,networkPathlossMap);
                utils.plotUtils.add_cell_centers(figure_handle,eNodeBs_sectors,networkPathlossMap);
            end
        end
        
        function [optimum_BFR optimum_FR_SINR_switching sim_results_all] = generate_FFR_optimum_BFRs_standalone(storage_mode,skip_simulations,plot_figures)
            %% Based on Martin Taranetz's 2010 Diploma Thesis:
            %  "Cell Capacity Optimization by Fractional Frequency Partitioning"
            
            %% Load Initial Parameters
            LTE_load_params_hex_grid_tilted_FFR;
            fprintf('FFR-Optimal B_FR search\n');
            
            BFR_stepsize     = 0.01;
            FFR_mode_vec     = {'shannon' '2x2CLSM' '4x2CLSM' '4x4CLSM'};
                   
            %% Network Generation
            
            if ~skip_simulations
                % Load and plot BLER curves
                [ BLER_curves CQI_mapper ] = LTE_init_load_BLER_curves(LTE_config);
                
                % Get the eNodeBs, the Macroscopic Pathloss and the Shadow Fading
                % No need to generate shadow fading maps when using network planning tool
                if strcmp(LTE_config.network_source, 'generated')
                    [eNodeBs eNodeBs_sectors networkPathlossMap networkShadowFadingMap throughput_mapper] = LTE_init_network_generation(LTE_config,CQI_mapper);
                else
                    [eNodeBs eNodeBs_sectors networkPathlossMap throughput_mapper] = LTE_init_network_generation(LTE_config,CQI_mapper);
                end
                
                %% Calculate target sector
                [LTE_config.target_sector LTE_config.target_cell] = utils.miscUtils.get_target_sector(LTE_config,eNodeBs_sectors,networkPathlossMap);
                
                %% Calculate frequency assignment
                LTE_config.scheduler_params.frequency_assignment = utils.ffrUtils.assign_frequencies_to_hex_grid(eNodeBs(LTE_config.target_sector(1)).sectors(LTE_config.target_sector(2)).eNodeB_id,eNodeBs_sectors,networkPathlossMap.sector_centers);
                
                %% Simulations
                fprintf('Starting simulation for sectorized case...\n');
                
                optimum_BFR = zeros(1,length(FFR_mode_vec));
                for FFR_mode_idx = 1:length(FFR_mode_vec)
                    %% Simulate
                    FFR_mode    = FFR_mode_vec{FFR_mode_idx};
                    sim_results = utils.ffrUtils.calculate_optimum_BFR(LTE_config,eNodeBs,eNodeBs_sectors,networkPathlossMap,BFR_stepsize,storage_mode,FFR_mode,throughput_mapper);
                    
                    optimum_BFR(FFR_mode_idx)               = sim_results.maximums.mean_BFR;
                    optimum_FR_SINR_switching(FFR_mode_idx) = sim_results.maximums.FR_SINR_threshold(2);
                    sim_results_all(FFR_mode_idx)           = sim_results;
                    
                    %% Save results
                    file = fullfile('./data_files',sprintf('FFR_results_%s_BFR_stepsize=%g.mat',FFR_mode,BFR_stepsize));
                    
                    save(file,'sim_results');
                end
                fprintf('Simulation finished\n');
                fprintf('\n');
            end
            
            %% Plot sector assignment
            fprintf('Plotting results\n');
            fprintf('\n');
            
            %%  Plot FR-Area in FR- and PR-SINR plot
            for FFR_mode_idx = 1:length(FFR_mode_vec)
                FFR_mode = FFR_mode_vec{FFR_mode_idx};
                file = fullfile('./data_files',sprintf('FFR_results_%s_BFR_stepsize=%g.mat',FFR_mode,BFR_stepsize));
                fprintf('Loading: %s\n',file);
                loaded_data = load(file);
                sim_results = loaded_data.sim_results;
                if plot_figures
                    if LTE_config.show_network > 1
                        % Create ECDF plot
                        utils.ffrUtils.plot_capacity_density_ECDF(sim_results);
                        
                        % Plot frequency assignments
                        utils.ffrUtils.plot_freq_assignment(eNodeBs_sectors,[],networkPathlossMap.sector_centers,LTE_config.scheduler_params.frequency_assignment,1000);
                        
                        % Plot results for each B_FR case
                        for BFR_idx = 1:length(sim_results.BFR_vect)
                            utils.ffrUtils.plot_capacity_maps(sim_results,BFR_idx);
                        end
                    end
                    
                    % Plot optimum FFR points
                    if LTE_config.show_network > 0
                        utils.ffrUtils.plot_FFR_sim_results_over_BFR(sim_results);
                    end
                end
                
                if skip_simulations
                    optimum_BFR(FFR_mode_idx)               = sim_results.maximums.mean_BFR;
                    optimum_FR_SINR_switching(FFR_mode_idx) = sim_results.maximums.FR_SINR_threshold(2);
                    sim_results_all(FFR_mode_idx)           = sim_results;
                end
            end
        end
        
        function [optimum_BFR optimum_FR_SINR_switching sim_results_all] = load_FFR_optimum_BFRs(FFR_mode)
            % Loads the optimum B_FR and switching points for FFR use.
            
            %% Load Initial Parameters
            BFR_stepsize     = 0.01;
            
            %%  Plot FR-Area in FR- and PR-SINR plot
            file = fullfile('./data_files',sprintf('FFR_results_%s_BFR_stepsize=%g.mat',FFR_mode,BFR_stepsize));
            loaded_data = load(file);
            sim_results = loaded_data.sim_results;
            
            optimum_BFR               = sim_results.maximums.mean_BFR;
            optimum_FR_SINR_switching = sim_results.maximums.FR_SINR_threshold(2);
            sim_results_all           = sim_results;
        end
        
        function main_generate_files_for_SL_simulator
            load_files   = true;
            plot_figures = false;
            [optimum_BFR optimum_FR_SINR_switching sim_results] = utils.ffrUtils.generate_FFR_optimum_BFRs_standalone('compact',load_files,plot_figures);
        end
        
    end
end

