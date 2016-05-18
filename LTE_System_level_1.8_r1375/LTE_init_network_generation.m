function [sites, eNodeBs, networkPathlossMap, networkShadowFadingMap] = LTE_init_network_generation(LTE_config,varargin)
% Network generation. Either from a file (cache) or calling the necessary function.
% (c) Josep Colom Ikuno, INTHFT, 2009
% www.nt.tuwien.ac.at

if ~LTE_config.trace_simulation_mode
    vararginIsEmpy = isempty(varargin);
    if ~vararginIsEmpy
        pathlossDataIsEmpty = isempty(varargin{1}{2});
    else
        pathlossDataIsEmpty = true;
    end
    
    if LTE_config.usePathlossMapFromMemory&&pathlossDataIsEmpty
        error('When setting LTE_config.usePathlossMapFromMemory to true, pathloss data must be input.\n');
    end
    
    if strcmp(LTE_config.network_source,'fixed pathloss')
        LTE_config.shadow_fading_type            = 'none';
        LTE_config.cache_network                 = false;
        LTE_config.add_femtocells                = false;
        LTE_config.macroscopic_pathloss_is_model = true;
        ignore_MCL = true;
    elseif isfield(LTE_config, 'apply_MCL')
        if ~LTE_config.apply_MCL
            ignore_MCL = true;
        else
            ignore_MCL = false;
        end
    else
        ignore_MCL = false;
    end
    
    
    % If cache is on and the cache file exists, load network from disk
    if LTE_config.cache_network && exist(LTE_config.network_cache,'file') || LTE_config.usePathlossMapFromMemory
        if ~LTE_config.usePathlossMapFromMemory
            if LTE_config.debug_level>=1
                fprintf('Loading network from %s\n',LTE_config.network_cache);
            end
            load(LTE_config.network_cache);
            % Change map resolution to match the loaded maps'
            LTE_config.map_resolution = networkPathlossMap.data_res; %#ok<NODEF>
        else
            if LTE_config.debug_level>=1
                fprintf('Using pathloss data input as variable\n');
            end
            [sites,eNodeBs]    = utils.miscUtils.copySitesAndeNodeBs(varargin{1}{2}.eNodeBs,varargin{1}{2}.eNodeBs_sectors);
            networkPathlossMap = varargin{1}{2}.networkPathlossMap.clone();
            if isfield(varargin{1}{2},'networkShadowFadingMap')
                networkShadowFadingMap = varargin{1}{2}.networkShadowFadingMap.clone();
            end
        end
        
        % Rewrite the neighbors, just in case (it may lead to
        % difficult-to-trace bugs due to loaded networks)
        % Store the other eNodeBs as (potential) interferers
        for s_ = 1:length(eNodeBs)
            eNodeBs(s_).neighbors_eNodeB = eNodeBs([1:(s_-1) (s_+1):length(eNodeBs)]);
        end
        
        % Add the power separation. X% to signaling/pilots (always on) and the rest for data
        set_signaling_power_and_number_of_TX_antennas(LTE_config,sites);
    else
        % Generate network (eNodeBs and macroscopic pathloss)
        if LTE_config.debug_level>=1
            sprintf('Generating network\n');
        end
        switch LTE_config.network_source
            case 'generated'
                [sites, eNodeBs, networkPathlossMap] = network_generation.generated_network(LTE_config);
                % Other case here -> other sources, eg. Network planning tool
            case 'capesso'
                [sites, eNodeBs, networkPathlossMap] = network_generation.capesso_network(LTE_config);
                if LTE_config.debug_level>=1
                    fprintf('Using data from planning tool Capesso');
                end
            case 'fixed pathloss'
                [sites, eNodeBs, networkPathlossMap] = network_generation.fixed_pathloss_network(LTE_config);
                if LTE_config.debug_level>=1
                    fprintf('Using fixed pathloss network (defined in the config parameters)');
                end
            otherwise
                error([LTE_config.network_source ' network source not supported\n']);
        end
        
        % Option to add extra eNodeBs with omnidirectional antennas
        if LTE_config.add_femtocells
            [sites, eNodeBs] = network_generation.add_femtocells(LTE_config,sites,eNodeBs,networkPathlossMap );
        end

        % Option to add extra Remote Radio Heads (RRHs)
        if LTE_config.RRHs_enabled
            if LTE_config.RRHs_enabled
                network_generation.add_RRHs(LTE_config,eNodeBs,networkPathlossMap );
            end
        end
        RRHs  = [eNodeBs.RRHs];
        nRRHs = length(RRHs);

        % Add the power separation. X% to signaling/pilots (always on) and the rest for data
        set_signaling_power_and_number_of_TX_antennas(LTE_config,sites);
        
        % Store the other eNodeBs as (potential) interferers
        for s_ = 1:length(eNodeBs)
            eNodeBs(s_).neighbors_eNodeB = eNodeBs([1:(s_-1) (s_+1):length(eNodeBs)]);
        end
        
        % Generate shadow fading
        if LTE_config.macroscopic_pathloss_is_model
            if LTE_config.debug_level>=1
                fprintf('Generating shadow fading\n');
            end
            switch LTE_config.shadow_fading_type
                case 'claussen'
                    [LTE_config.roi_x,LTE_config.roi_y] = networkPathlossMap.valid_range;
                    if LTE_config.debug_level>=1
                        fprintf('Generating Claussen space-correlated shadow fading map, ');
                    end
                    
                    if ~LTE_config.decouple_site_shadow_fading_maps
                        fprintf('one map per site or RRH\n');
                        nShadowFadingMaps = length(sites)+nRRHs;
                    else
                        fprintf('one map per cell or RRH\n');
                        nShadowFadingMaps = length(eNodeBs)+nRRHs;
                    end
                    networkShadowFadingMap = channel_gain_wrappers.shadowFadingMapClaussen(...
                        LTE_config.shadow_fading_map_resolution,...
                        LTE_config.roi_x,...
                        LTE_config.roi_y,...
                        LTE_config.shadow_fading_n_neighbors,...
                        nShadowFadingMaps,...
                        LTE_config.shadow_fading_mean,...
                        LTE_config.shadow_fading_sd,...
                        LTE_config.r_eNodeBs,...
                        LTE_config.deactivate_claussen_spatial_correlation);
                    networkShadowFadingMap.oneMapPerSite = ~LTE_config.decouple_site_shadow_fading_maps;
                    
                    for rrh_=1:nRRHs
                        RRHs(rrh_).site_id = rrh_+length(sites);
                    end
                case 'none'
                    [LTE_config.roi_x, LTE_config.roi_y] = networkPathlossMap.valid_range;
                    if LTE_config.debug_level>=1
                        fprintf('Generating dummy shadow fading map (i.e. no shadow fading');
                    end
                    networkShadowFadingMap = channel_gain_wrappers.shadowFadingDummyMap(...
                        LTE_config.roi_x,...
                        LTE_config.roi_y,...
                        length(eNodeBs));
                otherwise
                    error('%s shadow fading type not supported. Only "claussen" and "none" supported',LTE_config.shadow_fading_type);
            end
        end
        
        %% Finally apply Minimum Coupling Loss (MCL) after having all of the generated pathlosses.
        %  After the antenna gain, as TS.942-900
        %  and according to the type of pathloss on each position (i.e.,
        %  femtos, macros, RRHs...
        %
        %  states: RX_PWR = TX_PWR – Max (pathloss – G_TX – G_RX, MCL)
        %  KNOWN ISSUE: this assumes that G_RX is 0dB, which will normally be the
        %  case for a mobile terminal. This would have to be moved to the link
        %  level model (UE) if G_RX is to be taken into account
        
        if ~ignore_MCL
            if LTE_config.debug_level>=1
                fprintf('Applying Minimum Coupling Loss to the final pathloss maps\n');
            end
            networkPathlossMap.apply_MCL(LTE_config);
        end
        
        %% Calculate the SINR for each sector based on the pathloss and maximum TX power alone
        if LTE_config.macroscopic_pathloss_is_model
            % With shadow fading
            [networkPathlossMap.capacity,...
                networkPathlossMap.SINR,...
                networkPathlossMap.sector_assignment,...
                networkPathlossMap.maxSINR_assignment,...
                networkPathlossMap.diff_SINR_dB,...
                networkPathlossMap.sector_sizes,...
                networkPathlossMap.sector_centers] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,sites,eNodeBs,networkShadowFadingMap);
            % Without shadow fading
            [networkPathlossMap.capacity2,...
                networkPathlossMap.SINR2,...
                networkPathlossMap.sector_assignment2,...
                networkPathlossMap.maxSINR_assignment,...
                networkPathlossMap.diff_SINR_dB2,...
                networkPathlossMap.sector_sizes2,...
                networkPathlossMap.sector_centers2] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,sites,eNodeBs);
        else
            % No shadow fading for the case where the pathloss does not
            % come from a model (e.g., capesso)
            [networkPathlossMap.capacity,...
                networkPathlossMap.SINR,...
                networkPathlossMap.sector_assignment,...
                networkPathlossMap.maxSINR_assignment,...
                networkPathlossMap.diff_SINR_dB,...
                networkPathlossMap.sector_sizes,...
                networkPathlossMap.sector_centers] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,sites,eNodeBs);
        end
        
        % Save network
        if LTE_config.cache_network
            try
                [pathstr, name, ext] = fileparts(LTE_config.network_cache);
                
                %% Calculate hashtag for the eNodeB positions and powers
                % data_hash.eNodeB_pos                  = [eNodeBs.pos];
                % data_hash.eNodeB_altitude             = [eNodeBs.altitude];
                % data_hash.eNodeB_site_type            = [eNodeBs.site_type];
                % data_hash.eNodeB_max_power            = [eNodeBs_sectors.max_power];
                % data_hash.eNodeB_tx_height            = [eNodeBs_sectors.tx_height];
                % data_hash.eNodeBs_electrical_downtilt = [eNodeBs_sectors.electrical_downtilt];
                % data_hash.eNodeBs_mechanical_downtilt = [eNodeBs_sectors.mechanical_downtilt];
                % hashtag_network_cache                 = sprintf('_%s',utils.hashing.DataHash(data_hash));
                hashtag_network_cache                 = '';
                
                LTE_config.network_cache = fullfile(pathstr,[name hashtag_network_cache ext]);
                
                if LTE_config.debug_level>=1
                    fprintf('Saving network to file: %s\n',LTE_config.network_cache);
                end
                
                if exist(LTE_config.network_cache,'file')
                    throw(MException('LTEsim:cacheExists', 'The cache file was concurrently generated during another simulation run'));
                end
                
                if exist('networkShadowFadingMap','var')
                    save(LTE_config.network_cache,'sites','networkPathlossMap','networkShadowFadingMap','eNodeBs');
                else
                    save(LTE_config.network_cache,'sites','networkPathlossMap','eNodeBs');
                end
            catch err
                fprintf('Network cache could not be saved. If needed, it will be generated again in the next run (%s).\n',err.message);
            end
        end
    end
    
    % Delete this variable to save memory
    networkPathlossMap.distances = [];
else
    % Generate dummy object from the trace information
    [sites, eNodeBs, networkPathlossMap] = network_generation.read_trace(LTE_config,varargin);
    
    % Add the power separation. X% to signaling/pilots (always on) and the rest for data
    set_signaling_power_and_number_of_TX_antennas(LTE_config,sites);
end

% Configure the case for zero-delay
if LTE_config.feedback_channel_delay==0
    for b_=1:length(sites)
        for s_=1:length(sites(b_).sectors)
            sites(b_).sectors(s_).zero_delay_feedback = true;
        end
    end
end

% Configure unquantized feedback
if LTE_config.unquantized_CQI_feedback
    for b_=1:length(sites)
        for s_=1:length(sites(1).sectors)
            sites(b_).sectors(s_).unquantized_CQI_feedback = true;
        end
    end
end

% To avoid error of missing return argument
if ~exist('networkShadowFadingMap','var')
    networkShadowFadingMap = [];
end

function set_signaling_power_and_number_of_TX_antennas(LTE_config,eNodeBs)
% Set signaling power
for b_=1:length(eNodeBs)
    for s_=1:length(eNodeBs(b_).sectors)
        % Power part assigned to data and signaling
        data_power      = eNodeBs(b_).sectors(s_).max_power * (1-LTE_config.signaling_ratio);
        signaling_power = eNodeBs(b_).sectors(s_).max_power * LTE_config.signaling_ratio;
        eNodeBs(b_).sectors(s_).max_power       = data_power;
        eNodeBs(b_).sectors(s_).signaling_power = signaling_power;
        LTE_config.scheduler_params.max_power   = data_power; % max data transmit power in Watts
        
        % Number of TX antennas
        eNodeBs(b_).sectors(s_).nTX       = LTE_config.nTX;
        eNodeBs(b_).sectors(s_).total_nTX = LTE_config.nTX;
        if ~isempty(eNodeBs(b_).sectors(s_).RRHs)
            RRHs = [eNodeBs(b_).sectors(s_).RRHs];
            eNodeBs(b_).sectors(s_).total_nTX = ...
                eNodeBs(b_).sectors(s_).total_nTX + sum([RRHs.nTX]);
        end
    end
end
