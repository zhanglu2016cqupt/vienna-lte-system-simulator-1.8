function [ UEs, extra_info ] = LTE_init_generate_users_and_add_schedulers(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,BLER_curves,networkClock)
% Get a cell assignment matrix. Here I will directly take it from the
% networkMacroscopicPathlossMap, but the following code could generate an equivalent
% one for any type of pathloss model.
%
% (c) Josep Colom Ikuno, Martin Taranetz, INTHFT, 2010

use_UE_cache      = LTE_config.UE_cache;
cache_file_exists = exist(LTE_config.UE_cache_file,'file');


%% Data needed also for plotting: generate every time
fprintf('Generating UEs: %s\n',LTE_config.UE_distribution);
switch LTE_config.UE_distribution
    case 'constant UEs per cell'
        UE_spatial_distribution = spatial_distributions.constantElementsPerCellSpatialDistribution(networkPathlossMap, LTE_config.UE_per_eNodeB);
    case 'traffic map'
        UE_spatial_distribution = spatial_distributions.UeTrafficMapSpatialDistribution(networkPathlossMap, LTE_config.rescale_factor, LTE_config.traffic_map_config);
    case 'radial'
        UE_spatial_distribution = spatial_distributions.radialSpatialDistribution(networkPathlossMap,sites,LTE_config.UE_distribution_radii,LTE_config.UE_distribution_nUEs,LTE_config.UE_distribution_overlap_radial_UEs);
    case 'predefined'
        UE_spatial_distribution = spatial_distributions.predefinedSpatialDistribution(networkPathlossMap, LTE_config.UE_positions);
    case 'constant UEs per ROI'
        UE_spatial_distribution = spatial_distributions.constantElementsPerRoiSpatialDistribution(networkPathlossMap, LTE_config.nUEs);
    case 'trace'
        UE_spatial_distribution = spatial_distributions.predefinedSpatialDistribution(zeros(length(networkPathlossMap.pathloss),2));
    otherwise
        error('UE distribution %s not supported.',LTE_config.UE_distribution);
end
UE_positions = UE_spatial_distribution.generate_positions;
% UE_positions
%% Creating or loading UE position, depending on the configuration
% Create UEs according to the previously generated positions
UEs = network_elements.UE;
if (~use_UE_cache) || (use_UE_cache&&~cache_file_exists)
    for u_ = 1:size(UE_positions,1)
        % General UE settings that can be saved and re-used
        UEs(u_)     = network_elements.UE;
        UEs(u_).id  = u_;
        UEs(u_).pos = UE_positions(u_,:);
        
        % Generate a walking model for the user
        switch LTE_config.UE_distribution
            case 'trace'
                % Do not add walking model and set UE to trace mode
                UEs(u_).trace_UE      = true;
                UEs(u_).deactivate_UE = true; % By default deactivate UE
            otherwise
                % Add walking model and set UE to NOT-trace mode
                UEs(u_).trace_UE = false;
                if isfield(LTE_config.UE,'walk') && strcmp(LTE_config.UE.walk,'SLvsLL')
                    UEs(u_).walking_model = walking_models.SLvsLLWalkingModel; % Since no angle is specified, a random one is chosen
                else
                    UEs(u_).walking_model = walking_models.straightWalkingModel(LTE_config.UE_speed*LTE_config.TTI_length); % Since no angle is specified, a random one is chosen
                end
        end
        
        UEs(u_).runtime_precoding = LTE_config.runtime_precoding;
    end
    if use_UE_cache
        try
            if exist(LTE_config.UE_cache_file,'file')
                throw(MException('LTEsim:cacheExists', 'The cache file was concurrently generated during another simulation run'));
            end
            if LTE_config.debug_level>=1
                fprintf('Saving UE positions to %s\n',LTE_config.UE_cache_file);
            end
            save(LTE_config.UE_cache_file,'UEs','UE_positions');
        catch err
            fprintf('UE cache could not be saved. If needed, it will be generated again in the next run (%s).\n',err.message);
        end
    end
else
    % Load UEs
    if LTE_config.debug_level>=1
        fprintf('Loading UE positions from %s\n',LTE_config.UE_cache_file);
    end
    UE_data_from_cache = load(LTE_config.UE_cache_file);
    UEs = UE_data_from_cache.UEs;
    
    % Get a struct name with everything except the UEs field
    names = fieldnames(UE_data_from_cache);
    if length(names) > 1 % If there is nothing more besides the 'UEs' data
        extra_info = struct;
        for field_idx=1:length(names)
            current_field = names{field_idx};
            switch current_field
                case 'UEs'
                    % Do nothing
                otherwise
                    extra_info.(current_field) = UE_data_from_cache.(current_field);
            end
        end
    end
end

%% Steps needed for FFR
if LTE_config.FFR_active
    switch LTE_config.shadow_fading_type
        case 'none'
            % OK
        otherwise
            error('Right now, due to how the R3 frequency assignment is realized, just simulations WITHOUT shadow fading are supported for FFR simulations');
    end
    
    site_pos                 = reshape([sites.pos],2,[])';
    ROI_center               = [mean(networkPathlossMap.roi_x) mean(networkPathlossMap.roi_y)];
    distance_to_center       = sqrt(sum([site_pos(:,1)-ROI_center(1) site_pos(:,2)-ROI_center(2)].^2,2));
    [~, min_distance_eNodeB] = min(distance_to_center);
    
    % Frequency assignment
    LTE_config.scheduler_params.frequency_assignment = utils.ffrUtils.assign_frequencies_to_hex_grid(sites(min_distance_eNodeB).sectors(1).eNodeB_id,eNodeBs,networkPathlossMap.sector_centers);
    
    if ~LTE_config.FFR_override
        tx_mode_string_long = utils.miscUtils.tx_mode_to_string_long(LTE_config.tx_mode,LTE_config.nTX,LTE_config.nRX);
        [optimum_BFR, optimum_FR_SINR_switching_dB, ~] = utils.ffrUtils.load_FFR_optimum_BFRs(tx_mode_string_long);
    else
        optimum_BFR                  = LTE_config.FFR_params.beta_FR;
        optimum_FR_SINR_switching_dB = LTE_config.FFR_params.SINR_threshold_value;
    end
    
    %% Calculate frequency assignment
    FFR_UE_mapping = utils.ffrUtils.assign_FFR_band_to_UEs(UEs,optimum_FR_SINR_switching_dB,networkPathlossMap);
    
    LTE_config.scheduler_params.FFR_UE_mapping = FFR_UE_mapping;
    LTE_config.scheduler_params.beta_FR        = optimum_BFR;
end

%% Initialise schedulers
LTE_init_add_schedulers(LTE_config,sites,UEs,CQI_mapper,BLER_curves);

%% Other UE initialization, including adding a downlink and uplink channel object to each user
% The downlink will contain pathloss maps, so depending on the user's position, it will 'see' a certain pathloss.
% Add also the penetration loss and noise figure.
% The uplink is simply a delay between the UE and the eNodeB.

for u_=1:length(UEs)
    
    % Add wall loss, which is experienced by interferers.
    if LTE_config.add_femtocells
        UEs(u_).wall_loss = LTE_config.femtocells_config.macroscopic_pathloss_model_settings.penetration_loss;
    else
        UEs(u_).wall_loss = 0;
    end
    
    % Add receiver antenna gain
    UEs(u_).antenna_gain = LTE_config.UE.antenna_gain;
    
    % Add noise figure
    UEs(u_).receiver_noise_figure = LTE_config.UE.receiver_noise_figure;
    
    % Thermal noise (receiver) for the link quality model (in linear: watts)
    UEs(u_).thermal_noise_W_RB = 10^(0.1*LTE_config.UE.thermal_noise_density)/1000 * LTE_config.RB_bandwidth * 10^(UEs(u_).receiver_noise_figure/10);
    
    % Default tx mode for feedback (for the old trace format -v1- this
    % sets the only tx mode that can be used)
    UEs(u_).default_tx_mode = LTE_config.tx_mode;
    
    % LTE precoding codebook for all TX modes
    UEs(u_).codebook = phy_modeling.miscUtils.get_all_precoding_combinations;
    
    % Set signaling channel (eNodeB to UE)
    UEs(u_).eNodeB_signaling = network_elements.eNodebSignaling;
    
    % Number of RX antennas
    UEs(u_).nRX = LTE_config.nRX;
    
    % Set BLER curves for ACK/NACK calculation
    UEs(u_).BLER_curves = BLER_curves;
    
    % Clock
    UEs(u_).clock = networkClock;
    
    % CQI mapper
    UEs(u_).CQI_mapper = CQI_mapper;
    
    % Configure unquantized feedback
    UEs(u_).unquantized_CQI_feedback = LTE_config.unquantized_CQI_feedback;
    
    % Configure extra tracing
    UEs(u_).trace_SINR = LTE_config.trace_SINR;
    
    % Adaptive RI
    UEs(u_).adaptive_RI = LTE_config.adaptive_RI;
end



%% Safeguard against having no UEs
if length(UEs)==1 && isempty(UEs(1).id)
    no_UEs = true;
else
    no_UEs = false;
end

if ~no_UEs
    % Assign the UEs to their nearest (in terms of SINR) eNodeB and assign some extra parameters
    sector_UE = false(1,3);
    for u_ = 1:length(UEs)        
        if ~UEs(u_).trace_UE
            % Attach UE to eNodeB
            [ site_id, sector_num, eNodeB_id] = networkPathlossMap.cell_assignment(UEs(u_).pos); %#ok<ASGLU>
            eNodeBs(eNodeB_id).attachUser(UEs(u_));
        else
            % Do not attach the UEs yet (do that when we move and/or activate them)
        end
        
        if LTE_config.add_femtocells
            % Macro covered users experience penetration loss (dB) from outdoor to indoor environment. Wall loss for femtocells is already incorporated  in the dual slope model
            if  strcmp(UEs(u_).attached_site.site_type,'macro')
                UEs(u_).penetration_loss = LTE_config.femtocells_config.macroscopic_pathloss_model_settings.penetration_loss;
            else
                % Non-femto case
                UEs(u_).penetration_loss = 0;
            end
        else
            UEs(u_).penetration_loss = LTE_config.additional_penetration_loss;
        end
        
        % Check whether this UE should be deactivated to speed-up simulation
        if ~isempty(LTE_config.compute_only_UEs_from_this_eNodeBs)
            if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id==LTE_config.compute_only_UEs_from_this_eNodeBs,1))
                % Deactivate UE
                UEs(u_).deactivate_UE = true;
            else
                % Activate UE (already activated by default, but just in case)
                UEs(u_).deactivate_UE = false;
            end
        end
        
        % Append traffic model to users
        if isfield(LTE_config.traffic_models,'type') && strcmp(LTE_config.traffic_models.type,'MLaner')
            if ~LTE_config.traffic_models.traffic_patterns || (UEs(u_).attached_site.id == 5 && ~sector_UE(UEs(u_).attached_sector)) % fullbuffer UE from center eNodeB
                UEs(u_).traffic_model = LTE_trafficmodel(LTE_config.traffic_models,UEs(u_),max(LTE_config.feedback_channel_delay,0),7,0);
                sector_UE(UEs(u_).attached_sector) = true;
            else
                UEs(u_).traffic_model = LTE_trafficmodel(LTE_config.traffic_models,UEs(u_),max(LTE_config.feedback_channel_delay,0),7,1,LTE_config.simulation_time_tti);
            end
        else
            UEs(u_).traffic_model = LTE_trafficmodel(LTE_config.traffic_models,UEs(u_),max(LTE_config.feedback_channel_delay,0));
        end
        end
    end
    
    % Choose as many points per cell as users
    if ~exist('extra_info','var')
        extra_info = [];
    end
end

function LTE_init_add_schedulers(LTE_config,sites,UEs,CQI_mapper,BLER_curves)
% Adds the needed scheduler type and resource block grid to each eNodeb's
% sector
% (c) Josep Colom Ikuno, INTHFT, 2008
% input:   eNodeBs  ... array of eNodeBs
%          UEs      ... array of UEs

% Check whether this scheduler exists
schedulers.schedulerFactory.check_whether_scheduler_is_defined(LTE_config.scheduler);

if LTE_config.debug_level>=1
    switch LTE_config.scheduler
        case 'FFR'
            fprintf('Creating %s scheduler (FR: %s, PR: %s) and resource block grids\n',LTE_config.scheduler,LTE_config.scheduler_params.FR_scheduler.scheduler,LTE_config.scheduler_params.PR_scheduler.scheduler);
        otherwise
            fprintf('Creating %s schedulers and resource block grids\n',LTE_config.scheduler);
    end
end

% No reason to use a different SINR averager instance for each scheduler, we can reuse the same one
switch LTE_config.SINR_averaging.algorithm
    case 'MIESM'
        the_SINR_averager = utils.miesmAveragerFast(LTE_config,LTE_config.SINR_averaging.BICM_capacity_tables,LTE_config.SINR_averaging.betas);
    case 'EESM'
        error('EESM SINR averaging is no longer supported supported');
    otherwise
        error('SINR averaging algorithm not supported');
end

% Add RB grid representation and scheduler to each sector.
% Set also homogeneous power load
for b_ = 1:length(sites)
    for s_=1:length(sites(b_).sectors)
        
        % Set whether the eNodeBs will always transmit, even if no UEs are attached.
        sites(b_).sectors(s_).always_on = LTE_config.always_on;
        
        max_data_power  = sites(b_).sectors(s_).max_power;
        
        LTE_config.scheduler_params.max_power       = max_data_power; % For backwards compatibility
        LTE_config.scheduler_params.CQI_params      = LTE_config.CQI_params;
        LTE_config.scheduler_params.default_tx_mode = LTE_config.tx_mode;
        if ~isfield(LTE_config.scheduler_params, 'cqi_reduction')
            LTE_config.scheduler_params.cqi_reduction = 0;
        end
        % RB grid creation and initialization
        switch LTE_config.tx_mode
            case 5
                sites(b_).sectors(s_).RB_grid = network_elements.resourceBlockGridMU(LTE_config.N_RB,LTE_config.sym_per_RB_nosync,LTE_config.sym_per_RB_sync, LTE_config.nTX);
            otherwise
                sites(b_).sectors(s_).RB_grid = network_elements.resourceBlockGrid(LTE_config.N_RB,LTE_config.sym_per_RB_nosync,LTE_config.sym_per_RB_sync);
        end
        sites(b_).sectors(s_).RB_grid.set_homogeneous_power_allocation(sites(b_).sectors(s_).max_power,sites(b_).sectors(s_).signaling_power);
        
        % Continue with Scheduler initialization
        sites(b_).sectors(s_).scheduler = schedulers.schedulerFactory.create_scheduler(LTE_config.scheduler,LTE_config.scheduler_params,sites(b_).sectors(s_));
        
        % Set scheduler SINR averaging algorithm
        sites(b_).sectors(s_).scheduler.set_SINR_averager(the_SINR_averager);

        % Other data required to perform SINR averaging at the transmitter side
        sites(b_).sectors(s_).scheduler.set_CQI_mapper(CQI_mapper);
        sites(b_).sectors(s_).scheduler.set_BLER_curves(BLER_curves);
        
        % Add genie information
        sites(b_).sectors(s_).scheduler.set_genie_UEs(UEs);
        sites(b_).sectors(s_).scheduler.set_genie_eNodeBs(sites);
        
        % Add TTI delay information
        sites(b_).sectors(s_).scheduler.set_feedback_delay_TTIs(LTE_config.feedback_channel_delay);
        
        % Add codebook
        sites(b_).sectors(s_).scheduler.runtime_precoding = LTE_config.runtime_precoding;
        if LTE_config.runtime_precoding
            sites(b_).sectors(s_).scheduler.Rel8_codebook     = phy_modeling.miscUtils.get_all_precoding_combinations;
        end
    end
end
