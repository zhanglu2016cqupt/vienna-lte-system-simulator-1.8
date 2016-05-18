% LTE System Level simulator
% 
% (c) Josep Colom Ikuno, INTHFT, 2008
%

fprintf('Vienna LTE System Level simulator\n');
fprintf('(c) 2008, INTHFT, TU Wien\n');
fprintf(' This work has been funded by A1 Telekom Austria AG and the Christian Doppler Laboratory for Design Methodology of Signal Processing Algorithms.\n\n');
fprintf('  By using this simulator, you agree to the license terms stated in the license agreement included with this work\n\n');

%% Load and plot BLER curves
[ BLER_curves CQI_mapper ] = LTE_init_load_BLER_curves;

%% Get the eNodeBs, the Macroscopic Pathloss and the Shadow Fading
% No need to generate shadow fading maps when using network planning tool
if strcmp(LTE_config.network_source, 'generated')
    [eNodeBs eNodeBs_sectors networkPathlossMap networkShadowFadingMap] = LTE_init_network_generation;  % Squeezed network generation code in this script for ease of reuse (cell capacity calculation)
else
    [eNodeBs eNodeBs_sectors networkPathlossMap] = LTE_init_network_generation;
end

%% Calculate target sector
LTE_config.target_sector = LTE_common_get_target_sector(eNodeBs,networkPathlossMap);

%% Calculate average cell capacity. NOTE: sector_SINR denotes the max SINR
%% in the ROI
if exist('networkShadowFadingMap','var')
    [LTE_config.capacity_no_shadowing sector_SINR_no_shadowing] = LTE_common_calculate_cell_capacity(networkPathlossMap,eNodeBs,CQI_mapper);
    [LTE_config.capacity              sector_SINR             ] = LTE_common_calculate_cell_capacity(networkPathlossMap,eNodeBs,CQI_mapper,networkShadowFadingMap);
else
    [LTE_config.capacity              sector_SINR             ] = LTE_common_calculate_cell_capacity(networkPathlossMap,eNodeBs,CQI_mapper);
end

%% Plot network
if LTE_config.show_network>0
    if exist('networkShadowFadingMap','var')
        LTE_plot_sector_SINR_cdfs(sector_SINR,sector_SINR_no_shadowing);
        LTE_plot_loaded_network(eNodeBs,networkPathlossMap,networkShadowFadingMap);
    else
        LTE_plot_sector_SINR_cdfs(sector_SINR);
        % When using Capesso - use own plot function for network
        if ~isfield(LTE_config,'show_capesso_network')
            LTE_plot_loaded_network(eNodeBs,networkPathlossMap);
        end
    end
end

%% Add a clock to each network element
networkClock = network_elements.clock(LTE_config.TTI_length);

for b_=1:length(eNodeBs)
    eNodeBs(b_).clock = networkClock;
end

%% Create users (UEs)
[UEs extra_UE_cache_info] = LTE_init_generate_users(eNodeBs,networkPathlossMap);

if isempty(UEs)
    warning('No UEs generated. Simulation will be skipped and no results file will be saved.');
else
    
    %% Generate/load the fast fading traces
    if isempty(strfind(LTE_config.pregenerated_ff_file,'.mat'))
        LTE_config.pregenerated_ff_file = [LTE_config.pregenerated_ff_file '.mat'];
    end
    ff_file_exists = exist(LTE_config.pregenerated_ff_file,'file');
    if LTE_config.recalculate_fast_fading || (~ff_file_exists && ~LTE_config.recalculate_fast_fading)
        % Generated UE fast fading
        print_log(1,sprintf('Generating UE fast fading and saving to %s\n',LTE_config.pregenerated_ff_file));
        pregenerated_ff = LTE_init_get_microscale_fading_SL_trace;
        save(LTE_config.pregenerated_ff_file,'pregenerated_ff','-v7.3');
    else
        % Load UE fast fading
        if LTE_config.reuse_pregenerated_ff_trace_from_last_run
            print_log(1,'Reusing UE fast fading trace from the previous run.\n');
            if ~exist('pregenerated_ff','var')
                error('pregenerated_ff variable should exist. Not found.\n');
            end
        else
            print_log(1,sprintf('Loading UE fast fading from %s \n.',LTE_config.pregenerated_ff_file));
            load(LTE_config.pregenerated_ff_file,'pregenerated_ff');
        end
        
        % Wrong number of nTX or nRX antennas
        if LTE_config.nTX~=pregenerated_ff.nTX || LTE_config.nRX~=pregenerated_ff.nRX
            error('Trace is for a %dx%d system. Config defines a %dx%d system.',pregenerated_ff.nTX,pregenerated_ff.nRX,LTE_config.nTX,LTE_config.nRX);
        end
        
        % Wrong bandwidth case
        if LTE_config.bandwidth ~= pregenerated_ff.system_bandwidth
            error('Loaded FF trace is not at the correct frequency: %3.2f MHz required, %3.2f MHz found',LTE_config.bandwidth/1e6,pregenerated_ff.system_bandwidth/1e6);
        end
        
        % Wrong UE speed case (not applicable if trace UE_speed is NaN->speed independent)
        if (pregenerated_ff.UE_speed~=LTE_config.UE_speed) && ~isnan(pregenerated_ff.UE_speed)
            error('Loaded FF trace is generated at %3.2f m/s. UE speed is %3.2f m/s. Trace cannot be used.',pregenerated_ff.UE_speed,LTE_config.UE_speed);
        end
        
        % Wrong channel type
        if ~strcmpi(pregenerated_ff.channel_type,LTE_config.channel_model.type)
            error('Loaded FF trace is not for the correct channel type. %s found, %s required.',pregenerated_ff.channel_type,LTE_config.channel_model.type);
        end
        
        % Print microscale fading trace speed
        if isnan(pregenerated_ff.UE_speed)
            print_log(1,sprintf('Microscale fading trace is speed-independent\n'));
        else
            print_log(1,sprintf('UE Fast fading trace at %3.2f m/s (%3.2f Km/h)\n',pregenerated_ff.UE_speed,pregenerated_ff.UE_speed*3.6));
        end
        
    end
    
    % The SINR averaging algorithm to be used. It does not make sense that
    % different UEs would use a different one, so it is the same one (saves
    % some memory)
    switch LTE_config.SINR_averaging.algorithm
        case 'MIESM'
            the_SINR_averager = utils.miesmAveragerFast(LTE_config.SINR_averaging.BICM_capacity_tables,LTE_config.SINR_averaging.betas);
        case 'EESM'
            the_SINR_averager = utils.eesmAverager(LTE_config.SINR_averaging.betas,LTE_config.SINR_averaging.MCSs);
        otherwise
            error('SINR averaging algorithm not supported');
    end
    
    % Add a downlink and uplink channel object to each user
    % The downlink will contain pathloss maps, so depending on the user's position, it will 'see' a certain pathloss.
    % Add also the penetration loss and noise figure.
    % The uplink is simply a delay between the UE and the eNodeB.
    for u_=1:length(UEs)
        
        % Add penetration loss
        UEs(u_).penetration_loss = LTE_config.additional_penetration_loss;
        
        % Add receiver antenna gain
        UEs(u_).antenna_gain = LTE_config.UE.antenna_gain;
        
        % Add noise figure
        UEs(u_).receiver_noise_figure = LTE_config.UE.receiver_noise_figure;
        
        % Add downlink channel (includes macroscopic pathloss, shadow fading and fast fading models)
        UEs(u_).downlink_channel = channel_models.downlinkChannelModel(UEs(u_));
        
        % Macroscopic pathloss
        UEs(u_).downlink_channel.set_macroscopic_pathloss_model(networkPathlossMap);
        
        % Shadow fading (data obtained from planning tools already have this information incorporated)
        if ~strcmp(LTE_config.network_source,'capesso')
            UEs(u_).downlink_channel.set_shadow_fading_model(networkShadowFadingMap);
        end
        
        % Set fast fading from the eNodeB to an attached UE.
        UEs(u_).downlink_channel.set_fast_fading_model_model(channel_gain_wrappers.fastFadingWrapper(pregenerated_ff,'random',length(eNodeBs),length(eNodeBs(1).sectors)));
        
        % Set UE SINR averaging algorithm
        UEs(u_).SINR_averager = the_SINR_averager;
        
        % Set signaling channel (eNodeB to UE)
        UEs(u_).eNodeB_signaling = network_elements.eNodebSignaling;
        
        % Number of RX antennas
        UEs(u_).nRX = LTE_config.nRX;
        
        % Thermal noise in dBm
        UEs(u_).downlink_channel.thermal_noise_watts_RB = 10^(0.1*LTE_config.UE.thermal_noise_density)/1000 * LTE_config.RB_bandwidth;
        UEs(u_).downlink_channel.thermal_noise_dBW_RB   = 10*log10(UEs(u_).downlink_channel.thermal_noise_watts_RB);
        
        % Set BLER curves for ACK/NACK calculation
        UEs(u_).BLER_curves = BLER_curves;
        
        % Uplink channel
        UEs(u_).uplink_channel = channel_models.uplinkChannelModel(...
            UEs(u_),...
            LTE_config.N_RB,...
            LTE_config.maxStreams,...
            LTE_config.feedback_channel_delay);
        UEs(u_).clock = networkClock;
        UEs(u_).CQI_mapper = CQI_mapper;
        
        % Configure unquantized feedback
        if LTE_config.unquantized_CQI_feedback
            UEs(u_).unquantized_CQI_feedback = true;
        end
    end
    
    %% Initialise schedulers
    LTE_init_add_schedulers(eNodeBs,UEs,CQI_mapper,BLER_curves);
    
    % Cache the RB_grid object in the UE object, so as to avoid too many
    % calls to the function. This would have to be taken into account when
    % implementing handover
    for u_=1:length(UEs)
        UEs(u_).RB_grid = UEs(u_).downlink_channel.RB_grid;
    end
    
    %% Initialise the tracing
    % Global traces
    simulation_traces = tracing.simTraces;
    % Traces from received UE feedbacks (eNodeB side)
    simulation_traces.eNodeB_rx_feedback_traces = tracing.receivedFeedbackTrace(...
        LTE_config.simulation_time_tti,...
        length(UEs),...
        LTE_config.N_RB,...
        LTE_config.maxStreams,...
        LTE_config.traces_config.unquantized_CQI_feedback);
    
    for b_=1:length(eNodeBs)
        for s_=1:length(eNodeBs(b_).sectors)
            eNodeBs(b_).sectors(s_).feedback_trace = simulation_traces.eNodeB_rx_feedback_traces;
            
            % Scheduler trace
            scheduler_trace = tracing.schedulerTrace(LTE_config.simulation_time_tti);
            eNodeBs(b_).sectors(s_).scheduler.trace   = scheduler_trace;
            simulation_traces.scheduler_traces{b_,s_} = scheduler_trace;
        end
    end
    % eNodeB traces
    simulation_traces.eNodeB_tx_traces = tracing.enodebTrace(eNodeBs(1),UEs(1).downlink_channel.RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti);
    for b_=2:length(eNodeBs)
        simulation_traces.eNodeB_tx_traces(b_) = tracing.enodebTrace(eNodeBs(b_),UEs(1).downlink_channel.RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti);
    end
    % UE traces
    for u_=1:length(UEs)
        UEs(u_).trace = tracing.ueTrace(LTE_config.simulation_time_tti,LTE_config.N_RB,LTE_config.maxStreams,LTE_config.traces_config,LTE_config.latency_time_scale,LTE_config.TTI_length);
        if u_==1
            simulation_traces.UE_traces = UEs(u_).trace;
        else
            simulation_traces.UE_traces(u_) = UEs(u_).trace;
        end
    end
    
    %% Give the schedulers access to the UE traces
    % Then they can make decisions base on their received throughput. More
    % complex and realistic solutions may be possible, but then the eNodeB
    % should dynamically allocate resources to store UE-related data (and then
    % release once the UE is not attached to it anymore). It is easier like this :P
    for b_=1:length(eNodeBs)
        for s_=1:length(eNodeBs(b_).sectors)
            eNodeBs(b_).sectors(s_).scheduler.UE_traces = simulation_traces.UE_traces;
        end
    end
    
    %% Print all the eNodeBs
    % Print the eNodeBs (debug)
    print_log(2,'eNodeB List\n');
    if LTE_config.debug_level >=2
        for b_=1:length(eNodeBs)
            eNodeBs(b_).print;
        end
    end
    print_log(2,'\n');
    
    %% Print all the Users
    print_log(2,'User List\n');
    if LTE_config.debug_level >=2
        for u_=1:length(UEs)
            UEs(u_).print;
        end
    end
    print_log(2,'\n');
    
    
    %% Main simulation loop
    print_log(1,['Entering main simulation loop, ' num2str(LTE_config.simulation_time_tti,'%5.0f') ' TTIs\n']);
    
    % Inititialize timer
    tic;
    starting_time = toc;
    
    % Network clock is initialised to 0
    while networkClock.current_TTI < LTE_config.simulation_time_tti
        % First of all, advance the network clock
        networkClock.advance_1_TTI;
        % Print all the eNodeBs and UEs from their absolute coordinates
        if LTE_config.show_network>1 || (networkClock.current_TTI==1 && LTE_config.show_network>0)
            LTE_plot_show_network(eNodeBs,UEs,LTE_config.map_resolution,networkClock.current_TTI);
            %         waitforbuttonpress
        end
        
        % Move users. To improve
        for u_ = 1:length(UEs)
            UEs(u_).move;
            [ x_range y_range ] = networkPathlossMap.valid_range;
            
            % If user went outside of ROI, relocate him somewhere else. Beam me
            % up, Scotty! Take me somewhere in the map!!
            if ~UEs(u_).is_in_roi(x_range,y_range)
                new_UE_position = networkPathlossMap.random_position;
                old_UE_position = UEs(u_).pos;
                
                old_eNodeB_id = UEs(u_).attached_eNodeB.id;
                % Actually it should not be done like this. Measure all the
                % neighboring cells' SNR and then decide which one is better
                [new_eNodeB_id new_eNodeB_sector] = networkPathlossMap.cell_assignment(new_UE_position);
                
                % Teleport UE
                UEs(u_).pos = new_UE_position;
                
                % Deattach UE from old eNodeB and reattach to new one
                UEs(u_).start_handover(eNodeBs(new_eNodeB_id),new_eNodeB_sector);
                
                if LTE_config.show_network>1
                    scatter(new_UE_position(1),new_UE_position(2),'Marker','o','MarkerEdgeColor','green');
                    pause(0.1)
                end
                
                % Print some debug
                print_log(2,['TTI ' num2str(networkClock.current_TTI) ': UE ' num2str(UEs(u_).id) ' going out of ROI, teleporting to ' num2str(new_UE_position(1)) ' ' num2str(new_UE_position(2)) '. eNodeB ' num2str(old_eNodeB_id) ' -> eNodeB ' num2str(new_eNodeB_id) '\n']);
            end
        end
        
        if  ~mod(networkClock.current_TTI,500)
            LTE_plot_show_network(eNodeBs,UEs,LTE_config.map_resolution,networkClock.current_TTI);
            %         waitforbuttonpress
        end
        
        % Placing the link measurement model here allows us to simulate transmissions with 0 delay
        if LTE_config.feedback_channel_delay==0
            for u_ = 1:length(UEs)
                % Measure SINR and prepare CQI feedback
                UEs(u_).link_quality_model(LTE_config);
            end
        end
        
        % The eNodeBs receive the feedbacks from the UEs
        for s_ = 1:length(eNodeBs_sectors)
            % Receives and stores the received feedbacks from the UEs
            eNodeBs_sectors(s_).receive_UE_feedback;
            % Schedule users
            eNodeBs_sectors(s_).schedule_users;
        end
        
        % For the non 0-delay case, call link quality model
        if LTE_config.feedback_channel_delay~=0
            for u_ = 1:length(UEs)
                % Measure SINR and prepare CQI feedback
                UEs(u_).link_quality_model(LTE_config);
            end
        end
        
        % Call link performance model (evaluates whether TBs are received
        % corretly according to the information conveyed by the link quality
        % model. Additionally, send feedback (channel quality indicator +
        % ACK/NACK)
        for u_ = 1:length(UEs)
            UEs(u_).link_performance_model;
            UEs(u_).send_feedback;
        end
        
        %fprintf('.');
        if mod(networkClock.current_TTI,50)==0
            elapsed_time = toc;
            time_per_iteration = elapsed_time / networkClock.current_TTI;
            estimated_time_to_finish = (LTE_config.simulation_time_tti - networkClock.current_TTI)*time_per_iteration;
            estimated_time_to_finish_h = floor(estimated_time_to_finish/3600);
            estimated_time_to_finish_m = estimated_time_to_finish/60 - estimated_time_to_finish_h*60;
            fprintf('Time to finish: %3.0f hours and %3.2f minutes\n',estimated_time_to_finish_h,estimated_time_to_finish_m);
        end
    end
    
    print_log(1,'Simulation finished\n');
    
    if ~isempty(extra_UE_cache_info)
        simulation_traces.extra_UE_info = extra_UE_cache_info;
    end
    
    print_log(1,['Saving results to ' LTE_config.results_file '\n']);
    
    if LTE_config.compact_results_file
        all_UE_traces = [simulation_traces.UE_traces];
        all_eNodeB_traces = [simulation_traces.eNodeB_tx_traces.sector_traces];
        extra_UE_info = simulation_traces.extra_UE_info;
        clear the_UE_traces;
        for s_=1:length(all_eNodeB_traces)
            the_eNodeB_traces(s_).acknowledged_data = all_eNodeB_traces(s_).acknowledged_data;
        end
        for u_=1:length(all_UE_traces)
            the_UE_traces(u_).TB_size     = all_UE_traces(u_).TB_size;
            the_UE_traces(u_).ACK         = all_UE_traces(u_).ACK;
            the_UE_traces(u_).TB_CQI      = all_UE_traces(u_).TB_CQI;
            the_UE_traces(u_).nCodewords  = all_UE_traces(u_).nCodewords;
            the_UE_traces(u_).CQI_sent    = all_UE_traces(u_).CQI_sent;
            the_UE_traces(u_).TB_SINR_dB  = all_UE_traces(u_).TB_SINR_dB;
            the_UE_traces(u_).position    = all_UE_traces(u_).position;
            the_UE_traces(u_).ACK         = all_UE_traces(u_).ACK;
            the_UE_traces(u_).SINR_RS_not = all_UE_traces(u_).SINR_RS_not;
        end
        save(LTE_config.results_file,'the_UE_traces','the_eNodeB_traces','extra_UE_info');
        print_log(1,'Only UE some UE and eNodeB traces saved (compact results file)\n');
    else
        % Some options to save space
        if LTE_config.delete_ff_trace_at_end
            pregenerated_ff.traces = [];
        end
        
        if LTE_config.delete_pathloss_at_end
            networkPathlossMap.pathloss = [];
            networkPathlossMap.sector_assignment = [];
            networkPathlossMap.sector_assignment_no_shadowing = [];
            if exist('networkShadowFadingMap','var')
                networkPathlossMap.pathloss = [];
            else
                % Do nothing
            end
        end
        save(LTE_config.results_file);
        print_log(1,'Traces saved in the standard format\n');
    end
end
