classdef naming
    % Groups some functions used to generate names.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2012 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
        the_date
        date_time_string
    end
    
    methods
        function obj = naming
            obj.the_date = clock;
            obj.date_time_string = sprintf('%04d%02d%02d_%02d%02d%02d',...
                obj.the_date(1),...                     % Date: year
                obj.the_date(2),...                     % Date: month
                obj.the_date(3),...                     % Date: day
                obj.the_date(4),...                     % Date: hour
                obj.the_date(5),...                     % Date: minutes
                floor(obj.the_date(6)));                % Date: seconds
        end
        
        function results_filename = results_file(obj,LTE_config)
            if strcmp(LTE_config.results_file,'auto')
                
                if LTE_config.frequency/1e9 >= 1
                    this_freq = sprintf('%3.2fGHz',LTE_config.frequency/1e9);
                else
                    this_freq = sprintf('%3.0fMHz',LTE_config.frequency/1e6);
                end
                
                this_bw = sprintf('%gfMHz',LTE_config.bandwidth/1e6);
                
                switch LTE_config.tx_mode
                    case 1
                        MIMO_mode = 'SISO';
                    case 2
                        MIMO_mode = sprintf('%dx%dTxD',LTE_config.nTX,LTE_config.nRX);
                    case 3
                        MIMO_mode = sprintf('%dx%dOLSM',LTE_config.nTX,LTE_config.nRX);
                    case 4
                        MIMO_mode = sprintf('%dx%dCLSM',LTE_config.nTX,LTE_config.nRX);
                    case 5
                        MIMO_mode = sprintf('%dx%dMU-MiMO',LTE_config.nTX,LTE_config.nRX);
                    case 6
                        MIMO_mode = sprintf('%dx%d-rank1-CLSM',LTE_config.nTX,LTE_config.nRX);
                    otherwise
                        MIMO_mode = sprintf('%dx%d_%d',LTE_config.nTX,LTE_config.nRX,LTE_config.tx_mode);
                end
                
                if LTE_config.runtime_precoding
                    runtime_precoding_string = '_runtime_precoding_';
                else
                    runtime_precoding_string = '_precomputed_precoding_';
                end
                if isfield(LTE_config, 'rep')
                    rep = sprintf('-rep%ld', LTE_config.rep);
                else
                    rep = '';
                end
                if isempty(LTE_config.output_filename_suffix)
                    suffix_to_use = [];
                else
                    suffix_to_use = ['_' LTE_config.output_filename_suffix];
                end
                if exist('labindex','builtin')    % does not exist in older Matlab versions
                    lab_ind = labindex;
                else
                    lab_ind = 1;
                end
                if LTE_config.parallel_network 
                    t = getCurrentTask();
                    if ~isempty(t)
                        par_ID = num2str(t.ID);
                    else
                        par_ID = '';
                    end
                else
                    par_ID = '';
                end
                
                % Output file name
                results_filename = sprintf('%s_freq_%s_bw_%s_%.1fKmph_%dTTIs_%s_lab%02.0f_%s_%s%s%s%s%s%s.mat',...
                    this_freq,...                             % Frequency
                    this_bw,...                               % System Bandwidth
                    LTE_config.channel_model.type,...         % Channel model userd
                    LTE_config.UE_speed*3.6,...               % Channel speed (Km/h)
                    LTE_config.simulation_time_tti,...        % Simulaton length
                    obj.date_time_string,...                  % Date string
                    lab_ind,...                               % Lab index (in order to avoid multiple labs writing to the same file)
                    strrep(LTE_config.scheduler,' ','_'),...  % Scheduler type
                    MIMO_mode,...                             % MIMO mode
                    runtime_precoding_string,...              % Whether runtime precoding is used. If not, it is included in the trace
                    LTE_config.release,...                    % Release number
                    suffix_to_use, ...                        % Optional suffix                    
                    rep,...
                    par_ID);                                  % ID of the current worker (for parallel simulations only) 
                
            else
                results_filename = LTE_config.results_file;
            end
            
            % Check for .mat extension
            if ~strcmp(results_filename((end-3):(end)),'.mat')
                results_filename = [results_filename '.mat'];
            end
        end
        
        function cache_filename = macroscopic_pathloss_cache(obj,LTE_config)
            if strcmp(LTE_config.network_cache,'auto')
                if LTE_config.frequency >= 1e9
                    this_freq = sprintf('%3.2fGHz',LTE_config.frequency/1e9);
                else
                    this_freq = sprintf('%3.0fMHz',LTE_config.frequency/1e6);
                end
                if strcmp(LTE_config.network_source, 'generated')
                    % Generated network layout
                    if isfield(LTE_config,'AntennaPattern3d') && LTE_config.AntennaPattern3d
                        downtilting = sprintf('_%d°-%d°',LTE_config.antenna.electrical_downtilt,LTE_config.antenna.mechanical_downtilt);
                    else
                        downtilting = [];
                    end
                    if LTE_config.add_femtocells
                        femto_string = '_plus_femtocells';
                    else
                        femto_string = [];
                    end
                    
                    if isfield(LTE_config,'macroscopic_pathloss_model_settings')
                        if isfield(LTE_config.macroscopic_pathloss_model_settings,'environment')
                            pathloss_model_string = sprintf('%s %s',LTE_config.macroscopic_pathloss_model,LTE_config.macroscopic_pathloss_model_settings.environment);
                        elseif isfield(LTE_config.macroscopic_pathloss_model_settings,'alpha')
                            pathloss_model_string = sprintf('%s %d',LTE_config.macroscopic_pathloss_model,LTE_config.macroscopic_pathloss_model_settings.alpha);
                        end
                    else
                        pathloss_model_string = LTE_config.macroscopic_pathloss_model;
                    end
                    pathloss_model_string = strrep(strtrim(pathloss_model_string),' ','_');
                    
                    if isfield(LTE_config,'nr_eNodeB_rings')
                        size_string = sprintf('%d_rings',LTE_config.nr_eNodeB_rings);
                    else
                        size_string = sprintf('size_%d', LTE_config.network_size);
                    end
                    
                    if LTE_config.RRHs_enabled
                        % Make it a bit shorter
                        hash_opts.Format = 'base64';
                        DAS_string       = sprintf('_DAS_%s',utils.hashing.DataHash(LTE_config.RRH,hash_opts));
                        
                        % To avoid path errors
                        DAS_string       = strrep(DAS_string,'/','f');
                        DAS_string       = strrep(DAS_string,'\','b');
                    else
                        DAS_string = [];
                    end
                    
                    if LTE_config.parallel_network
                        t = getCurrentTask();
                        if ~isempty(t)
                            par_ID = num2str(t.ID);
                        else
                            par_ID = '';
                        end
                    else
                        par_ID = '';
                    end
                    
                    cache_filename = fullfile(LTE_config.default_network_cache_folder,...
                        sprintf('network_%s_%d_sectors_%d°_offset_%dm_res_%s_%s_antenna%s_%s_freq%s%s',...
                        size_string,...
                        length(LTE_config.sector_azimuths),...
                        LTE_config.antenna_azimuth_offsett,...
                        LTE_config.map_resolution,...
                        pathloss_model_string,...
                        LTE_config.antenna.antenna_gain_pattern,...
                        downtilting,...
                        this_freq,...
                        femto_string,...
                        DAS_string,...
                        par_ID));
                else
                    % When using planning tool Capesso
                    cache_filename = fullfile(LTE_config.default_network_cache_folder,...
                        sprintf('network_%dm_res_%s_%s_freq_capesso',...
                        LTE_config.map_resolution,...
                        strrep(strtrim([LTE_config.macroscopic_pathloss_model ' ' LTE_config.macroscopic_pathloss_model_settings.environment]),' ','_'),...
                        this_freq));
                end
                
                % Add shadow fading part to the filename
                switch LTE_config.shadow_fading_type
                    case 'none'
                        cache_filename = sprintf('%s_no_shadow_fading',cache_filename);
                    case 'claussen'
                        cache_filename = sprintf('%s_claussen%gdB_shadow_fading',cache_filename,LTE_config.shadow_fading_sd);
                end
                % Finish network cache file generation
                cache_filename = sprintf('%s.mat',cache_filename);
            else
                cache_filename = LTE_config.network_cache;
            end
        end
        
        function cache_filename = UE_cache(obj,LTE_config)
            if LTE_config.parallel_network
                t = getCurrentTask();
                if ~isempty(t)
                    par_ID = num2str(t.ID);
                else
                    par_ID = '';
                end
            else
                par_ID = '';
            end
            if strcmp(LTE_config.UE_cache_file,'auto')
                switch LTE_config.network_source
                    case 'capesso'
                        cache_filename = fullfile(LTE_config.default_UE_cache_folder,...
                            sprintf('UE_cache_capesso_%s_%s_%s.mat',strrep(LTE_config.UE_distribution,' ','_'),obj.date_time_string,par_ID));
                    otherwise
                        cache_filename = fullfile(LTE_config.default_UE_cache_folder,...
                            sprintf('UE_cache_%s_%s_%s.mat',strrep(LTE_config.UE_distribution,' ','_'),obj.date_time_string,par_ID));
                end
            else
                cache_filename = LTE_config.UE_cache_file;
            end            
        end
        
        function LTE_config = channel_trace_cache(obj,LTE_config)
            tx_mode_string = utils.miscUtils.tx_mode_to_string(LTE_config.tx_mode);
            if isempty(tx_mode_string)
                error('TX mode %d not supported',LTE_config.tx_mode);
            end
            
            if isempty(LTE_config.channel_trace_id)
                LTE_config.extra_params_trace_id = [];
            else
                LTE_config.extra_params_trace_id = ['_' sprintf('%g',LTE_config.channel_trace_id)];
            end
            
            % Check if this is the Winner II channel model
            if strcmp(LTE_config.channel_model.type,'winner+')
                % Check if there are extra winner config params
                if ~isempty(LTE_config.winner_antenna_params)
                    LTE_config.extra_params_hash = ['_' utils.hashing.DataHash(LTE_config.winner_antenna_params)];
                else
                    LTE_config.extra_params_hash = [];
                end
            else
                LTE_config.extra_params_hash = [];
            end
            % This only makes sense for v2. However, it is pulled out, so
            % the config-string is always set
            if LTE_config.RRHs_enabled
                LTE_config.RRH_string = sprintf('_RRH-%dTX',LTE_config.RRH.nTX);
            else
                LTE_config.RRH_string = '';
            end
            
            if LTE_config.parallel_network
                t = getCurrentTask();
                if ~isempty(t)
                    par_ID = num2str(t.ID);
                else
                    par_ID = '';
                end
            else
                par_ID = '';
            end
            
            if strcmp(LTE_config.pregenerated_ff_file,'auto')
                switch LTE_config.trace_version
                    case 'v1'
                        cache_filename = fullfile(LTE_config.default_channel_trace_folder,...
                            sprintf('%s_%dx%d_%s_%gMHz_%.1fKmph_%ds_delay_%d_%s%s%s_%s',...
                            LTE_config.channel_model.type,...
                            LTE_config.nTX,...
                            LTE_config.nRX,...
                            tx_mode_string,...
                            LTE_config.bandwidth/1e6,...
                            LTE_config.UE_speed*3.6,...
                            LTE_config.channel_model.trace_length,...
                            LTE_config.feedback_channel_delay,...
                            LTE_config.trace_version,...
                            LTE_config.extra_params_hash,...
                            LTE_config.extra_params_trace_id,...
                            par_ID));
                    otherwise
                        
                        
                        cache_filename = fullfile(LTE_config.default_channel_trace_folder,...
                            sprintf('%s_BS-%dTX%s_UE-%dRX_%dMHz_%.1fKmph_%ds_%s%s%s_%s',...
                            LTE_config.channel_model.type,...
                            LTE_config.nTX,...
                            LTE_config.RRH_string,...
                            LTE_config.nRX,...
                            LTE_config.bandwidth/1e6,...
                            LTE_config.UE_speed*3.6,...
                            LTE_config.channel_model.trace_length,...
                            LTE_config.trace_version,...
                            LTE_config.extra_params_hash,...
                            LTE_config.extra_params_trace_id,...
                            par_ID));
                end
                
                if LTE_config.tx_mode==4
                    if LTE_config.wideband_precoding
                        cache_filename = [cache_filename 'WB_precoding'];
                    end
                end
            end
            LTE_config.pregenerated_ff_file = cache_filename;
        end
    end
end

