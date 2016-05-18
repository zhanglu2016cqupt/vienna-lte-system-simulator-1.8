function LTE_config = LTE_load_params(varargin)

%% Default simulation option
default_simulation = 'tri_sector_tilted';

%% If no simulation option is defined, use the default one
if ~isempty(varargin)
    simulation_type = varargin{1};
else
    simulation_type = default_simulation; %  Default value
end

fprintf('Using "%s" simulation configuration.\n',simulation_type)

%% Load the corresponding simulation parameters
switch simulation_type
    case 'tri_sector'
        LTE_config = simulation_config.tri_sector.apply_parameters;
    case {'tri_sector_tilted', 'tri_sector_tilted_4x2', 'tri_sector_tilted_4x4', 'stochastic_tri_sector_tilted'}
        LTE_config = simulation_config.hex_grid_tilted.apply_parameters;
    case 'tri_sector_plus_femtocells'
        LTE_config = simulation_config.hex_grid_tilted_with_femtocells.apply_parameters;
    case 'six_sector_tilted'
        LTE_config = simulation_config.hex_grid_sixsectors.apply_parameters;
    case 'capesso_pathlossmaps'
        LTE_config = simulation_config.example_capesso.apply_parameters;
    case 'omnidirectional_eNodeBs'
        LTE_config = simulation_config.hex_grid_omnidirectional.apply_parameters;
    case 'tri_sector_tilted_traffic'
        LTE_config = simulation_config.hex_grid_tilted_traffic.apply_parameters;
    case 'LLvsSL'
        LTE_config = simulation_config.LLvsSL.apply_parameters;
    case 'trace'
        LTE_config = simulation_config.trace.apply_parameters;
    case 'TPvsSNR'
        LTE_config = simulation_config.TPvsSNR.apply_parameters;       
    otherwise
        warning('Simulation type not defined: using default one instead');
        simulation_type = default_simulation;
        LTE_config = simulation_config.hex_grid_tilted.apply_parameters;
end

%% Some adjustments to the loaded simulation parameters
switch simulation_type
    case 'tri_sector'
        LTE_config.results_file                = 'auto';
        LTE_config.output_filename_suffix      = 'tri_sector';
    case 'tri_sector_tilted'
        LTE_config.results_file                       = 'auto';
        LTE_config.output_filename_suffix             = 'tri_sector_tilted';
        % Uncommenting these lines would select the innermost site and its
        % adjacent ring of sites only. Of course, in the case where no extra
        % sites (e.g., femtocells) are present. In general, it saves time
        % in hex-grid scenarios to have it uncommented, but it may lead to
        % confusion in non-hex-grid scenarios, so it is thus left
        % commented.
        %if LTE_config.support_MBSFN == true
            %LTE_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
            %LTE_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
            %LTE_config.compute_only_UEs_from_this_eNodeBs = [13 20 31 32 33 48];
            %LTE_config.default_shown_GUI_cells            = [13 20 31 32 33 48];
        %end
    case 'stochastic_tri_sector_tilted'
        LTE_config.results_file                = 'auto';
        LTE_config.output_filename_suffix      = 'stochastic_tri_sector_tilted';
        % UEs to compute and show are determined at run-time.
    case 'tri_sector_plus_femtocells'
        LTE_config.results_file                = 'auto';
        LTE_config.output_filename_suffix      = 'tri_sector_plus_femtocells';
    case 'six_sector_tilted'
        LTE_config.results_file                = 'auto';
        LTE_config.outpu5t_filename_suffix     = 'six_sector';
    case 'capesso_pathlossmaps'
        LTE_config.results_file                = 'auto';
        LTE_config.output_filename_suffix      = 'capesso_pathloss';
    case 'tri_sector_tilted_4x2'
        LTE_config.results_file  = 'auto';
        LTE_config.nTX           = 4;
        LTE_config.nRX           = 2;
        LTE_config.channel_model.trace_length  = 10;
        LTE_config.output_filename_suffix      = 'tri_sector_tilted_4x2';
        LTE_config.pregenerated_ff_file        = 'auto';
    case 'tri_sector_tilted_4x4'
        LTE_config.results_file  = 'auto';
        LTE_config.nTX           = 4;
        LTE_config.nRX           = 4;
        LTE_config.channel_model.trace_length  = 5;
        LTE_config.output_filename_suffix      = 'tri_sector_tilted_4x4';
        LTE_config.compute_only_UEs_from_this_eNodeBs = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
        LTE_config.default_shown_GUI_cells            = [13 14 15 16 17 18 19 20 21 28 29 30 31 32 33 34 35 36 46 47 48];
    case 'omnidirectional_eNodeBs'
        LTE_config.results_file                = 'auto';
        LTE_config.output_filename_suffix      = 'omnidirectional_eNodeBs';
end