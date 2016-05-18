classdef capessoFilesReader < handle
    % Wraps functionalities needed to read Capesso pathloss maps
    % (c) Josep Colom Ikuno , INTHFT, 2012
    
    properties
        capesso_params
    end    
    
    methods
        function obj = capessoFilesReader(capesso_params)
            obj.capesso_params  = capesso_params;
        end
        function [sites eNodeBs] = read_cell_data(obj)
            % Reads out the specified folder and searchs for the following files:
            %   - *.par files
            %   - LTE_BTS_trial_cluster_20100809.txt
            %   - LTE_sites_trial_cluster.txt
            %   - LTE_cell_trial_cluster.txt
            %
            %   Text files have a different format in newer Atoll version
            %   For files of old Atoll version use function
            %   LTE_init_read_atoll_cell_data
            %
            % (c) Josep Colom Ikuno, Martin Taranetz, INTHFT, 2010
            
            capesso_params = obj.capesso_params;
            
            %% Parameter definitions
            folder                        = capesso_params.pathloss_data_folder;
            cell_atoll_filename           = capesso_params.cell_atoll_filename;
            site_atoll_filename           = capesso_params.site_atoll_filename;
            transmitter_atoll_filename    = capesso_params.transmitter_atoll_filename;
            kathrein_antenna_folder       = capesso_params.kathrein_antenna_folder;
            
            % Workaround for using frequency defined in LTE_load_params
            if capesso_params.frequency == 2.14e+9
                frequency = 2140;
            else
                warning('Selected frequency for the antenna gain pattern has been substituted for convinience for 2.14 GHz');
            end
            
            %% Read out LTE_cell_trial_cluster
            % The fields are
            % - Name
            % - Transmitter
            % - Max Power (dBm)
            % - Frequency Band
            % - Channel Number
            % - Traffic Load (UL) (%)
            % - Traffic Load (DL) (%)
            % - LTE Equipment
            % - Comments
            % - Max Number of Users
            % - UL Noise Rise (dB)
            % - Max number of intra-technology neighbours
            % - Max number of inter-technology neighbours
            % - Active	AMS & MU-MIMO Threshold (dB)
            % - Scheduler
            % - Reference Signal C/N Threshold (dB)
            % - Diversity Support (DL)
            % - TDD Frame Configuration
            % - Physical Cell ID
            % - Min Reuse Distance (m)
            % - Physical Cell ID Status
            % - Channel Allocation Status
            % - Max Traffic Load (UL) (%)
            % - Max Traffic Load (DL) (%)
            % - Inter-technology UL Noise Rise (dB)
            % - Inter-technology DL Noise Rise (dB)
            % - Diversity Support (UL)
            % - MU-MIMO Capacity Gain (UL)
            % - SS & PBCH EPRE Offset / RS (dB)
            % - PDSCH & PDCCH EPRE Offset / RS (dB)
            % - Layer (Lowest Layer = Highest Priority)
            % - Interference Coordination Support
            % - ICIC Ratio (DL) (%)
            % - ICIC Delta Path Loss Threshold (dB)
            % - Fractional Power Control Factor
            % - Max UL Noise Rise (dB)
            % - Max PUSCH C/(I+N) (dB)
            % - ICIC UL Noise Rise (dB)
            % - RS_EPRE	PBCH_POWER_OFFSET
            % - PDCCH_POWER_OFFSET
            % - PSS ID
            % - SSS ID
            % - Instantaneous Reference Signal Power (dBm)
            % - Instantaneous SS & PBCH Power (dBm)
            % - Average PDSCH & PDCCH Power (dBm)
            
            
            % Since the numerical fields are separated by a COMMA instead of a point, I
            % will have to do the 4,7 -> 4.7 conversion by hand
            fid_ = fopen (fullfile(folder,cell_atoll_filename),'r');
            % M = textscan(fid_,'%s %s %n %s %n %n %n %n %n %s %n %n %s %n %n %n %n %n %s %s %n %n %n %n %s %s %n %n %n %s %s %n %n %s %s %n %n %n %n %n %n','Delimiter','\t','Headerlines',1);
            cell_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
            fclose(fid_);
            
            %% Read LTE_sites_trial_cluster
            
            % The fields are:
            % - Name
            % - X
            % - Y
            % - Altitude (m)
            % - Pylon Height (m)
            % - Support Type
            % - SITENAME
            % - PlanningRegion
            % - RMRegion
            % - SITEARCHITECTURE
            % - Master
            % - SPATIAL_REGION
            % - MCC
            % - MNC
            
            fid_ = fopen (fullfile(folder,site_atoll_filename),'r');
            site_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
            fclose(fid_);
            
            %% Read LTE_BTS_trial_cluster_20100809
            %The fields are:
            % - Site
            % - Transmitter
            % - Active
            % - DX (m)
            % - DY (m)
            % - Polarisation
            % - Antenna
            % - Height (m)
            % - Azimuth (°)
            % - Mechanical Downtilt (°)
            % - Main Calculation Radius (m)
            % - Main Propagation Model
            % - Transmission Loss (dB)
            % - Reception Loss (dB)
            % - Noise Figure (dB)
            % - Hexagon groups
            % - Hexagon radius (m)
            % - Comments
            % - TMA Equipment
            % - Feeder Equipment
            % - BTS Equipment
            % - Transmission Feeder Length (m)
            % - Reception Feeder Length (m)
            % - Receiver antenna diversity gain (dB)
            % - Miscellaneous Transmission Losses (dB)
            % - Miscellaneous Reception Losses (dB)
            % - Extended Calculation Radius (m)
            % - Extended Propagation Model
            % - Main Resolution (m)
            % - Extended Resolution (m)
            % - Additional Electrical Downtilt (°)
            % - Number of Transmission Antenna Ports
            % - Number of Reception Antenna Ports
            % - Transmitter Type
            % - COVERAGELAYER
            % - BORDERCOORDINATION
            % - Diversity
            % - CoverageReported
            % - PatternType
            % - Optical	PARAMITRISATION_CLASS
            
            fid_ = fopen (fullfile(folder,transmitter_atoll_filename),'r');
            transmitter_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
            fclose(fid_);
            
            %% Generate eNodeBs
            
            % Now read the site and transmitter description and put it in the corresponding objects
            % TODO : Change to fprintf
            display('Creating eNodeBs');
            sites = network_elements.eNodeB;
            
            for site_idx = 1:size(site_atoll{1},1)
                
                id        = site_idx;
                name      = site_atoll{1}{site_idx};
                pos       = [ str2double(strrep(site_atoll{2}{site_idx},',','.')) str2double(strrep(site_atoll{3}{site_idx},',','.')) ]; % pos in meters (x,y)
                altitude  = str2double(site_atoll{4}{site_idx}(2:end-1));
                site_name = site_atoll{7}{site_idx};
                
                % Search for this eNodeB's sectors
                sectors_idx = search_string_in_textscan_output(transmitter_atoll{1},name);
                num_sectors = length(sectors_idx);
                
                % Fill in eNode object data
                sites(site_idx)           = network_elements.eNodeB;
                sites(site_idx).id        = id;
                sites(site_idx).name      = name;
                sites(site_idx).pos       = pos;
                sites(site_idx).altitude  = altitude;
                sites(site_idx).site_name = site_name;
                sites(site_idx).site_type = 'macro';
                % neighbor info will be filled in a posteriori
                
                for sector_idx = 1:num_sectors
                    if sector_idx==1
                        sites(site_idx).sectors             = network_elements.eNodeB_sector;
                    else
                        sites(site_idx).sectors(sector_idx) = network_elements.eNodeB_sector;
                    end
                    sites(site_idx).sectors(sector_idx).parent_eNodeB = sites(site_idx);
                    
                    transmitter_atoll_idx = sectors_idx(sector_idx);
                    transmitter           = transmitter_atoll{2}{transmitter_atoll_idx};
                    cell_idx              = find(strcmp(transmitter,cell_atoll{1}),1,'first');
                    if isempty(cell_idx)
                        error('The cell and transmitter Atoll files do not match. %s not found',transmitter);
                    end
                    switch capesso_params.planning_tool
                        case 'atoll'
                            % Do nothing
                        case 'capesso'
                            % Change the transmitter to read always the first sector
                            % (we are using an isotrop antenna, so very probably we
                            % will only have the pathloss map for one sector: no need
                            % for the others).
                            transmitter = regexprep(transmitter, '/U.*', '/U1');
                        otherwise
                            error('Only "atoll" and "capesso" supported');
                    end
                    
                    antenna_name          = transmitter_atoll{7}{transmitter_atoll_idx};
                    height                = str2double(strrep(transmitter_atoll{8}{transmitter_atoll_idx},',','.'));
                    azimuth               = str2double(strrep(transmitter_atoll{9}{transmitter_atoll_idx},',','.'));
                    %azimuth = 0;
                    mechanical_downtilt   = str2double(strrep(transmitter_atoll{10}{transmitter_atoll_idx},',','.'));
                    % Additional information contained in antenna_name
                    antenna_name_split = textscan(antenna_name, '%s %s %f %f %f %s', 'delimiter', '|');
                    antenna_type              = antenna_name_split{1}{1};
                    % antenna_polarization    = antenna_name_split{2};
                    % antenna_3dB_beamwidth_h = antenna_name_split{3};
                    % antenna_3dB_beamwidth_v = antenna_name_split{4};
                    electrical_downtilt       = antenna_name_split{5};
                    % antenna_frequency_band  = antenna_name_split{6};
                    % Kathrein antenna used for this sector
                    antenna               = antennas.kathreinTSAntenna(kathrein_antenna_folder, antenna_type, frequency);
                    power_dBm             = str2double(cell_atoll{3}{cell_idx});
                    
                    % Fill in eNodeB Sector Object data
                    sites(site_idx).sectors(sector_idx).parent_eNodeB       = sites(site_idx);
                    sites(site_idx).sectors(sector_idx).id                  = sector_idx;
                    sites(site_idx).sectors(sector_idx).antenna             = antenna;
                    sites(site_idx).sectors(sector_idx).max_power           = 10^(power_dBm/10) / 1000;
                    sites(site_idx).sectors(sector_idx).transmitter         = transmitter;
                    
                    sites(site_idx).sectors(sector_idx).antenna_name        = antenna_name;
                    sites(site_idx).sectors(sector_idx).antenna_type        = antenna_type;
                    sites(site_idx).sectors(sector_idx).tx_height           = height;
                    sites(site_idx).sectors(sector_idx).azimuth             = azimuth;
                    sites(site_idx).sectors(sector_idx).mechanical_downtilt = mechanical_downtilt;
                    sites(site_idx).sectors(sector_idx).electrical_downtilt = electrical_downtilt;
                end
            end
            
            eNodeBs = sites(1).sectors(1); % Initialization
            eNodeBs_idx = 1;
            for b_=1:length(sites)
                for s_=1:length(sites(b_).sectors)
                    eNodeBs(eNodeBs_idx) = sites(b_).sectors(s_);
                    eNodeBs(eNodeBs_idx).eNodeB_id = eNodeBs_idx;
                    eNodeBs_idx = eNodeBs_idx + 1;
                end
            end
            
            %% Nested utility function
            function indexes = search_string_in_textscan_output(a_column,a_string)
                indexes = [];
                for i_=1:size(a_column,1)
                    if strcmp(a_column{i_},a_string)
                        indexes = [indexes i_];
                    end
                end
            end
        end
        
        function DTM = read_DTM(obj)
            % Martin Taranetz, INTHFT 2009
            
            % Function returns elevation map and description data with coordinates,
            % resolution and matrix size
            
            %% Initial Parameters
            enable_plotting     = obj.capesso_params.debug_plotting;
            dtm_folder          = obj.capesso_params.dtm_folder;
            dtm_file_name       = obj.capesso_params.dtm_file_name;
            dtm_hdr_file_name   = obj.capesso_params.dtm_hdr_file_name;
            
            %% Read out the DTM .bil file and generate elevation map
            dtmfile_    = fopen(fullfile(dtm_folder, dtm_file_name), 'r');
            dtmhdrfile_ = fopen(fullfile(dtm_folder, dtm_hdr_file_name), 'r');
            
            dtm_hdr = textscan(dtmhdrfile_, '%s %f');
            
            DTM.description.NWxmap = dtm_hdr{2}(1);
            DTM.description.NWymap = dtm_hdr{2}(2);
            DTM.description.xdim   = dtm_hdr{2}(3);
            DTM.description.ydim   = dtm_hdr{2}(4);
            DTM.description.ncols  = dtm_hdr{2}(5);
            DTM.description.nrows  = dtm_hdr{2}(6);
            DTM.description.nbits  = dtm_hdr{2}(7);
            DTM.description.nbands = dtm_hdr{2}(8);
            
            DTM.description.SWxmap = DTM.description.NWxmap;                                                      % Lower-leftmost corner (x) -> SW
            DTM.description.SWymap = DTM.description.NWymap - DTM.description.ydim*(DTM.description.nrows-1/100); % Lower-leftmost corner (y) -> SW
            
            DTM.description.NExmap = DTM.description.NWxmap + DTM.description.xdim*(DTM.description.ncols-1/100);
            DTM.description.NEymap = DTM.description.NWymap;
            
            DTM.description.SExmap = DTM.description.NExmap;
            DTM.description.SEymap = DTM.description.SWymap;
            
            DTM.description.roi_x = [DTM.description.SWxmap DTM.description.SExmap];
            DTM.description.roi_y = [DTM.description.SWymap DTM.description.NWymap];
            
            if DTM.description.nbits == 16
                % Read in Matlabs columnwise manner and transpose
                dtm_map_t = fread(dtmfile_,[DTM.description.ncols, DTM.description.nrows],'int16');
            end
            DTM.data  = flipud(dtm_map_t');
            
            if enable_plotting
                figure;
                imagesc(DTM.description.roi_x,DTM.description.roi_y,DTM.data);
                hold on;
                set(gca,'YDir','normal');
                xlabel('x pos [m]');
                ylabel('y pos [m]');
                colorbar;
                title('Elevation map (m)');
            end
            
            fclose(dtmfile_);
            fclose(dtmhdrfile_);
        end
        
        function [pathlossmaps unused_maps_filenames] = read_pathloss_data(obj,sites,elevation_map_full)
            % Reads out data files from capesso. Outputs a matrix containing the
            % pathloss maps (eNodeB_idx,sector_idx) in the same order as in the eNodeBs
            % object matrix. For capesso file the matrix is of size
            % (length(eNodeBs),1), as only one pathloss file per eNodeB is needed. The
            % second output logs interpreted pathloss files that were not output (the
            % eNodeB/site was not on the list or only one sector was needed, thus the
            % others were discarded).
            %
            % (c) Martin Taranetz, Josep Colom Ikuno INTHFT, 2010
            
            %% Parameters
            capesso_params = obj.capesso_params;
            folder         = capesso_params.pathloss_data_folder;
            
            % Enable / disable plots of the pathloss maps
            enable_debug_plotting     = obj.capesso_params.debug_plotting;
            
            %% Read Capesso Files and generate pathloss map
            
            % .clos files are renamed .zip files containing the .los files
            % If there are any .clos files in the directory, they are shifted to the
            % zip_files directory and renamed to .zip
            % In a second step, the .zip files are extracted back in the folder of
            % origin, which also contains the .par description files
            
            clos_files = dir(fullfile(folder,'*.clos'));
            if ~isempty(clos_files)
                mkdir(folder, 'zip_files');
                for i_=1:1:length(clos_files)
                    zipfilename = [strtok(clos_files(i_).name,'.'),'.zip'];
                    movefile(fullfile(folder,clos_files(i_).name),[folder,'/zip_files/',zipfilename]);
                    unzip([folder,'/zip_files/',zipfilename], folder);
                end
            end
            
            % unzipped .los files are ready for readout
            los_files = dir(fullfile(folder,'*.los'));
            
            unused_maps_filenames = [];
            
            for i_=1:1:length(los_files)
                file_name = los_files(i_).name;
                plfile_ = fopen(fullfile(folder, file_name),'r');
                %Search for .par file and get description for .los file
                pathlossmap.description = get_pathloss_files_description(file_name, folder);
                
                PATHLOSS_MAP_A = fread(plfile_,'int16');
                PATHLOSS_MAP_A = reshape(PATHLOSS_MAP_A,pathlossmap.description.ncols, pathlossmap.description.nrows);
                
                PATHLOSS_MAP_B = PATHLOSS_MAP_A.';
                PATHLOSS_MAP_B = PATHLOSS_MAP_B./16;
                
                pathlossmap.pathloss_map(:,:)= PATHLOSS_MAP_B;
                
                % if enable_debug_plotting
                if enable_debug_plotting
                    coloraxis = [100,200];
                    plot_pathloss_map(pathlossmap,coloraxis);
                end
                
                % Search for the index where to save this pathloss map (it is the
                % same order as how the eNodeBs are stored)
                filename_to_compare = strtok(strrep(pathlossmap.description.filename,'#','/'),'.');
                eNodeB_idx = [];
                sector_idx = [];
                for b_=1:length(sites)
                    switch capesso_params.planning_tool
                        case 'atoll'
                            for s_=1:length(sites(b_).sectors)
                                site_name      = sites(b_).sectors(1).transmitter;
                                site_name_file = sites(b_).sectors(1).transmitter;
                                if strcmp(site_name,site_name_file)
                                    eNodeB_idx = b_;
                                    sector_idx = s_;
                                    break
                                end
                            end
                        case 'capesso'
                            site_name      = strtok(sites(b_).sectors(1).transmitter,'/');
                            site_name_file = strtok(filename_to_compare,'/');
                            if strcmp(site_name,site_name_file)
                                eNodeB_idx = b_;
                                sector_idx = 1;
                                break
                            end
                    end
                    
                end
                
                % Save the pathloss file in the corresponding position
                if isempty(eNodeB_idx) || isempty(sector_idx)
                    % Do not assign this pathloss to the output (no eNodeB/sector is using it)
                    unused_maps_filenames{length(unused_maps_filenames)+1} = pathlossmap.description.filename;
                else
                    % This site is not used
                    if ~exist('pathlossmaps')
                        pathlossmaps(eNodeB_idx,sector_idx) = pathlossmap;
                    else
                        % This sector is not used
                        if size(pathlossmaps,1)>=eNodeB_idx && size(pathlossmaps,2)>sector_idx
                            if ~isempty(pathlossmaps(eNodeB_idx,sector_idx).description.filename)
                                unused_maps_filenames{length(unused_maps_filenames)+1} = pathlossmap.description.filename;
                            else
                                pathlossmaps(eNodeB_idx,sector_idx) = pathlossmap;
                            end
                        else
                            pathlossmaps(eNodeB_idx,sector_idx) = pathlossmap;
                        end
                    end
                end
                
            end
            
            fclose('all');
            
            % Assign the eNodeBs' sectors to the pathloss maps. Since these are isotrop maps, if
            % we have maps for more than one sector (eg. W85A_32_t2#2FU1 and
            % W85A_32_t2#2FU2), we will just take the #2FU (which means '/') and discard the rest.
            
            function color_axis = plot_pathloss_map(pathlossmap_plot,varargin)
                % In :
                % pathlossmap_: containing data and description
                % i_ for assignment to description data and subplot
                %
                % The extra argument allows you to specify the scale of the plot.
                % Nevertheless, the used scale is returned by the function so you can
                % use the same for all of your plots
                
                figure_steps_col = 5;
                figure_steps_row = 5;
                pncols_ = pathlossmap_plot.description.ncols;
                pnrows_ = pathlossmap_plot.description.nrows;
                pxdim_  = pathlossmap_plot.description.xdim;
                pydim_  = pathlossmap_plot.description.ydim;
                
                % Attention : 3 and 4 are chosen static
                % change, if more plots required
                
                figure_capesso_los_file = figure('Colormap',[1 0 0;1 0.06275 0;1 0.1294 0;1 0.1922 0;1 0.2588 0;1 0.3216 0;1 0.3882 0;1 0.451 0;1 0.5176 0;1 0.5804 0;1 0.6471 0;1 0.7098 0;1 0.7725 0;1 0.8392 0;1 0.902 0;1 0.9686 0;0.9686 1 0;0.902 1 0;0.8392 1 0;0.7725 1 0;0.7098 1 0;0.6471 1 0;0.5804 1 0;0.5176 1 0;0.451 1 0;0.3882 1 0;0.3216 1 0;0.2588 1 0;0.1922 1 0;0.1294 1 0;0.06275 1 0;0 1 0;0 1 0.06275;0 1 0.1255;0 1 0.1882;0 1 0.251;0 1 0.3137;0 1 0.3765;0 1 0.4392;0 1 0.502;0 1 0.5608;0 1 0.6235;0 1 0.6863;0 1 0.749;0 1 0.8118;0 1 0.8745;0 1 0.9373;0 1 1;0 0.9373 1;0 0.8745 1;0 0.8118 1;0 0.749 1;0 0.6863 1;0 0.6235 1;0 0.5608 1;0 0.498 1;0 0.4392 1;0 0.3765 1;0 0.3137 1;0 0.251 1;0 0.1882 1;0 0.1255 1;0 0.06275 1;0 0 1]);
                subplot(1,1, 1,'Parent',figure_capesso_los_file);
                
                imagesc(pathlossmap_plot.description.roi_x,pathlossmap_plot.description.roi_y,pathlossmap_plot.pathloss_map(:,:));
                set(gca,'YDir','normal');
                xlabel('x pos [m]');
                ylabel('y pos [m]');
                if length(varargin)>=1
                    color_axis = varargin{1};
                    caxis(color_axis);
                else
                    color_axis = caxis;
                end
                colorbar;
                title(['Capesso .LOS File: ',strrep(pathlossmap_plot.description.filename,'_','\_')]);
            end
            
            function pathloss_files_description = get_pathloss_files_description(los_file_name, folder)
                % Description of .LOS stored in .PAR Files
                
                par_files = dir(fullfile(folder,'*.par'));
                par_index = 0;
                
                % Search for corresponding .PAR File in Folder
                for j_=1:1:length(par_files)
                    if strcmp(strtok(los_file_name,'.'),strtok(par_files(j_).name,'.'))
                        par_index = j_;
                        break;
                    end
                end
                
                % Read out corresponding .PAR file
                if(par_index)
                    par_fid = fopen(fullfile(folder,par_files(par_index).name),'r');
                    c = textscan(par_fid, '%s %f', 'delimiter','=');
                    
                    descr.filename = los_file_name;
                    descr.NWxmap   = c{2}(1);
                    descr.NWymap   = c{2}(2);
                    descr.nrows    = c{2}(3);
                    descr.ncols    = c{2}(4);
                    descr.xdim     = c{2}(5);
                    descr.ydim     = c{2}(6);
                    descr.masked   = c{2}(7);
                    
                    descr.SWxmap = descr.NWxmap;                            % Lower-leftmost corner (x) -> SW
                    descr.SWymap = descr.NWymap - descr.ydim*(descr.nrows-1/100); % Lower-leftmost corner (y) -> SW
                    
                    descr.NExmap = descr.NWxmap + descr.xdim*(descr.ncols-1/100);
                    descr.NEymap = descr.NWymap;
                    
                    descr.SExmap = descr.NExmap;
                    descr.SEymap = descr.SWymap;
                    
                    descr.roi_x = [descr.SWxmap descr.SExmap];
                    descr.roi_y = [descr.SWymap descr.NWymap];
                    
                    fclose(par_fid);
                    
                    pathloss_files_description = descr;
                else
                    % TODO ... Return value
                    if LTE_config.debug_level>=1
                        fprintf('Could not find .par file for %s\n', los_file_name);
                    end
                end
            end
        end
    end
end

