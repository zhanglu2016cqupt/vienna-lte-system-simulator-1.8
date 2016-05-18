classdef antennaFileImporter
    % Reads out a given folder/s containing antenna .txap OR .msi files and stores
    % them in a searchable way. Basically the same as .msiFileImporter but
    % considering the new header of the .txap files
    % (c) Josep Colom Ikuno, Martin Taranetz, INTHFT, 2010
    % www.nt.tuwien.ac.at
    
    properties
        antenna_data
    end
    
    methods
        function obj=antennaFileImporter(dirname,varargin)
            data_idx = 1;
            for folder_idx = 1:1+length(varargin)
                if folder_idx==1
                    dirname_to_dir = dirname;
                else
                    dirname_to_dir = varargin{folder_idx-1}
                end
                
                % Read .txap files
                antenna_files = dir(fullfile(dirname_to_dir,'*.txap'));
                for i_=1:length(antenna_files)
                    file_name_ = antenna_files(i_).name;
                    full_filename = fullfile(dirname_to_dir,file_name_);
                    fid_ = fopen(full_filename,'r');
                    
                    % Read header of antenna file
                    M1 =  textscan(fid_,'%s %s',1,'Delimiter',' ');
                    M2 =  textscan(fid_,'%s %s',2,'Delimiter',' ');
                    M3 =  textscan(fid_,'%s %s',1,'Delimiter',' ');
                    M4 =  textscan(fid_,'%s %n',1,'Delimiter',' ');
                    M5 =  textscan(fid_,'%s %n %s',1,'Delimiter',' ');
                    
                    % Store data into matrices
                    antenna_name_split = textscan(M1{2}{1}, '%n %s %f %f %f %s', 'delimiter', '|');
                    
                    % Horizontal antenna gain pattern
                    MP = textscan(fid_,'%f %f',360,'Delimiter',' ','Headerlines',2);
                    
                    % Vertical antenna gain pattern
                    MP2 = textscan(fid_,'%f %f',360,'Delimiter',' ','Headerlines',2);
                    
                    obj.antenna_data(data_idx).id               = antenna_name_split{1};
                    obj.antenna_data(data_idx).name             = M1{2};
                    % obj.antenna_data(data_idx).frequency       = M3{2};
                    obj.antenna_data(data_idx).frequency        = 2140; % STATIC - FURTHER DATA NEEDED !!!
                    obj.antenna_data(data_idx).max_antenna_gain = M5{2};
                    obj.antenna_data(data_idx).electrical_tilt  = M4{2};
                    
                    % Store different gain patterns columnwise
                    obj.antenna_data(data_idx).horizontal_degree = MP{1};
                    obj.antenna_data(data_idx).horizontal_gain_pattern = MP{2};
                    obj.antenna_data(data_idx).vertical_degree = MP2{1};
                    obj.antenna_data(data_idx).vertical_gain_pattern = MP2{2};
                    
                    fclose(fid_);
                    data_idx = data_idx+1;
                end
                
                % Read .msi files
                antenna_files = dir(fullfile(dirname_to_dir,'*.msi'));
                for i_=1:length(antenna_files)
                    file_name_ = antenna_files(i_).name;
                    full_filename = fullfile(dirname_to_dir,file_name_);
                    fid_ = fopen(full_filename,'r');
                    
                    % Read header of antenna file
                    M1 =  textscan(fid_,'%s %n',1,'Delimiter',' ');
                    M2 =  textscan(fid_,'%s %s',1,'Delimiter',' ');
                    M3 =  textscan(fid_,'%s %n',1,'Delimiter',' ');
                    M4 =  textscan(fid_,'%s %n',1,'Delimiter',' ');
                    M5 =  textscan(fid_,'%s %n %s',1,'Delimiter',' ');
                    M6 =  textscan(fid_,'%s %s %s %s %s %s %s  %s %s',1,'Delimiter',' ');
                    M7 =  textscan(fid_,'%s %n %s',1,'Delimiter',' ');
                    
                    % Horizontal antenna gain pattern
                    MP = textscan(fid_,'%f %f',360,'Delimiter',' ','Headerlines',2);
                    
                    % Vertical antenna gain pattern
                    MP2 = textscan(fid_,'%f %f',360,'Delimiter',' ','Headerlines',2);
                    
                    % Store data into matrices
                    obj.antenna_data(data_idx).id               = M1{2};
                    obj.antenna_data(data_idx).name             = M2{2};
                    obj.antenna_data(data_idx).frequency        = M3{2};
                    obj.antenna_data(data_idx).max_antenna_gain = M5{2};
                    obj.antenna_data(data_idx).electrical_tilt  = M7{2};
                    
                    % Store different gain patterns columnwise
                    obj.antenna_data(data_idx).horizontal_degree       = MP{1};
                    obj.antenna_data(data_idx).horizontal_gain_pattern = MP{2};
                    obj.antenna_data(data_idx).vertical_degree         = MP2{1};
                    obj.antenna_data(data_idx).vertical_gain_pattern   = MP2{2};
                    
                    fclose(fid_);
                    data_idx = data_idx+1;
                end
            end
        end
        
        % Shows the available antenna ids
        function ids = antenna_ids(obj)
            ids = unique([obj.antenna_data.id]);
        end
        
        % Shows the available frequencies for (optionally specified) antenna type
        function freqs = frequencies(obj,varargin)
            if length(varargin)>0
                id = varargin{1};
                subset = obj.id_subset(id);
                freqs = unique([subset.antenna_data.frequency]);
            else
                freqs = unique([obj.antenna_data.frequency]);
            end
        end
        
        % For a given frequency and antenna id, return the available
        % frequency closest to it. freq in MHz
        function frequency = closest_frequency(obj,id,freq)
            id_subset = obj.id_subset(id);
            if isempty(id_subset)
                error('antenna type %d not found',id);
            end
            subset_freqs = id_subset.frequencies();
            diffs = abs(subset_freqs-freq);
            [C I] = min(diffs);
            frequency = subset_freqs(I);
        end
        
        % Get txapFileImporter object containing a subset with the chosen properties
        function obj = id_subset(obj,id)
            ids = [obj.antenna_data.id];
            idxs = (ids==id);
            obj.antenna_data = obj.antenna_data(idxs);
        end
        
        function obj = frequency_subset(obj,frequency)
            ids = [obj.antenna_data.frequency];
            idxs = (ids==frequency);
            obj.antenna_data = obj.antenna_data(idxs);
        end
        
        function obj = id_frequency_subset(obj,id,frequency,varargin) % accepts a "closest" parameter for the frequency
            if length(varargin)>0
                if strcmp(varargin{1},'closest')
                    closest_mode = true;
                else
                    error('option not recognized');
                end
            else
                closest_mode = false;
            end
            subset = obj.id_subset(id);
            
            if closest_mode
                the_frequency = obj.closest_frequency(id,frequency);
                if the_frequency~=frequency
                    warning('Closest frequency to %d MHz has been chosen: %d MHz',frequency,the_frequency);
                end
            else
                the_frequency = frequency;
            end
            obj = subset.frequency_subset(the_frequency);
        end
    end
end