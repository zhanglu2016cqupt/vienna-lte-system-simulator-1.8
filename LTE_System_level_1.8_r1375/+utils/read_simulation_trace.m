function [UE_traces,cell_mapping] = read_simulation_trace( LTE_config )
% Reads the pathloss (power) data from a trace file.
% (c) Josep Colom Ikuno, INTHFT, 2013
% www.nt.tuwien.ac.at

%% Read trace data
pathloss_filename = [LTE_config.trace_filename '_losses.txt'];
zipfile           = [LTE_config.trace_filename '.zip'];
trace_scale_s     = 0.1;

% Unzip trace file if needed
if ~exist(pathloss_filename,'file')
    if exist(zipfile,'file')
        fprintf('Unzipping %s\n',zipfile);
        [pathstr, ~, ~] = fileparts(zipfile);
        unzip(zipfile,pathstr);
    else
        error('Trace file "%s" of atrace file zip "%s" not found.',pathloss_filename,zipfile);
    end
end

% Read data
if LTE_config.debug_level>=1
    fprintf('Reading %s\n\n',pathloss_filename);
end
pathloss_data = dlmread(pathloss_filename,';');

% Extract basic info
idx_time           = 1;
idx_UEid           = 2;
idx_cellids        = 3:2:size(pathloss_data,2);
idx_cellpathlosses = 4:2:size(pathloss_data,2);

% Optionally delete first row
delete_first_sample = true;
if delete_first_sample
    non_first_positions = pathloss_data(:,idx_time) ~= 0;
    pathloss_data           = reshape(pathloss_data(non_first_positions(:,ones(1,size(pathloss_data,2)))),[],size(pathloss_data,2));
    
    % Adjust the time index column
    pathloss_data(:,1) = pathloss_data(:,1)-min(pathloss_data(:,1));
end

% Convert cell ID to one-indexed
pathloss_data(:,idx_cellids) = pathloss_data(:,idx_cellids);

% General info
time_samples = unique(pathloss_data(:,idx_time));
UE_ids       = unique(pathloss_data(:,idx_UEid));
N_cells      = length(idx_cellids);
t_length     = length(time_samples);
N_UEs        = length(UE_ids);

% All cell IDs
all_cell_ids_trace                           = unique(reshape(pathloss_data(:,idx_cellids),[],1));
all_cell_ids                                 = 1:length(all_cell_ids_trace);
traceCellId_to_simCellId(all_cell_ids_trace) = all_cell_ids;

cell_mapping.trace2sim = [(1:max(all_cell_ids_trace))',traceCellId_to_simCellId'];
cell_mapping.sim2trace = [all_cell_ids',all_cell_ids_trace];

UE.cellsIds          = NaN(t_length,N_cells);
UE.pathloss          = NaN(t_length,N_cells);
UE.unique_pathlosses = [];

% Translate from time to time Idx (change is TTI sampling changes)
time2timeIdx = @(s) round(s.*1/trace_scale_s+1);

UE_traces = UE(ones(1,N_UEs));
for u_=1:length(UE_ids)
    UE_rows       = pathloss_data(:,idx_UEid)==UE_ids(u_);
    
    % Rows where this UE is mentioned
    time_idxs     = time2timeIdx(pathloss_data(UE_rows,idx_time));
    UE_cell_idxs  = pathloss_data(UE_rows,idx_cellids);
    UE_pathloss   = pathloss_data(UE_rows,idx_cellpathlosses);
    NaN_rows      = sum(isnan(UE_pathloss),2)==N_cells;
    
    % Filter out NaN rows
    time_idxs    = time_idxs(~NaN_rows);
    UE_cell_idxs = UE_cell_idxs(~NaN_rows,:);
    UE_pathloss  = UE_pathloss(~NaN_rows,:);
    
    UE_traces(u_).trace_id              = UE_ids(u_);
    UE_traces(u_).time_idxs             = time_idxs;
    UE_traces(u_).cellsIds(time_idxs,:) = traceCellId_to_simCellId(UE_cell_idxs);
    UE_traces(u_).pathloss(time_idxs,:) = UE_pathloss;
    
    % Check which UEs do not move or have just NaNs
    if ~isempty(UE_pathloss)
        [~,I] = min(UE_pathloss,[],2);

        UE_cell_idxs_T = UE_cell_idxs';
        UE_traces(u_).attached_cell(1:max(time_idxs)) = 0;
        UE_traces(u_).attached_cell(time_idxs)        = traceCellId_to_simCellId(UE_cell_idxs_T(I'+N_cells*(0:(length(time_idxs)-1))));

        UE_traces(u_).UE_isOnlyNaNs = false;
        [unique_rows,ia,~]         = unique(UE_pathloss,'rows','R2012a');
        
        if size(unique_rows,1)==1
            UE_traces(u_).UE_didnt_move     = true;
            UE_traces(u_).nDifferentRows    = 1;
            UE_traces(u_).unique_pathlosses = unique_rows;
        else
            UE_traces(u_).UE_didnt_move         = false;
            UE_traces(u_).nDifferentRows        = size(unique_rows,1);
            UE_traces(u_).rowsWherePathlosJumps = sort(ia,'ascend');
            UE_traces(u_).unique_pathlosses     = UE_pathloss(UE_traces(u_).rowsWherePathlosJumps,:);
        end
    else
        UE_traces(u_).UE_isOnlyNaNs = true;
    end
end

% Delete NaN-only UEs from the trace
UE_traces = UE_traces(~[UE_traces.UE_isOnlyNaNs]);

% Ignores NaN-only users
UEs_with_changing_pathloss = ~[UE_traces.UE_didnt_move];
[C,~,~] = unique([UE_traces.nDifferentRows]);

% Fill in extra information in the first UE
UE_traces(1).all_cell_ids = all_cell_ids;

all_connected_cells          = unique([UE_traces.attached_cell]);
all_connected_cells          = all_connected_cells(all_connected_cells~=0);
UE_traces(1).connected_cells = all_connected_cells;

if LTE_config.debug_level>=1
    fprintf('%d/%d UEs had a changing pathloss over time.\n\n',sum(UEs_with_changing_pathloss),length(UE_traces));
    fprintf('Number of different rows (pathloss change):\n');
end
for i_=1:length(C)
    UE_subset   = [UE_traces.nDifferentRows]==C(i_);
    cardinality = sum(UE_subset);
    if C(i_)>1
        pathloss_change_positions = reshape([UE_traces(UE_subset).rowsWherePathlosJumps],1,[]);
        pathloss_change_positions = pathloss_change_positions(pathloss_change_positions~=1);
        [n,xout] = hist(pathloss_change_positions,5000);
        
        % Filter zero values
        filter_out = n==0;
        n    = round(n(~filter_out));
        xout = round(xout(~filter_out));
        
        if LTE_config.debug_level>=1
            fprintf('  -%d: %d/%d (%.2f%%), ',C(i_),cardinality,length(UE_traces),cardinality/length(UE_traces)*100);
            for j_=1:length(n)
                fprintf('%d->%.2f%%; ',xout(j_),n(j_)/sum(n)*100);
            end
            fprintf('\n');
        end
    else
        if LTE_config.debug_level>=1
            fprintf('  -%d: %d/%d (%.2f%%)\n',C(i_),cardinality,length(UE_traces),cardinality/length(UE_traces)*100);
        end
    end
end

end

