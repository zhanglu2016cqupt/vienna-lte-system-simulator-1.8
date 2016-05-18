function varargout = LTE_GUI_show_UEs_and_cells(varargin)
% LTE_GUI_SHOW_UES_AND_CELLS MATLAB code for LTE_GUI_show_UEs_and_cells.fig
%      LTE_GUI_SHOW_UES_AND_CELLS, by itself, creates a new LTE_GUI_SHOW_UES_AND_CELLS or raises the existing
%      singleton*.
%
%      H = LTE_GUI_SHOW_UES_AND_CELLS returns the handle to a new LTE_GUI_SHOW_UES_AND_CELLS or the handle to
%      the existing singleton*.
%
%      LTE_GUI_SHOW_UES_AND_CELLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LTE_GUI_SHOW_UES_AND_CELLS.M with the given input arguments.
%
%      LTE_GUI_SHOW_UES_AND_CELLS('Property','Value',...) creates a new LTE_GUI_SHOW_UES_AND_CELLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LTE_GUI_show_UEs_and_cells_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LTE_GUI_show_UEs_and_cells_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LTE_GUI_show_UEs_and_cells

% Last Modified by GUIDE v2.5 17-Jan-2012 19:00:31

% Plots eNodeB and UE positions as read from the laoded simulation results
%
% (c) Josep Colom Ikuno, INTHFT, 2012

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LTE_GUI_show_UEs_and_cells_OpeningFcn, ...
                   'gui_OutputFcn',  @LTE_GUI_show_UEs_and_cells_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LTE_GUI_show_UEs_and_cells is made visible.
function LTE_GUI_show_UEs_and_cells_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LTE_GUI_show_UEs_and_cells (see VARARGIN)

% Choose default command line output for LTE_GUI_show_UEs_and_cells
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Save user data (the simulation results)
simulation_data = varargin{1};

if length(varargin) > 1
    simulation_data.sim_results_GUI_handle     = varargin{2};
    set(handles.change_results_GUI_checkbox,'Value',true);
else
    simulation_data.sim_results_GUI_handle     = [];
    set(handles.change_results_GUI_checkbox,'Value',false);
end

set(hObject,'UserData',simulation_data);

UE_list = cell(1,length(simulation_data.UEs));
for u_=1:length(simulation_data.UEs)
    UE_list{u_} = sprintf('%g',u_);
end

cell_list = cell(1,length(simulation_data.sites));
for c_=1:length(simulation_data.eNodeBs)
    cell_list{c_} = sprintf('%g',c_);
end

simulation_data.UE_list   = UE_list;
simulation_data.cell_list = cell_list;

% Initialize checkboxes
set(handles.cell_show_cell_areas,'Value',true);
set(handles.cell_all,'Value',true);
set(handles.cell_show_ids,'Value',true);
set(handles.UE_all,'Value',true);

% Fill list boxes
set(handles.UE_listbox,'String',UE_list);
set(handles.UE_listbox,'Max',length(simulation_data.UEs));
set(handles.UE_listbox,'Min',0);

set(handles.cell_listbox,'String',cell_list);
set(handles.cell_listbox,'Max',length(simulation_data.eNodeBs));
set(handles.cell_listbox,'Min',0);

if isfield(simulation_data.LTE_config,'default_shown_GUI_cells') && ~isempty(simulation_data.LTE_config.default_shown_GUI_cells)
    cells_to_plot = simulation_data.LTE_config.default_shown_GUI_cells;
    if simulation_data.LTE_config.compact_results_file
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.the_eNodeB_traces)); % Filtero out possible out of range values
    else
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.simulation_traces.eNodeB_tx_traces)); % Filter out possible out of range values
    end
    set(handles.cell_listbox,'Value',cells_to_plot);
    if isfield(simulation_data,'simulation_traces')
        % non-compact results
        the_UE_traces = [simulation_data.simulation_traces.UE_traces];
    else
        the_UE_traces = simulation_data.the_UE_traces;
    end
    UE_in_this_cells = find(utils.resultsFileReader.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
    set(handles.UE_listbox,'Value',UE_in_this_cells);
    set(handles.cell_default,'Value',true);
else
    set(handles.cell_listbox,'Value',1:length(simulation_data.eNodeBs));
    set(handles.UE_listbox,'Value',1:length(simulation_data.UEs));
    set(handles.cell_default,'Value',false);
end

% Plot main plot
plot_cells_and_UEs(handles)

% UIWAIT makes LTE_GUI_show_UEs_and_cells wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = LTE_GUI_show_UEs_and_cells_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in UE_listbox.
function UE_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to UE_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UE_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UE_listbox
set(handles.UE_all,'Value',false);
set(handles.UE_none,'Value',false);
plot_cells_and_UEs(handles);

% --- Executes during object creation, after setting all properties.
function UE_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UE_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cell_listbox.
function cell_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to cell_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell_listbox
set(handles.cell_all,'Value',false);
set(handles.cell_none,'Value',false);
set(handles.cell_default,'Value',false);
simulation_data = get(handles.figure1,'UserData');
if isfield(simulation_data,'simulation_traces')
    % non-compact results
    the_UE_traces = [simulation_data.simulation_traces.UE_traces];
else
    the_UE_traces = simulation_data.the_UE_traces;
end
UE_in_this_cells = find(utils.resultsFileReader.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
set(handles.UE_listbox,'Value',UE_in_this_cells);
plot_cells_and_UEs(handles);

% --- Executes during object creation, after setting all properties.
function cell_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UE_show_ids.
function UE_show_ids_Callback(hObject, eventdata, handles)
% hObject    handle to UE_show_ids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UE_show_ids
plot_cells_and_UEs(handles);


% --- Executes on button press in cell_show_ids.
function cell_show_ids_Callback(hObject, eventdata, handles)
% hObject    handle to cell_show_ids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_show_ids
plot_cells_and_UEs(handles);


% --- Executes on button press in UE_all.
function UE_all_Callback(hObject, eventdata, handles)
% hObject    handle to UE_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UE_all
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.UE_none,'Value',false);
    set(handles.UE_listbox,'Value',1:length(simulation_data.UEs));
    plot_cells_and_UEs(handles);
end


% --- Executes on button press in UE_none.
function UE_none_Callback(hObject, eventdata, handles)
% hObject    handle to UE_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UE_none
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.UE_all,'Value',false);
    set(handles.UE_listbox,'Value',[]);
    plot_cells_and_UEs(handles);
end


% --- Executes on button press in cell_all.
function cell_all_Callback(hObject, eventdata, handles)
% hObject    handle to cell_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_all
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.cell_none,'Value',false);
    set(handles.cell_default,'Value',false);
    set(handles.cell_listbox,'Value',1:length(simulation_data.eNodeBs));
    if isfield(simulation_data,'simulation_traces')
        % non-compact results
        the_UE_traces = [simulation_data.simulation_traces.UE_traces];
    else
        the_UE_traces = simulation_data.the_UE_traces;
    end
    UE_in_this_cells = find(utils.resultsFileReader.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
    set(handles.UE_listbox,'Value',UE_in_this_cells);
    plot_cells_and_UEs(handles);
end


% --- Executes on button press in cell_none.
function cell_none_Callback(hObject, eventdata, handles)
% hObject    handle to cell_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_none
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.cell_all,'Value',false);
    set(handles.cell_default,'Value',false);
    set(handles.cell_listbox,'Value',[]);
    if isfield(simulation_data,'simulation_traces')
        % non-compact results
        the_UE_traces = [simulation_data.simulation_traces.UE_traces];
    else
        the_UE_traces = simulation_data.the_UE_traces;
    end
    UE_in_this_cells = find(utils.resultsFileReader.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
    set(handles.UE_listbox,'Value',UE_in_this_cells);
    plot_cells_and_UEs(handles);
end

function plot_cells_and_UEs(handles,varargin)

if isempty(varargin)
    plot_in_GUI = true;
else
    plot_in_GUI = varargin{1};
end

% Load data
simulation_data    = get(handles.figure1,'UserData');
eNodeB_sites       = simulation_data.sites;
eNodeBs            = simulation_data.eNodeBs;
UEs                = simulation_data.UEs;
networkPathlossMap = simulation_data.networkPathlossMap;

UEs_to_plot_lin    = get(handles.UE_listbox,'Value');
cells_to_plot      = get(handles.cell_listbox,'Value');
plot_cell_ids      = get(handles.cell_show_ids,'Value');
plot_cell_area     = get(handles.cell_show_cell_areas,'Value');
plot_UE_ids        = get(handles.UE_show_ids,'Value');
change_results_GUI = get(handles.change_results_GUI_checkbox,'Value');

if isfield(simulation_data,'simulation_traces')
    % non-compact results
    the_UE_traces = [simulation_data.simulation_traces.UE_traces];
else
    the_UE_traces = simulation_data.the_UE_traces;
end

if change_results_GUI && ~isempty(simulation_data.sim_results_GUI_handle) && ishandle(simulation_data.sim_results_GUI_handle)
    sim_results_GUI_handle_guidata = guidata(simulation_data.sim_results_GUI_handle);
    set(sim_results_GUI_handle_guidata.cell_listbox,'Value',cells_to_plot);
    set(sim_results_GUI_handle_guidata.cell_all,'Value',false);
    set(sim_results_GUI_handle_guidata.cell_none,'Value',false);
end

UE_pos             = reshape([UEs.pos],2,length(UEs))';
site_pos           = reshape([eNodeB_sites.pos],2,length(eNodeB_sites))';

if plot_in_GUI
    the_axes = handles.map_axes;
else
    the_new_figure = figure;
    the_axes       = axes('Parent',the_new_figure);
end

cla(the_axes);
hold(the_axes,'all');
set(the_axes,'YDir','normal');

% Shade the cell area of the specified cell/s
shaded_DotColor = 0.85*[1 1 1];
if plot_cell_area
    pos_to_shade = false(size(networkPathlossMap.sector_assignment));
    for c_ = cells_to_plot
        pos_to_shade = pos_to_shade | (networkPathlossMap.sector_assignment==c_);
    end
    [row, col] = find(pos_to_shade);
    shaded_points = LTE_common_pixel_to_pos([col row],networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
    scatter(the_axes,shaded_points(:,1),shaded_points(:,2),'Marker','.','MarkerFaceColor',shaded_DotColor,'MarkerEdgeColor',shaded_DotColor);
end

% Add the position of each UE
UE_DotColor_blue = [0 0 1];
dark_blue        = [0.7 0.7 1];
UE_DotColor_red  = [1 0 0];
dark_red         = [1 0.7 0.7];

UEs_to_plot = false(length(UEs),1);
UEs_to_plot(UEs_to_plot_lin) = true;
rest_UEs = true(length(UEs),1);
rest_UEs(UEs_to_plot) = false;
if simulation_data.LTE_config.compact_results_file
    completely_disabled_UEs = isnan([simulation_data.the_UE_traces.average_throughput_Mbps])';
else
    completely_disabled_UEs = isnan([simulation_data.simulation_traces.UE_traces.average_throughput_Mbps])';
end

if ~isempty(simulation_data.FFR_UE_mapping)
    % FFR simulation
    all_FR_UEs = simulation_data.FFR_UE_mapping.FR_assignment;
    all_PR_UEs = simulation_data.FFR_UE_mapping.PR_assignment;
    
    UEs_to_plot_FR_enabled  = all_FR_UEs & UEs_to_plot & ~completely_disabled_UEs;
    UEs_to_plot_PR_enabled  = all_PR_UEs & UEs_to_plot & ~completely_disabled_UEs;
    UEs_to_plot_FR_disabled = all_FR_UEs & UEs_to_plot &  completely_disabled_UEs;
    UEs_to_plot_PR_disabled = all_PR_UEs & UEs_to_plot &  completely_disabled_UEs;
    
    rest_UEs_FR_enabled      = all_FR_UEs & rest_UEs & ~completely_disabled_UEs;
    rest_UEs_PR_enabled      = all_PR_UEs & rest_UEs & ~completely_disabled_UEs;
    rest_UEs_FR_disabled     = all_FR_UEs & rest_UEs &  completely_disabled_UEs;
    rest_UEs_PR_disabled     = all_PR_UEs & rest_UEs &  completely_disabled_UEs;
    
    % unselected UEs (enabled)
    scatter(the_axes,UE_pos(rest_UEs_FR_enabled,1),UE_pos(rest_UEs_FR_enabled,2),'Marker','.','MarkerFaceColor',dark_blue,'MarkerEdgeColor',dark_blue);
    scatter(the_axes,UE_pos(rest_UEs_PR_enabled,1),UE_pos(rest_UEs_PR_enabled,2),'Marker','+','MarkerFaceColor',dark_blue,'MarkerEdgeColor',dark_blue);
    % unselected UEs (disaabled)
    scatter(the_axes,UE_pos(rest_UEs_FR_disabled,1),UE_pos(rest_UEs_FR_disabled,2),'Marker','.','MarkerFaceColor',dark_red,'MarkerEdgeColor',dark_red);
    scatter(the_axes,UE_pos(rest_UEs_PR_disabled,1),UE_pos(rest_UEs_PR_disabled,2),'Marker','+','MarkerFaceColor',dark_red,'MarkerEdgeColor',dark_red);
    
    % FR UEs
    scatter(the_axes,UE_pos(UEs_to_plot_FR_enabled,1),UE_pos(UEs_to_plot_FR_enabled,2),'Marker','.','MarkerFaceColor',UE_DotColor_blue,'MarkerEdgeColor',UE_DotColor_blue);
    scatter(the_axes,UE_pos(UEs_to_plot_FR_disabled,1),UE_pos(UEs_to_plot_FR_disabled,2),'Marker','.','MarkerFaceColor',UE_DotColor_red,'MarkerEdgeColor',UE_DotColor_red);
    
    % PR UEs
    scatter(the_axes,UE_pos(UEs_to_plot_PR_enabled,1),UE_pos(UEs_to_plot_PR_enabled,2),'Marker','+','MarkerFaceColor',UE_DotColor_blue,'MarkerEdgeColor',UE_DotColor_blue);
    scatter(the_axes,UE_pos(UEs_to_plot_PR_disabled,1),UE_pos(UEs_to_plot_PR_disabled,2),'Marker','+','MarkerFaceColor',UE_DotColor_red,'MarkerEdgeColor',UE_DotColor_red);
else
    % no FFR in use
    UEs_to_plot_enabled  = UEs_to_plot & ~completely_disabled_UEs;
    UEs_to_plot_disabled = UEs_to_plot &  completely_disabled_UEs;
    
    rest_enabled  = rest_UEs & ~completely_disabled_UEs;
    rest_disabled = rest_UEs &  completely_disabled_UEs;
    
    % unselected UEs (enabled)
    scatter(the_axes,UE_pos(rest_enabled,1),UE_pos(rest_enabled,2),'Marker','.','MarkerFaceColor',dark_blue,'MarkerEdgeColor',dark_blue);
    
    % unselected UEs (disabled)
    scatter(the_axes,UE_pos(rest_disabled,1),UE_pos(rest_disabled,2),'Marker','.','MarkerFaceColor',dark_blue,'MarkerEdgeColor',dark_red);
    
    % selected UEs (enabled)
    scatter(the_axes,UE_pos(UEs_to_plot_enabled,1),UE_pos(UEs_to_plot_enabled,2),'Marker','.','MarkerFaceColor',UE_DotColor_blue,'MarkerEdgeColor',UE_DotColor_blue);
    
    % selected UEs (disabled)
    scatter(the_axes,UE_pos(UEs_to_plot_disabled,1),UE_pos(UEs_to_plot_disabled,2),'Marker','.','MarkerFaceColor',UE_DotColor_blue,'MarkerEdgeColor',UE_DotColor_red);
end

if plot_UE_ids
    for u_=UEs_to_plot_lin
        text(UEs(u_).pos(1)+15*1,UEs(u_).pos(2),num2str(UEs(u_).id),'FontSize',8);
    end
end

% Plot a line that tells where the antennas are pointing
antenna_LineColor = 'blue';
vector_length     = 40;
for b_=1:length(eNodeB_sites)
    switch eNodeB_sites(b_).site_type
        case 'femto'
            % Do nothing
        otherwise
            % Plot the antenna/s
            origin = eNodeB_sites(b_).pos;
            for s_=1:length(eNodeB_sites(b_).sectors)
                angle = wrapTo360(-eNodeB_sites(b_).sectors(s_).azimuth+90);
                vector = vector_length*[ cosd(angle) sind(angle) ];
                destiny = vector + origin;
                plot([origin(1) destiny(1)],[origin(2) destiny(2)],antenna_LineColor);
            end
    end
end

% Plot all of the sites' positions
site_MarkerEdgeColor = 'black';
site_MarkerFaceColor = 'red';

the_xlims          = xlim;
total_xaxis        = diff(the_xlims);
small_hor_distance = 0.015;
small_distance     = small_hor_distance*total_xaxis;

for b_ = 1:length(eNodeB_sites)
    switch eNodeB_sites(b_).site_type
        case 'femto'
            % Triangular marker
            the_marker = 'diamond';
        otherwise
            % round marker
            the_marker = 'o';
    end
    scatter(eNodeB_sites(b_).pos(1),eNodeB_sites(b_).pos(2),'Marker',the_marker,'MarkerEdgeColor',site_MarkerEdgeColor,'MarkerFaceColor',site_MarkerFaceColor);
    
    % Add RRHs
    site_eNodeBs = [eNodeB_sites(b_).sectors];
    site_RRHs    = [site_eNodeBs.RRHs];
    
    if ~isempty(site_RRHs)
        graycolor = [0.5 0.5 0.5];
        RRH_pos = vertcat(site_RRHs.pos);
        for rrh_=1:length(site_RRHs)
            plot([eNodeB_sites(b_).pos(1) RRH_pos(rrh_,1)],[eNodeB_sites(b_).pos(2) RRH_pos(rrh_,2)],':','Color',graycolor);
        end
        scatter(RRH_pos(:,1),RRH_pos(:,2),'Marker',the_marker,'MarkerEdgeColor',graycolor,'MarkerFaceColor',graycolor);
    end
    
    switch eNodeB_sites(b_).site_type
        case 'femto'
            % No site number printed
        otherwise
            % Print site number
            text(eNodeB_sites(b_).pos(1)+small_distance,eNodeB_sites(b_).pos(2),num2str(b_),'Color','k');
    end
end

% Mark with a text the center of the cells
if plot_cell_ids
    for c_ = 1:length(eNodeBs) % cells_to_plot
        switch eNodeBs(c_).parent_eNodeB.site_type
            case 'femto'
                % Print cell id to the right of the cell
                text(eNodeBs(c_).parent_eNodeB.pos(1)+small_distance,eNodeBs(c_).parent_eNodeB.pos(2),num2str(c_),'Color',0.50*[1 1 1]);
            otherwise
                % Print cell id in the center of the cell
                text(networkPathlossMap.sector_centers(c_,1),networkPathlossMap.sector_centers(c_,2),num2str(c_),'HorizontalAlignment','center','Verticalalignment','middle','Color',0.50*[1 1 1]);
        end
        
    end
end

xlabel(the_axes,'x pos [m]');
ylabel(the_axes,'y pos [m]');
xlim(the_axes,networkPathlossMap.roi_x);
ylim(the_axes,networkPathlossMap.roi_y);
title(the_axes,sprintf('eNodeB and UE positions'));

if ~plot_in_GUI
    axis(the_axes,'equal');
end

% Release the hold on the axes
hold(the_axes,'off');

if isempty(UEs_to_plot_lin)
    mean_average_throughput   = NaN;
    mean_average_spectral_eff = NaN;
    mean_average_RBs_per_TTI  = NaN;
    RI_ratio                  = NaN;
else
    to_average                = [the_UE_traces(UEs_to_plot_lin).average_throughput_Mbps];
    to_average                = to_average(isfinite(to_average));
    mean_average_throughput   = mean(to_average);
    to_average                = [the_UE_traces(UEs_to_plot_lin).average_spectral_efficiency_bit_per_cu];
    to_average                = to_average(isfinite(to_average));
    mean_average_spectral_eff = mean(to_average);
    to_average                = [the_UE_traces(UEs_to_plot_lin).average_RBs_per_TTI];
    to_average                = to_average(isfinite(to_average));
    mean_average_RBs_per_TTI  = mean(to_average);
    
    RI_vect    = [the_UE_traces(UEs_to_plot_lin).RI];
    unique_RIs = unique(RI_vect);
    RI_hist    = hist(double(RI_vect),length(unique_RIs));
    RI_ratio   = RI_hist / sum(RI_hist);
end

sep = '------------------------------------';
statistics_text = sprintf([
    '%s\n',...
    'Simulations statistics:\n\n',...
    '%g UEs\n',...
    'Avg. UE throughput: %3.2f Mb/s\n',...
    'Avg. UE spectral eff.: %3.2f bit/cu\n',...
    'Avg. RBs/TTI/UE: %3.2f RBs\n',...
    'Rank Indicator (RI) distribution:',...
    ],...
    sep,...
    length(UEs_to_plot_lin),...
    mean_average_throughput,...
    mean_average_spectral_eff,...
    mean_average_RBs_per_TTI);
for r_=1:length(RI_ratio)
    statistics_text = sprintf('%s\n  rank %g: %3.2f%%',statistics_text,r_,RI_ratio(r_)*100);
end
statistics_text = sprintf('%s\n%s\n',statistics_text,sep);
if plot_in_GUI
    set(handles.UE_statistics_text,'String',statistics_text);
else
    fprintf('%s\n',statistics_text);
end


% --- Executes on button press in cell_show_cell_areas.
function cell_show_cell_areas_Callback(hObject, eventdata, handles)
% hObject    handle to cell_show_cell_areas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_show_cell_areas
plot_cells_and_UEs(handles);


% --- Executes on button press in open_plot_in_new_figure.
function open_plot_in_new_figure_Callback(hObject, eventdata, handles)
% hObject    handle to open_plot_in_new_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_cells_and_UEs(handles,false);


% --- Executes on button press in change_results_GUI_checkbox.
function change_results_GUI_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to change_results_GUI_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of change_results_GUI_checkbox
simulation_data    = get(handles.figure1,'UserData');
current_status     = get(hObject,'Value');
if ~isempty(simulation_data.sim_results_GUI_handle) && ishandle(simulation_data.sim_results_GUI_handle)
    if current_status
        cells_to_plot = get(handles.cell_listbox,'Value');
        sim_results_GUI_handle_guidata = guidata(simulation_data.sim_results_GUI_handle);
        set(sim_results_GUI_handle_guidata.cell_listbox,'Value',cells_to_plot);
        set(sim_results_GUI_handle_guidata.cell_all,'Value',false);
        set(sim_results_GUI_handle_guidata.cell_none,'Value',false);
    end
else
    set(handles.change_results_GUI_checkbox,'Value',~current_status);
end


% --- Executes on button press in cell_default.
function cell_default_Callback(hObject, eventdata, handles)
% hObject    handle to cell_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_default
simulation_data = get(handles.figure1,'UserData');

if isfield(simulation_data.LTE_config,'default_shown_GUI_cells') && ~isempty(simulation_data.LTE_config.default_shown_GUI_cells)
    if get(hObject,'Value')
        % UE_in_this_cells = UE_in_this_cells(UE_in_this_cells<=length(simulation_data.the_eNodeB_traces)); % Filtero out possible out of range values
        cells_to_plot = simulation_data.LTE_config.default_shown_GUI_cells;
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.the_eNodeB_traces)); % Filtero out possible out of range values
        set(handles.cell_listbox,'Value',cells_to_plot);
        if isfield(simulation_data,'simulation_traces')
            % non-compact results
            the_UE_traces = [simulation_data.simulation_traces.UE_traces];
        else
            the_UE_traces = simulation_data.the_UE_traces;
        end
        UE_in_this_cells = find(utils.resultsFileReader.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
        set(handles.UE_listbox,'Value',UE_in_this_cells);
        set(handles.cell_all,'Value',false);
        set(handles.cell_none,'Value',false);
        
        plot_cells_and_UEs(handles);
    end
else
    set(handles.cell_default,'Value',~get(hObject,'Value'));
end


% if get(hObject,'Value')
%     set(handles.cell_all,'Value',false);
%     set(handles.cell_default,'Value',false);
%     set(handles.cell_listbox,'Value',[]);
%     if isfield(simulation_data,'simulation_traces')
%         % non-compact results
%         the_UE_traces = [simulation_data.simulation_traces.UE_traces];
%     else
%         the_UE_traces = simulation_data.the_UE_traces;
%     end
%     UE_in_this_cells = find(utils.resultsFileReader.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
%     set(handles.UE_listbox,'Value',UE_in_this_cells);
%     plot_cells_and_UEs(handles);
% end
