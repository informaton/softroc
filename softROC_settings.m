function varargout = softROC_settings(varargin)
% SOFTROC_SETTINGS M-file for softroc_settings.fig
%
% softROC was written by Hyatt Moore IV at Stanford University
% questions or comments can be sent to him at sleepmoore@stanford.edu

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @softROC_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @softROC_settings_OutputFcn, ...
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


% --- Executes just before softroc_settings is made visible.
function softROC_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to softroc_settings (see VARARGIN)
%init_settings = varargin{1} = a struct containing the settings to initialize the GUI to
%if it is passed in.


settings = [];
%if the user has entered their own settings then
if(numel(varargin)==1) 
    settings = varargin{1};
else
    
    pfile = 'softROC_settings.plist';
    
    if(exist(pfile,'file'))
        %load it
        settings = plist.loadXMLPlist(pfile);
        
        %this can occur when exclude_string_value is left blank, and
        %loadXMLPlist does not return the field in this case
        if(~isfield(settings,'exclude_string_value'))
            settings.exclude_string_value = '';
        end
    end
end

if(~isempty(settings))
    
    
    %% handle the numeric exclusion criteria/settings first
    set(handles.radio_exclude_number,'value',settings.exclude_number_data);

    
    binops = get(handles.pop_binop_number,'string');
    binop_choice = find(strcmp(settings.binop_number,binops));
    
    if(~isempty(binop_choice))
        set(handles.pop_binop_number,'value',binop_choice);
    end

    if(isnumeric(settings.exclude_number_value))
        settings.exclude_number_value = num2str(settings.exclude_number_value);
    end
    
    set(handles.edit_exclude_number_value,'string',settings.exclude_number_value);

    if(settings.exclude_number_data)
        set(handles.edit_exclude_number_value,'enable','on');
        set(handles.pop_binop_number,'enable','on');
    else
        set(handles.edit_exclude_number_value,'enable','off');
        set(handles.pop_binop_number,'enable','off');
    end
    
    
    %% handle the string exclusion criteria/settings second

    set(handles.radio_exclude_string,'value',settings.exclude_string_data);

    binops = get(handles.pop_binop_string,'string');
    binop_choice = find(strcmp(settings.binop_string,binops));
    
    if(~isempty(binop_choice))
        set(handles.pop_binop_string,'value',binop_choice);
    end    
    
    set(handles.edit_exclude_string_value,'string',settings.exclude_string_value);

    if(settings.exclude_string_data)
        set(handles.edit_exclude_string_value,'enable','on');
        set(handles.pop_binop_string,'enable','on');
    else
        set(handles.edit_exclude_string_value,'enable','off');
        set(handles.pop_binop_string,'enable','off');
    end

    %% handle the interval method settings third
    set(handles.radio_steps_by_count,'value', settings.steps_by_count);
    set(handles.radio_steps_by_delta,'value', settings.steps_by_delta);
    
    
    %% handle the plot selection type fourth    
    set(handles.radio_plot_tpr_fpr,'value',settings.plot_tpr_vs_fpr);
    set(handles.radio_plot_se_sp,'value',settings.plot_se_vs_sp);

    %%handle the decimal places
    set(handles.menu_decimals,'value',settings.decimal_places+1);
    
    plot_options = {'All ROC results';
        'Optimal ROC results';
        'Training-Validation results'};
    
    plot_options_index = find(strcmp(plot_options,settings.plot_option));
    if(isempty(plot_options_index))
        plot_options_index = 1;
    end
    set(handles.menu_plot_options,'value',plot_options_index);
    
    training_validation_splits = [2/3,1/2]; %{'2/3-1/3','1/2-1/2'}
    training_validation_split_index = find(training_validation_splits==settings.training_validation_split);
    if(isempty(training_validation_split_index))
        training_validation_split_index=1;
    end
    set(handles.menu_trainingvalidationsplit,'value',training_validation_split_index);
    set(handles.edit_bootstrap_iterations,'string',num2str(settings.bootstrap_iterations));

    %update the gui based on the plot_option settings passed in
    menu_plot_options_Callback(handles.menu_plot_options,[],handles);
    
    set(handles.slider_se_vs_sp,'value',settings.se_vs_sp_alpha);
    slider_se_vs_sp_motion_callback(handles.slider_se_vs_sp,[]);  %update the text fields
    
    try
        addlistener(handles.slider_se_vs_sp,'ContinuousValueChange',@slider_se_vs_sp_motion_callback);
    catch me
       % fail silently; 
    end
       
    % Update handles structure
    guidata(hObject, handles);
    
    % UIWAIT makes softroc_settings wait for user response (see UIRESUME)
    uiwait(handles.softROC_settings_fig);
else
    errordlg({'Default settings not found!'},'Error loading');
    
    delete(hObject);
    
end



% --- Outputs from this function are returned to the command line.
function varargout = softROC_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%init_settings is a 

% Get default command line output from handles structure
if(ishandle(hObject))
    try
        output.exclude_number_data = get(handles.radio_exclude_number,'value');
        output.decimal_places = get(handles.menu_decimals,'value')-1;
%         binops = get(handles.pop_binop_number,'string');        

        binops = {'lt','le','==','gt','ge','~='};  %have to do it this way because the < cause the xmlplist function to fail since it does not like unequal numbers of < or > brackets
        binop_choice = get(handles.pop_binop_number,'value');
        output.binop_number =  binops{binop_choice};
        exclude_value = str2num(get(handles.edit_exclude_number_value,'string'));
        
        if(isempty(exclude_value))
            warndlg('A non numerical value was entered in the numerical exclusion criteria');
            exclude_value = get(handles.edit_exclude_number_value,'string');
        end
        output.exclude_number_value = exclude_value;
        
        output.exclude_string_data = get(handles.radio_exclude_string,'value');
        
        binops = get(handles.pop_binop_string,'string');
        binop_choice = get(handles.pop_binop_string,'value');
        output.binop_string =  binops{binop_choice};
        output.exclude_string_value = get(handles.edit_exclude_string_value,'string');
        
        
        output.steps_by_count = get(handles.radio_steps_by_count,'value');
        output.steps_by_delta = get(handles.radio_steps_by_delta,'value');
        
        output.plot_tpr_vs_fpr = get(handles.radio_plot_tpr_fpr,'value');
        output.plot_se_vs_sp = get(handles.radio_plot_se_sp,'value');
        
%         contents_cell = cellstr(get(Object,'string'));
%         output.plot_option = contents_cell{get(hObject,'value')};

        plot_options = {'All ROC results';
            'Optimal ROC results';
            'Training-Validation results'};
        
        output.plot_option = plot_options{get(handles.menu_plot_options,'value')};
        
        training_validation_splits = [2/3,1/2]; %{'2/3-1/3','1/2-1/2'}

        output.training_validation_split = training_validation_splits(get(handles.menu_trainingvalidationsplit,'value'));
        output.bootstrap_iterations = str2double(get(handles.edit_bootstrap_iterations,'string'));
        output.se_vs_sp_alpha = round(100*get(handles.slider_se_vs_sp,'value'))/100; %retrieve the percent value being shown to the user as a fraction
        varargout{1} = output;
    catch ME
        disp('Did not exit correctly');
        varargout{1} = [];
        
    end
    
    delete(handles.softROC_settings_fig);
else
    varargout{1} = [];
    
end

% --- Executes when selected object is changed in uipanel_numerical_exclusion.
function uipanel_numerical_exclusion_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_numerical_exclusion 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if(eventdata.NewValue==handles.radio_exclude_number)
    set(handles.pop_binop_number,'enable','on');
    set(handles.edit_exclude_number_value,'enable','on');
else
    set(handles.pop_binop_number,'enable','off');
    set(handles.edit_exclude_number_value,'enable','off');
end

guidata(hObject,handles);


% --- Executes when user attempts to close softROC_settings_fig.
function softROC_settings_fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to softROC_settings_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

uiresume(hObject);


% --- Executes when selected object is changed in uipanel_plot_settings.
function uipanel_plot_settings_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_settings 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)




% --- Executes when selected object is changed in uipanel_string_exclusion.
function uipanel_string_exclusion_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_string_exclusion 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if(eventdata.NewValue==handles.radio_exclude_string)
    set(handles.pop_binop_string,'enable','on');
    set(handles.edit_exclude_string_value,'enable','on');
else
    set(handles.pop_binop_string,'enable','off');
    set(handles.edit_exclude_string_value,'enable','off');
end

guidata(hObject,handles);



% --- Executes on selection change in menu_plot_options.
function menu_plot_options_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plot_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_plot_options contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_plot_options

contents_cell = cellstr(get(hObject,'string'));
option = contents_cell{get(hObject,'value')};

if(strcmp(option,'All ROC results')||strcmp(option,'Optimal ROC results'))
    set(handles.menu_trainingvalidationsplit,'enable','off');
    set(handles.edit_bootstrap_iterations,'enable','on');
elseif(strcmp(option,'Training-Validation results'));
    set(handles.menu_trainingvalidationsplit,'enable','on');
    set(handles.edit_bootstrap_iterations,'enable','off');
end

function edit_bootstrap_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bootstrap_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bootstrap_iterations as text
%        str2double(get(hObject,'String')) returns contents of edit_bootstrap_iterations as a double
value = str2double(get(hObject,'string'));
if(isempty(value) || value<1)
    value = 100;
end
set(hObject,'string',num2str(value));

guidata(hObject,handles);

function slider_se_vs_sp_motion_callback(hObject,eventdata)
%this function updates the Sensitivity and Speciificity percentage text
%fields based on the slider's position/value
    %se_pct_h is an alias for handles.text_se_pct_h
    %sp_pct_h is an alias for handles.text_sp_pct_h
    handles = guidata(hObject);    
    se_vs_sp_alpha = get(hObject,'value');
    set(handles.text_se_pct,'string',sprintf('SE (%i%%)',round((1-se_vs_sp_alpha)*100)));
    set(handles.text_sp_pct,'string',sprintf('SP (%i%%)',round((se_vs_sp_alpha)*100)));
    drawnow();
