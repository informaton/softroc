function varargout = softROC(varargin)
% SOFTROC M-file for softroc.fig
% 
% softROC was written by Hyatt Moore IV at Stanford University
% questions or comments can be sent to him at sleepmoore@stanford.edu
% Responses may even be sent back...


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @softROC_OpeningFcn, ...
                   'gui_OutputFcn',  @softROC_OutputFcn, ...
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


% --- Executes just before softroc is made visible.
function softROC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to softroc (see VARARGIN)

% Choose default command line output for softroc
handles.output = hObject;
handles.user.current_row = 1;

%remove the warning that gets thrown with 'basic' mode.  Basic mode is
%chosen in our case and is not a result of an undesired lack of
%performance.
handles.user.xlsread_warning_state = warning('off','MATLAB:xlsread:Mode');


reset_gui(hObject,handles);

% Update handles structure
% guidata(hObject, handles);

% UIWAIT makes softroc wait for user response (see UIRESUME)
% uiwait(handles.fig1);

function closeAllplots()
close all;

function removeRow(handles,row_index)
row_h = findobj(allchild(handles.pan_config),'-regexp','tag',['.*',num2str(row_index)]);
if(row_index > 1 && row_index <=3)
    set(row_h,'visible','off');
    %remove the problem of having a checkbox on, but then removing the file
    %and not seeing it or having the abilitiy to modify it. %softROC 1.6 fix
    for k=1:numel(row_h)
        if(strcmp(get(row_h(k),'style'),'checkbox'))
            set(row_h(k),'value',0);
        end
    end
    
elseif(row_index>3)
    config_pos = get(handles.pan_config,'position');
    config_pos(4) = config_pos(4)-handles.user.new_row_y_delta;
    set(handles.pan_config,'position',config_pos);
    
    config_children = allchild(handles.pan_config);
    
    config_positions = get(config_children,'position');
    
    config_positions = cell2mat(config_positions);
    % config_positions = mat2cell(config_positions,ones(numel(config_children),1));
    
    config_positions(:,2) = config_positions(:,2)-handles.user.new_row_y_delta;
    
    for k=1:numel(config_children)
        if(all(config_positions(k,:)>0))
            set(config_children(k),'position',config_positions(k,:));
        end
    end
    
    
    file_pos = get(handles.pan_file,'position');
    controls_pos =get(handles.pan_controls,'position');
    
    file_pos(2) = file_pos(2)-handles.user.new_row_y_delta;
    controls_pos(2) = controls_pos(2)-handles.user.new_row_y_delta;
    set(handles.pan_file,'position',file_pos);
    set(handles.pan_controls,'position',controls_pos);
    
    
    fig_pos = get(handles.fig1,'position');
    fig_pos(4) = fig_pos(4)-handles.user.new_row_y_delta;
    fig_pos(2) = fig_pos(2)+handles.user.new_row_y_delta;
    set(handles.fig1,'position',fig_pos);

    delete(row_h);
end

function handles = reset_gui(hObject,handles)

% closeAllplots();

        set(handles.check_start1,'enable','off');
        set(handles.check_end1,'enable','off');
        set(handles.push_reset,'enable','off');
        set(handles.push_delete_entry,'enable','off');
        set(handles.pop_control_label,'string','Label','value',1,'enable','off');
        set(handles.pop_control_value,'string','Value','value',1,'enable','off');

        
        set(handles.push_add_entry,'enable','off');
        set(handles.pop_binop1,'enable','off');
        set(handles.pop_key1,'string','Key 1','enable','off');
        set(handles.edit_start1,'string','1','enable','off');
        set(handles.edit_stop1,'string','1','enable','off');
        set(handles.edit_steps1,'string','1','enable','off');
        
        set(handles.pop_key2,'string','Key 2','enable','off');
        set(handles.edit_start2,'string','1','enable','off');
        set(handles.edit_stop2,'string','1','enable','off');
        set(handles.edit_steps2,'string','1','enable','off');
        
        set(handles.pop_key3,'string','Key 3','enable','off');
        set(handles.edit_start3,'string','1','enable','off');
        set(handles.edit_stop3,'string','1','enable','off');
        set(handles.edit_steps3,'string','1','enable','off');
        
        set(handles.push_save,'enable','off');
        set(handles.push_plot,'enable','off');
        
        set(handles.text_filename,'string','','enable','off');

for k=2:handles.user.current_row
    removeRow(handles,k);
end

handles.user.grouping_cell = {};
posRow1 = get(handles.pop_key1,'position');
posRow2 = get(handles.pop_logic2,'position');
handles.user.new_row_pos = posRow2;
handles.user.new_row_y_delta = posRow1(2)-posRow2(2); 
handles.user.current_row = 1;

pfile = 'softROC_settings.plist';
    
if(exist(pfile,'file'))
    %load it
    settings = plist.loadXMLPlist(pfile);
else
    %make it and save it for the future with these defaults
    settings.training_validation_split = 2/3; %split between training data and validation data
    settings.bootstrap_iterations = 100;
    settings.se_vs_sp_alpha = 0.5;
    
    plot_options = {'All ROC results';
        'Optimal ROC results';
        'Training-Validation results'};
    settings.plot_option = plot_options{1};
    
    settings.decimal_places = 2;
    
    settings.exclude_number_data = true;    
    settings.binop_number =  'lt';
    settings.exclude_number_value = 0;
    
    settings.exclude_string_data = false;    
    settings.binop_string =  '==';
    settings.exclude_string_value = '';
    
    settings.steps_by_count = false;
    settings.steps_by_delta = true;
    
    settings.plot_tpr_vs_fpr = false;
    settings.plot_se_vs_sp = true;

    settings.current_pathname = pwd;
    
    settings.pfile = pfile;
    plist.saveXMLPlist(pfile,settings);

end

    settings.log_filename = 'softROC_file_log.txt';
%this does not load correctly when the value is left blank
if(~isfield(settings,'exclude_string_value'))
    settings.exclude_string_value = '';
end
handles.user.settings = settings;

updateStepsText(handles);

guidata(hObject,handles);



function updateStepsText(handles)
%updates the gui to reflect the current settings used for step count/range
if(handles.user.settings.steps_by_count)
    set(handles.text_steps,'string',{'Steps by','Count'},'tooltipstring','This sets the total number|of steps to take from|start to stop values');
    
else
    set(handles.text_steps,'string',{'Steps by','Delta'},'tooltipstring','This sets the interval to take|between successive steps|starting at startst and continuing|while less than or equal to stop');
end
    


% --- Outputs from this function are returned to the command line.
function varargout = softROC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in push_select_file.
function push_select_file_Callback(hObject, eventdata, handles)
% hObject    handle to push_select_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

suggested_filename = '*.txt;*.xls;*.xlsx';

[filename,pathname]=uigetfile({suggested_filename,'Simple';...
    '*.*','All Files (*.*)'},'Data Loader',...
    handles.user.settings.current_pathname);
if(filename~=0)
     handles.user.settings.current_pathname = pathname;
     full_filename = fullfile(pathname,filename);
     try
         
        [number,txt,raw]=xlsread(full_filename,1,'','basic');
        
        settings = handles.user.settings;        

        %exclude bad rows on launch?     
        if(settings.exclude_number_data)
            numRows = size(raw,1)-1;
            numCols = size(raw,2);
            bad_rows =  false(numRows,1);
            for k=1:numCols                
                data = raw(2:end,k);
                if ~iscellstr(data)                     %#ok<ISCLSTR>
                    bad_rows =  bad_rows | feval(settings.binop_number,cell2mat(data),settings.exclude_number_value);
                end
            end
            
            bad_rows = [false;bad_rows];
            raw(bad_rows,:) = [];
        end
        
        handles.user.filename = filename;
        
        numCols = size(raw,2);
        unique_ind = false(1,numCols);
        unique_vals = cell(1,numCols);
        for k=1:numCols
            
            data = raw(2:end,k);
            ischarData = cellfun(@ischar, data);
            b = 0;
            if all(ischarData)
                b = unique(data);
            elseif ~any(ischarData)
                b = unique(cell2mat(data));
            else
                fprintf('Warning: Column %d contains a mix of numeric and character data and should not be used for testing.\n', k);
            end
            if(numel(b)==2)
                unique_ind(k)=true;
                unique_vals{k} = b;
            end            
        end
        
        handles.data.raw = raw(2:end,:);
        handles.data.labels = raw(1,:);

        handles.controls.indices = find(unique_ind);
        handles.controls.values = unique_vals(unique_ind);
        if(~iscell(handles.controls.values))
           handles.controls.values = {handles.controls.values}; 
        end
        
%         handles.controls.values = cell2mat(unique_vals(unique_ind));
%         handles.controls.raw = cell2mat(raw(2:end,handles.controls.indices));
        handles.controls.labels = handles.data.labels(handles.controls.indices); %leaves labels as a cell
        
        if(~iscell(handles.controls.labels))
            handles.controls.labels = {handles.controls.labels};
        end

        set(handles.check_start1,'enable','on');
        set(handles.check_end1,'enable','on');
        set(handles.push_reset,'enable','on');
        
        set([handles.check_start1,handles.check_start2,handles.check_start3],'callback',@check_start_Callback);
        set([handles.check_end1,handles.check_end2,handles.check_end3],'callback',@check_end_Callback);
        
        set(handles.pop_control_label,'string',handles.controls.labels,'enable','on');

        %populate the selection based on the current control (==1)
        pop_control_label_Callback(handles.pop_control_label,[],handles);
        
        set(handles.push_add_entry,'enable','on');
        set(handles.pop_binop1,'enable','on');
        set(handles.pop_key1,'string',handles.data.labels,'enable','on');
        set(handles.edit_start1,'string','1','enable','on');
        set(handles.edit_stop1,'string','1','enable','on');
        set(handles.edit_steps1,'string','1','enable','on');
        
        set(handles.pop_key2,'string',handles.data.labels,'enable','on');
        set(handles.edit_start2,'string','1','enable','on');
        set(handles.edit_stop2,'string','1','enable','on');
        set(handles.edit_steps2,'string','1','enable','on');
        
        set(handles.pop_key3,'string',handles.data.labels,'enable','on');
        set(handles.edit_start3,'string','1','enable','on');
        set(handles.edit_stop3,'string','1','enable','on');
        set(handles.edit_steps3,'string','1','enable','on');
        
        set(handles.push_save,'enable','on');
        set(handles.push_plot,'enable','on');
        
        
        set(handles.text_filename,'string',filename,'enable','inactive');
        
     catch ME
         % showME(ME);
         fprintf(newline);
         if isfield(handles, 'controls') && isempty(handles.controls.indices)
             warningMsg = sprintf('There was an error loading the file.\n\nAt least one column must contain two, and only two unique indices to represent the gold standard evaluation');             
         else
             warningMsg = sprintf('There was an error loading the file.\n\nVerify its integrity and Excel format. If the error involves ''biffparse'' try saving the XLS file using Microsoft Excel 95 format.');             
         end         
         warningMsg = sprintf('%s\n\nError: %s', warningMsg, ME.message);
         fprintf(1, '\n%s\n\n',warningMsg);
         warndlg(warningMsg, 'Warning','modal');
         set(handles.text_filename,'string','The file could not be loaded - please verify its format','enable','inactive');
     end
     plist.saveXMLPlist(handles.user.settings.pfile,handles.user.settings); %do this to save the most recent change to the pathname...
     
end

guidata(hObject,handles);



% --- Executes on selection change in pop_control_label.
function pop_control_label_Callback(hObject, eventdata, handles)
% hObject    handle to pop_control_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_control_label contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_control_label
value_options = handles.controls.values(:,get(hObject,'value'));

index = 1;

if(isnumeric(value_options))
    string = num2str(value_options);

    new_index = find(value_options==1,1);
    if(~isempty(new_index))
        index = new_index;
    end
else
    string = value_options;
    if(iscell(value_options))
       string = value_options{1};
       new_index = find(strcmp('1',value_options{1}));
       if(~isempty(new_index))
           index = new_index;
       end;
    end
end

set(handles.pop_control_value,'string',string,'value',index,'enable','on'); 



% --- Executes on button press in push_add_entry.
function push_add_entry_Callback(hObject, eventdata, handles)
% hObject    handle to push_add_entry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.user.current_row = handles.user.current_row+1;
set(handles.push_delete_entry,'enable','on');
if(handles.user.current_row<=3)
    next_row_h = findobj(allchild(handles.pan_config),'-regexp','tag',['.*',num2str(handles.user.current_row)]);

    set(next_row_h,'visible','on');
else
    
    last_row_h = findobj(allchild(handles.pan_config),'-regexp','tag',['.*',num2str(handles.user.current_row-1)]);
    props = get(last_row_h);

    config_pos = get(handles.pan_config,'position');
    config_pos(4) = config_pos(4)+handles.user.new_row_y_delta;
    set(handles.pan_config,'position',config_pos);
    
    config_children = allchild(handles.pan_config);
    
    config_positions = get(config_children,'position');
    
    config_positions = cell2mat(config_positions);
    % config_positions = mat2cell(config_positions,ones(numel(config_children),1));
    
    config_positions(:,2) = config_positions(:,2)+handles.user.new_row_y_delta;
    
    for k=1:numel(config_children)
        if(all(config_positions(k,:)>0))
            set(config_children(k),'position',config_positions(k,:));
        end
    end
    
    
    file_pos = get(handles.pan_file,'position');
    controls_pos =get(handles.pan_controls,'position');
    
    file_pos(2) = file_pos(2)+handles.user.new_row_y_delta;
    controls_pos(2) = controls_pos(2)+handles.user.new_row_y_delta;
    set(handles.pan_file,'position',file_pos);
    set(handles.pan_controls,'position',controls_pos);
    
    
    fig_pos = get(handles.fig1,'position');
    fig_pos(4) = fig_pos(4)+handles.user.new_row_y_delta;
    fig_pos(2) = fig_pos(2)-handles.user.new_row_y_delta;
    set(handles.fig1,'position',fig_pos);
    
    props = rmfield(props,'Selected');
    props = rmfield(props,'Extent');
    props = rmfield(props,'Type');

    pos = cell(numel(props),1);
    [pos{:}] = props.Position;
    props = rmfield(props,'Position');
    props = rmfield(props,'BeingDeleted');
    new_tag_suffix = num2str(handles.user.current_row);
    for k=1:numel(props)
       props(k).Position = pos{k};
       
       if(handles.user.current_row>10)
           props(k).Tag = [props(k).Tag(1:end-2),num2str(new_tag_suffix)]; %update the tag name
       else
           props(k).Tag = [props(k).Tag(1:end-1),num2str(new_tag_suffix)]; %update the tag name
       end
       

       
       handles.(props(k).Tag) = uicontrol(props(k));  %bang out the new row
    end

end

guidata(hObject,handles);


function controlStruct = getControlConfiguration(handles)
%returns structure of control information

controlStruct.index = get(handles.pop_control_label,'value');
controlStruct.label = handles.controls.labels{controlStruct.index};

value_index = get(handles.pop_control_value,'value');
values = get(handles.pop_control_value,'string');
if(iscell(values))
    diagnosis_value = str2num(values{value_index});
else
    diagnosis_value = str2num(values(value_index));
end
if(isempty(diagnosis_value))
    controlStruct.positive_diagnosis_value = values{value_index};
else
    controlStruct.positive_diagnosis_value = diagnosis_value;
end
if(ischar(handles.data.raw{1,handles.controls.indices(controlStruct.index)}))
    controlStruct.raw =strcmp(controlStruct.positive_diagnosis_value,handles.data.raw(:,handles.controls.indices(controlStruct.index)));
    controlStruct.positive_diagnosis_value = 1;
%     controlStruct.raw = char(handles.data.raw(:,handles.controls.indices(controlStruct.index)));
else
    controlStruct.raw = cell2mat(handles.data.raw(:,handles.controls.indices(controlStruct.index)));
end   


function configStruct = getParameterConfigurations(handles)
%this function extracts the parameters and builds a configuration for them
%to be returned as a configuration structure

numEntries = handles.user.current_row; %number of entries are defined by the user, who adds a row for each additional entry desired

num_configurations = 1; %the number of configurations grows as a cartesian product.  Make the first value be 1, so we can multiply to itself later with each additional configuration

for k = 1:numEntries
   binop = ['pop_binop',num2str(k)];
   binop_cell = get(handles.(binop),'string');
   params(k).binop = binop_cell{get(handles.(binop),'value')};
   
   key = ['pop_key',num2str(k)];
   key_cell = get(handles.(key),'string');
   params(k).data_column = get(handles.(key),'value'); %this index matches, the column of the header data
   params(k).label = key_cell{params(k).data_column}; % .. associated wih the label here
   
   
   if(k==1)
       params(k).logic = '   ';  %nothing to combine on the first entry/estimate made
   else
       logic = ['pop_logic',num2str(k)];
       logic_cell = get(handles.(logic),'string');
       params(k).logic = logic_cell{get(handles.(logic),'value')};
   
   end
   start = str2double(get(handles.(['edit_start',num2str(k)]),'string')); 
   stop = str2double(get(handles.(['edit_stop',num2str(k)]),'string')); 
   steps = str2double(get(handles.(['edit_steps',num2str(k)]),'string'));

   params(k).raw = cell2mat(handles.data.raw(:,params(k).data_column));
   
   if(handles.user.settings.steps_by_count)
       params(k).range = linspace(start,stop,steps)';
       num_steps = steps;
   else
       params(k).range = start:steps:stop;
       num_steps = numel(params(k).range);
   end
   
   num_configurations = num_configurations*num_steps;
   
   
end

config_cell = cell(numEntries,1);
if(numEntries == 1)
    config_cell{1} = params.range;
else
    [config_cell{:}] = ndgrid(params.range);
end


config_mat = zeros(num_configurations,numEntries);
for k = 1:numEntries
   config_mat(:,k) = config_cell{k}(:); 
end
configStruct.params = params;
configStruct.config_mat = config_mat;

% handles.user.params = params;
% handles.user.configuration_mat = config_mat;

% guidata(hObject,handles);


function okay = areGroupingsOkay(handles)
okay = true;
if(isempty(handles.user.grouping_cell))
    all_checkbox_h = findobj(allchild(handles.pan_config),'style','checkbox');
    values = get(all_checkbox_h,'value');
    if(iscell(values))
        values = cell2mat(values);
    end
    if(any(values))
        okay = false;
    end
end

% --- Executes on button press in push_plot.
function push_plot_Callback(hObject, eventdata, handles)
% hObject    handle to push_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.user.grouping_cell = getGroupings(handles);

% close all;
if(~areGroupingsOkay(handles))
    warndlg('Incorrect groupings found - will not continue');
else
    handles.user.paramStruct = getParameterConfigurations(handles);
    handles.user.controlStruct = getControlConfiguration(handles);
        
    %     plot_options = {'All ROC results';
    %         'Optimal ROC results';
    %         'Training-Validation results'};
    plot_option = handles.user.settings.plot_option;
    if(strcmp(plot_option,'Training-Validation results'))
        [confusion_training, sample_size_training, confusion_validation, sample_size_validation, config_mat_trainingvalidation] = getConfusionCrossValidation(handles);
        
        [Ttroc_fig, Tqroc_fig] = plotResults(handles,confusion_training,sample_size_training, config_mat_trainingvalidation);
        set(Ttroc_fig,'name','test Receiving Operating Characteristics (Training set)');        
        set(Tqroc_fig,'name','quality Receiving Operating Characteristics (Training set)');
        Tuserdata= get(Ttroc_fig,'userdata');
        title(Tuserdata.roc_axes_h,'Test Receiver Operating Characteristics (tROC) (Training set)')
        title(Tuserdata.qroc_axes_h,'Quality Receiver Operating Characteristics (qROC) (Training set)')
        set(Tuserdata.roc_h,'marker','o','markersize',4);
        set(Tuserdata.qroc_h,'marker','o','markersize',4);
     
        [Vtroc_fig, Vqroc_fig] = plotResults(handles,confusion_validation,sample_size_validation, config_mat_trainingvalidation);
        set(Vtroc_fig,'name','test Receiving Operating Characteristics (Validation set)');
        set(Vqroc_fig,'name','quality Receiving Operating Characteristics (Validation set)');
        Vuserdata= get(Vtroc_fig,'userdata');
        
        %update fields that should have combined handles when working with
        %cross-validation method
        combined_tROC_fig = [Ttroc_fig;Vtroc_fig];
        Vuserdata.tROC_fig = combined_tROC_fig;
        Tuserdata.tROC_fig = combined_tROC_fig;
        
        set(Ttroc_fig,'userdata',Tuserdata);        
        set(Vtroc_fig,'userdata',Vuserdata);        
        
        title(Vuserdata.roc_axes_h,'Test Receiver Operating Characteristics (tROC) (Validation set)')
        title(Vuserdata.qroc_axes_h,'Quality Receiver Operating Characteristics (qROC) (Validation set)')
        set(Vuserdata.roc_h,'marker','+','markersize',4);
        set(Vuserdata.qroc_h,'marker','+','markersize',4);
        
        handles_to_switch = [Tuserdata.handles_to_switch(:),Vuserdata.handles_to_switch(:)];
        
        %update the callbacks to use the combined_userdata 
        set(Tuserdata.roc_h,'buttondownfcn',{@ROC_Callback,Tuserdata}); 
        set(Vuserdata.roc_h,'buttondownfcn',{@ROC_Callback,Vuserdata}); 
        set(Tuserdata.qroc_h,'buttondownfcn',{@ROC_Callback,Tuserdata}); 
        set(Vuserdata.qroc_h,'buttondownfcn',{@ROC_Callback,Vuserdata}); 

        set(Tuserdata.roc_axes_h,'xlim',[0,1],'ylim',[0,1],'buttondownfcn',{@text_off_Callback,handles_to_switch});
        set(Tuserdata.qroc_axes_h,'xlim',[0,1],'ylim',[0,1],'buttondownfcn',{@text_off_Callback,handles_to_switch});
        set(Vuserdata.roc_axes_h,'xlim',[0,1],'ylim',[0,1],'buttondownfcn',{@text_off_Callback,handles_to_switch});
        set(Vuserdata.qroc_axes_h,'xlim',[0,1],'ylim',[0,1],'buttondownfcn',{@text_off_Callback,handles_to_switch});
    else
        [handles.user.confusion, sample_size] =  getConfusionAll(handles);
        if(strcmp(plot_option,'Optimal ROC results'))
            [TPR,FPR,~,~, ~,~,~, ~,~,~, ~] = confusion2roc(handles.user.confusion, sample_size);
            %now want to obtain indices of unique FPR with matching maximum TPR
            optimal_indices = getOutermostROCIndices(TPR,FPR);
            config_mat = handles.user.paramStruct.config_mat(optimal_indices,:);
            confusion = handles.user.confusion(optimal_indices,:);  %only get the best data points to validate-train with 
        elseif(strcmp(plot_option,'All ROC results'))
            confusion = handles.user.confusion;
            config_mat = handles.user.paramStruct.config_mat;
        end
        roc_fig=plotResults(handles,confusion,sample_size,config_mat);
        userdata = get(roc_fig,'userdata');

        %context menu for the figure
        uicontextmenu_handle = uicontextmenu('parent',roc_fig);
        uimenu(uicontextmenu_handle,'Label','Show optimal point','separator','off','callback',{@contextmenu_show_optimal_point_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Bootstrap optimal point','separator','off','callback',{@contextmenu_bootstrap_optimal_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Get unbiased ROC for optimal point','separator','off','callback',{@contextmenu_unbiased_optimal_point_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Show AUC for outermost curve','separator','off','callback',{@contextmenu_outer_AUC_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Bootstrap AUC 95% CI for outermost curve','separator','off','callback',{@contextmenu_bootstrap_outer_AUC_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Show outermost curve','separator','off','callback',{@contextmenu_show_outermost_curve_callback,handles,userdata});
        
        set(userdata.roc_axes_h,'uicontextmenu',uicontextmenu_handle);
        
        %context menu for the line
        uicontextmenu_handle = uicontextmenu('parent',roc_fig,'callback',{@contextmenu_parent,userdata});
        uimenu(uicontextmenu_handle,'Label','Bootstrap this point','separator','off','callback',{@contextmenu_bootstrap_this_point_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Get unbiased ROC for this point','separator','off','callback',{@contextmenu_unbiased_this_point_callback,handles,userdata});
        uimenu(uicontextmenu_handle,'Label','Bootstrap 95% ROC for this configuration','separator','off','callback',{@contextmenu_bootstrap_this_configuration_callback,handles,userdata});

        set(userdata.roc_h,'uicontextmenu',uicontextmenu_handle);
    end
end

guidata(hObject,handles);

function alpha = getAlphaFromPoint(index, TPR, FPR)
%obtain the sensitivity vs specificity alpha using supplied TPR and FPR
%values and the index of the point chosen

%obtain first convex hull of ROC values
selected_SE = TPR(index);
selected_SP = 1-FPR(index);
outer_indices = getOutermostROCIndices(TPR,FPR);
outerTPR = TPR(outer_indices);
outerFPR = FPR(outer_indices);


% model_fcn = (1-se_vs_sp_alpha)*TPR+(se_vs_sp_alpha).*(1-FPR); 
SE=outerTPR;
SP=1-outerFPR;

%sort the values in ascending order
[SP, s_ind] = sort(SP);
SE = SE(s_ind);

% %inflate our SE and SP to include 1,0 and 0,1 which is always possible
% if(SP(1)~=0 && SE(1)~=1)
%     SP = [0;SP(:)];
%     SE = [1;SE(:)];
% end
% if(SP(end)~=1 && SE(end)~=0)
%     SP = [SP(:);1];
%     SE = [SE(:);0];
% end

alpha_range = ones(numel(SE),2);   %last alpha is 1.0 - all weight toward SP
alpha_range(1) = 0;         %first alpha is 0.0 - meaning all weight toward SE

optimal_indices = false(size(SP));
optimal_indices(1,end) = true; %the first and last points are always in the optimal set

%obtain the slopes m=y(2)-y(1)/x(2)-x(1)
%for optimal_indices
for k=1:numel(optimal_indices)-1
   if(optimal_indices(k))
       slopes = (SE(k+1:end)-SE(k))./(SP(k+1:end)-SP(k));
       
       %this method is a bit more robust as it avoids the positive infinite
       %slope in cases where the added measures are already on the slope...
       [~, opt_index] = min(abs(slopes));
       
       % Here is the original way of doing it, knowing slopes are negative
       % [max_slope, max_index] = max(slopes);  %get the slope closest to zero

       opt_slope = slopes(opt_index);  %get the correct sign, which *should* be negative
       %solve for alpha using slope = -alpha/(1-alpha) ->  slope*(1-alpha) = -alpha
       % ->  alpha+slope*(1-alpha)=0  ->  alpha+slope-slope*alpha = 0 ->
       % alpha*(1-slope) = -slope  ->  alpha = -slope/(1-slope) = slope/(slope-1)
       alpha = opt_slope/(opt_slope-1);  %as slope increases, the alpha gets smaller which is what we expect/want
       alpha_range(k,2)=alpha;
       
       %update the next optimal points starting value
       alpha_range(k+opt_index,1)=alpha;
       
       optimal_indices(k+opt_index)=true;
       
   end
end

SE = SE(optimal_indices);
SP = SP(optimal_indices);
alpha_range = alpha_range(optimal_indices,:);

matched_ind = find(SE==selected_SE&SP==selected_SP);

%1. they selected a point that was considered optimal and has a range
%associated on each side of it. 
if(~isempty(matched_ind))
    alpha = alpha_range(matched_ind,2);

    %2. they selected a point in between the optimal points found, so we'll
%interpolate their alpha for them
else
    right_SP_ind = find(SP>selected_SP,1); %the right side of the range their selected point fell in
    if(isempty(right_SP_ind))
       right_SP_ind = numel(SP);  %give them the last/rightmost index available if there was a problem          
    end
    if(right_SP_ind ==1)
        left_SP = 0; 
    else
        left_SP = SP(right_SP_ind-1);
    end
    right_SP = SP(right_SP_ind);
    
    %if after all this, there is a problem, then set alpha to 0 and avoid
    %the divide by zero condition that lies ahead otherwise
    if(right_SP==left_SP)
        alpha = 0;
    else
        alpha = alpha_range(right_SP_ind,1)+diff(alpha_range(right_SP_ind,:))*(selected_SP-left_SP)/diff([right_SP,left_SP]);  %y=mx+b interpolation here
    end

end
function contextmenu_parent(hObject,eventdata,userdata)

set(gcf,'selectiontype','normal');
config_ind = ROC_Callback(userdata.roc_h,[],userdata);
childmenu_h = get(hObject,'children');
if(isempty(config_ind))
    set(childmenu_h,'visible','off')
else
    set(childmenu_h,'visible','on');
    set(childmenu_h,'userdata',config_ind);
end

function contextmenu_bootstrap_outer_AUC_callback(hObject,eventdata,handles,fig_userdata)
%calculates the area under the curve of the outermost curve....
getAUC_flag = true;
bootstrapFcn(hObject,eventdata,handles,[],getAUC_flag);

function contextmenu_show_outermost_curve_callback(hObject,eventdata,handles,fig_userdata)
%draw a curve from the given points

optimal_indices = getOutermostROCIndices(fig_userdata.TPR,fig_userdata.FPR);
outer_TPR = [0;fig_userdata.TPR(optimal_indices);1];
outer_FPR = [0;fig_userdata.FPR(optimal_indices);1];
userdata = get(gcf,'userdata');
set(userdata.roc_axes_h,'nextplot','add');
if(strcmpi(get(userdata.roc_h,'tag'),'se_vs_sp'))
    plot(userdata.roc_axes_h,1-outer_FPR,outer_TPR,'k:','hittest','off');
else
    plot(userdata.roc_axes_h,outer_FPR,outer_TPR,'k:','hittest','off');
end




function contextmenu_outer_AUC_callback(hObject,eventdata,handles,fig_userdata)
%calculates the area under the curve of the outermost curve....

optimal_indices = getOutermostROCIndices(fig_userdata.TPR,fig_userdata.FPR);
outer_TPR = [0;fig_userdata.TPR(optimal_indices);1];
outer_FPR = [0;fig_userdata.FPR(optimal_indices);1];
AUC = trapz(outer_FPR,outer_TPR);
decimal_format = ['%0.',num2str(handles.user.settings.decimal_places,'%d'),'f'];
message = ['Area Under Curve for outermost curve is ',num2str(AUC,decimal_format)];
msgbox(message,'Area Under Curve','modal');

function contextmenu_bootstrap_this_configuration_callback(hObject,eventdata,handles,fig_userdata)
%bootstrap the current points configuration (i.e. get(hObject,'userdata'))
%to obtain the 95% CI of sensitivity and specificity...
bootstrapFcn(hObject,eventdata,handles);

function contextmenu_bootstrap_this_point_callback(hObject,eventdata,handles,fig_userdata)
%find alpha based on mouse click selected point
config_ind = get(hObject,'userdata');
se_vs_sp_alpha = getAlphaFromPoint(config_ind,fig_userdata.TPR,fig_userdata.FPR);
bootstrapFcn(hObject,eventdata,handles,se_vs_sp_alpha);

function contextmenu_bootstrap_optimal_callback(hObject,eventdata,handles,fig_userdata)
%use preset-alpha from settings
showOptimalPoint(handles,fig_userdata);
se_vs_sp_alpha = handles.user.settings.se_vs_sp_alpha;  %relative weighting established by user
bootstrapFcn(hObject,eventdata,handles,se_vs_sp_alpha);

function contextmenu_unbiased_this_point_callback(hObject,eventdata,handles,fig_userdata)
%find alpha based on mouse click selected point
config_ind = get(hObject,'userdata');
se_vs_sp_alpha = getAlphaFromPoint(config_ind,fig_userdata.TPR,fig_userdata.FPR);
trainingvalidationFcn(hObject,eventdata,handles,se_vs_sp_alpha);

function contextmenu_unbiased_optimal_point_callback(hObject,eventdata,handles,fig_userdata)
%use preset-alpha from settings
showOptimalPoint(handles,fig_userdata);
se_vs_sp_alpha = handles.user.settings.se_vs_sp_alpha;  %relative weighting established by user
trainingvalidationFcn(hObject,eventdata,handles,se_vs_sp_alpha);

function trainingvalidationFcn(hObject,eventdata,handles,se_vs_sp_alpha)
%output the unbiased estimate obtained using a training-validation split
%on the weighting factor supplied (se_vs_sp_alpha)

set(gcf,'pointer','watch');
drawnow();

decimal_format = ['%0.',num2str(handles.user.settings.decimal_places,'%d'),'f'];

%badRows have already been excluded at this point
% handles = excludeBadRows(handles); %remove unwanted rows based on exclusion criteria established in settings

%obtain the training and validation data
[trainingControlStruct, validationControlStruct, trainingParamsStruct, validationParamsStruct] =getTrainingAndValidationSets(handles);

%evaluate the training set
config_mat = handles.user.paramStruct.config_mat;
[confusion_training,sample_size_training] = getConfusion(handles, trainingControlStruct, trainingParamsStruct, config_mat);
[TPR,FPR,~,~, ~,~,~, ~,~,~, ~] = confusion2roc(confusion_training, sample_size_training);


model_fcn = (1-se_vs_sp_alpha)*TPR+(se_vs_sp_alpha).*(1-FPR);    
optimal_ind = max(model_fcn)==model_fcn;
% optimal_config_mat = config_mat(optimal_ind,:);

optimal_config_mat = config_mat(optimal_ind,:);
[confusion_validation,sample_size_validation] = getConfusion(handles, validationControlStruct, validationParamsStruct, optimal_config_mat);
[unbiasedTPR,unbiasedFPR,~,~, ~,~,~, ~,~,~, ~] = confusion2roc(confusion_validation, sample_size_validation);

%sometimes we'll have repetitions or multiple entries due to the confusion
%matrix -> writing unbiasedTPR = unbiasedTPR(1) would be an 'equivalent'
%solution.
unbiasedTPR = mean(unbiasedTPR);
unbiasedFPR = mean(unbiasedFPR);

% handles.user.settings.training_validation_split;
message = cell(4,1);
message{2} = sprintf('Training sample size = %i samples\nValidation sample size = %i',sample_size_training,sample_size_validation);
message{1} = sprintf('Unbiased ROC estimate using training-validation split based on');

if(handles.user.settings.plot_tpr_vs_fpr)    
    message{3} = ['TPR = ',num2str(unbiasedTPR,decimal_format)];
    message{4} = ['FPR = ',num2str(unbiasedFPR,decimal_format)];
else
    message{3} = ['SE = ',num2str(unbiasedTPR,decimal_format)];
    message{4} = ['SP = ',num2str(1-unbiasedFPR,decimal_format)];
end

set(gcf,'pointer','arrow');
drawnow();

msgbox(message,'modal');
% for k = 1:numel(all_userdata.param_labels)
%     message{k+1} = [all_userdata.param_labels{k},' [ ',num2str(config_interval(1,k),decimal_format),', ',num2str(config_interval(2,k),decimal_format),']'];
% end


function bootstrapFcn(hObject,eventdata,handles,se_vs_sp_alpha,getAUC_flag)
%hObject's userdata will contain the config_ind of the selected point.
%config_ind is the index of the point selected by the user which indexes
%into the TPR and FPR field values in all_userdata which can be used to
%create a model for the rule/weight of this point in choosing the matching
%points from the bootstrap runs.
%
%If se_vs_sp_alpha is empty, then a bootstrap of the current configuration
%is performed in order to obtain 95% CI of sensitivity and specficity of
%the current configuration index instead.
%
%if getAUC is included and true then compute the 95% AUC CI as well

%retrieve the userdata from the parent figure....
all_userdata = get(get(get(hObject,'parent'),'parent'),'userdata');

decimal_format = ['%0.',num2str(handles.user.settings.decimal_places,'%d'),'f'];


% 1. get all bootstrap indices for the number of bootstraps desired
% (settings)
sample_size = all_userdata.sample_size;
numBootstraps = handles.user.settings.bootstrap_iterations;
boot_ind = randi(sample_size,[numBootstraps,sample_size]);
 
% 2. run the bootstrap function using these samples instead
handles = excludeBadRows(handles); %remove unwanted rows based on exclusion criteria established in settings

params = handles.user.paramStruct.params;
controls = handles.user.controlStruct;
config_mat = handles.user.paramStruct.config_mat;

if(nargin<4)
    %determine ROC confidence interval for one configuration
    config_ind = get(hObject,'userdata');
    se_vs_sp_alpha = [];
    config_mat = config_mat(config_ind,:);  %only use the selected configuration...
    bootstrapped_TPR = zeros(numBootstraps,1);
    bootstrapped_FPR = zeros(numBootstraps,1); 
    getAUC_flag = false;
else
    %otherwise determine the configuration Confidence Intervals
    bootstrapped_config = zeros(numBootstraps,size(config_mat,2));
    
    %if getAUC_flag is included and true, compute the 95% AUC CI as well
    if(nargin<5 || isempty(getAUC_flag))
        getAUC_flag = false;
    elseif(getAUC_flag)
        bootstrapped_AUC = zeros(numBootstraps,1);  
    end
end;

%could estimate the time required for one iteration ... 
%[TPR,FPR,K_1_0,K_0_0, CohensKappa,PPV,NPV, EFF,chi_square,P, Q] = confusion2roc(handles.user.confusion, sample_size);
 

h=waitbar(0,'Bootstrapping');


num_symptoms = numel(params);
for n = 1:numBootstraps
    
    %resample the raw symptoms with matching resample of
    %controls/positivity (i.e. controls.raw)
    for k=1:num_symptoms
        params(k).raw = handles.user.paramStruct.params(k).raw(boot_ind(n,:),:);
    end    
    controls.raw = handles.user.controlStruct.raw(boot_ind(n,:),:);
    
    %get the confusion values for the resampled data and obtain ROC results
    [confusion,sample_size] = getConfusion(handles, controls, params, config_mat);
    
    [TPR,FPR,K_1_0,K_0_0, CohensKappa,PPV,NPV, EFF,chi_square,P, Q] = confusion2roc(confusion, sample_size);

    if(getAUC_flag)
        optimal_indices = getOutermostROCIndices(TPR,FPR);
        outer_TPR = [0;TPR(optimal_indices);1];
        outer_FPR = [0;FPR(optimal_indices);1];
        AUC = trapz(outer_FPR,outer_TPR);
        bootstrapped_AUC(n) = AUC;
    else
        if(~isempty(se_vs_sp_alpha))
            %get the optimal point based on the rule established via slider or
            %model parameter extraction
            model_fcn = (1-se_vs_sp_alpha)*TPR+(se_vs_sp_alpha).*(1-FPR);
            bootstrap_config_ind = max(model_fcn)==model_fcn;
            bootstrapped_config(n,:) = mean(config_mat(bootstrap_config_ind,:),1); %take the mean of the columns
        else
            bootstrapped_TPR(n) = TPR(1); %put the (1) here in case duplicate entries returned above
            bootstrapped_FPR(n) = FPR(1);
        end
    end
    
    
    if(mod(n,5)==0)
        waitbar(n/numBootstraps,h);
    end
end

%calculate 95% confidence intervals using (1) basic percentile and (2)
%normal/gaussian approximations

%1 basic percentile
CI = 95;  %i.e. we want the 95% confidence interval
CI_alpha = 100-CI;
CI_range = [CI_alpha/2, 100-CI_alpha/2];  %[2.5, 97.5]

if(getAUC_flag)
    ROC_label = 'AUC';
    AUC_CI_percentile = prctile(bootstrapped_AUC,CI_range);
    message = cell(2,1); %header+AUC
    message{1} = sprintf('95%% Confidence Interval of ROC''s AUC using\n\tbootstrap iterations = %i\n\tsample size = %i\n',numBootstraps,sample_size);
    message{2} = sprintf(['%s [',decimal_format,', ',decimal_format,']'],ROC_label,AUC_CI_percentile(1)*100,AUC_CI_percentile(2)*100);
else
    if(isempty(se_vs_sp_alpha))
        ROC_CI_percentile = prctile([bootstrapped_TPR,bootstrapped_FPR],CI_range); %==> [SE_lower,SP_lower; SE_higher, SP_higher]
        message = cell(num_symptoms+3,1); %symptoms+header+TPR+FPR
        message{1} = sprintf('95%% Confidence Interval(s) using\n\tbootstrap iterations = %i\n\tsample size = %i\n',numBootstraps,sample_size);
        
        if(handles.user.settings.plot_tpr_vs_fpr)
            ROC_label1 = 'TPR';
            ROC_label2 = 'FPR';
        else
            ROC_label1 = 'Sensitivity';
            ROC_label2 = 'Specificity';
            ROC_CI_percentile(:,2) = 1-flipud(ROC_CI_percentile(:,2)); %need to flip this so that lowever value is on top after the subtraction
        end
        paramCell = handles.user.paramStruct.params;
        for k = 1:numel(paramCell)
            message{k+1} = [paramCell(k).label,'(',num2str(config_mat(k),decimal_format),')'];
        end
        message{end-1} = sprintf(['%s [',decimal_format,', ',decimal_format,']'],ROC_label1,ROC_CI_percentile(1,1)*100,ROC_CI_percentile(2,1)*100);
        message{end} = sprintf(['%s [',decimal_format,', ',decimal_format,']'],ROC_label2,ROC_CI_percentile(1,2)*100,ROC_CI_percentile(2,2)*100);
    else
        config_CI_percentile = prctile(bootstrapped_config,CI_range);
        %handle the case when only one row/configuration is made in which case the
        %prctile function outputs a (1x2) vector that is incompatible with the
        %message labeling system used below
        if(numel(config_CI_percentile)==2)
            config_CI_percentile = config_CI_percentile';
        end
        message = cell(num_symptoms+1,1);
        message{1} = sprintf('95%% Confidence Interval(s) using\n\tbootstrap iterations = %i\n\tsample size = %i\n',numBootstraps,sample_size);
        
        paramCell = handles.user.paramStruct.params;
        for k = 1:numel(paramCell)
            message{k+1} = [paramCell(k).label,' [ ',num2str(config_CI_percentile(1,k),decimal_format),', ',num2str(config_CI_percentile(2,k),decimal_format),']'];
        end
        
    end
end
delete(h);
msgbox(message,'Confidence','modal');
    

%in retrospect - gaussian based confidence intervals were removed as the
%configurations are uniformly distributed across user established ranges

%2 normal/gaussian approxmation
% observed_config = config_mat(handles.user.config_mat(config_ind));
%get the optimal point based on the rule established via slider or
%model parameter extraction

% model_fcn = (1-se_vs_sp_alpha)*all_userdata.TPR+(se_vs_sp_alpha).*(1-all_userdata.FPR);
% original_config_ind = max(model_fcn)==model_fcn;
% original_config = mean(config_mat(original_config_ind,:),1); %take the mean of the columns for the original dataset...
% 
% stdev = std(bootstrapped_config,0,1);  %estimate standard deviation down the columns
% % bias = mean(bootstrapped_config-repmat(original_config,size(bootstrapped_config,1),1),1);
% m = mean(bootstrapped_config,1); %similarly estimate the mean
% bias = m-
% Z = norminv(CI_range/100);  %get the z-values for our tails
% config_CI_normal = [m+stdev*Z(1); m+stdev*Z(2)];



% if(numBootstraps>=100)
%     confid_int = round(numBootstraps*[.025, 0.975]);
%     config_sorted = sort(bootstrapped_config,1); %sort down columns
%     config_interval = config_sorted(confid_int,:);
% end


function contextmenu_show_optimal_point_callback(hObject,~,handles,userdata)
showOptimalPoint(handles,userdata);

function showOptimalPoint(handles,userdata)
%selects the optimal point as determined by user setting of se_vs_sp_alpha
se_vs_sp_alpha = handles.user.settings.se_vs_sp_alpha;  %relative weighting established by user

model_fcn = (1-se_vs_sp_alpha)*userdata.TPR+(se_vs_sp_alpha).*(1-userdata.FPR);
optimal_ind = find(max(model_fcn)==model_fcn,1);
% userdata.tROC_fig;
ROC_Callback(userdata.roc_h,[],userdata,optimal_ind);


function makeTitle(handles)

params = handles.user.paramStruct.params;
group_cell = handles.user.grouping_cell;
group_mat = zeros(numel(params),2);
for k = 1:numel(group_cell)
    group_mat(group_cell{k}(1),1) = 1;
    group_mat(group_cell{k}(2),2) = 1;
end
titleStr = '';
for k=1:numel(params)

    group_label = '  ';
    if(group_mat(k,1))
        group_label(1) = '(';
    end
    if(group_mat(k,2))
        group_label(2) = ')';
    end
    
    titleStr = sprintf('%s%s %c%s%c ',titleStr,params(k).logic,group_label(1),params(k).label,group_label(2));
    
end
title(titleStr,'fontsize',12);

function varargout = plotResults(handles,confusion, sample_size, optional_config_mat)
%plotResults(handles,confusion, sample_size, optional_config_mat)
% generates two interactive figures from the supplied confusion matrix and
% sample_size variables
% if optional_config_mat is not included then handles.user.config_mat is
%  used
%
%[tROC_fig] = plotResults(handles,confusion, sample_size, optional_config_mat)
%[tROC_fig, qROC_fig] = plotResults(handles,confusion, sample_size, optional_config_mat)
% tROC_fig is the handle to the figure with tROC plot
% qROC_fig is the handle to the figure with the qROC plot

    %%plot the data here
    if(nargin>3 && isempty(optional_config_mat))
        config_mat = optional_config_mat;
    else
        config_mat = handles.user.paramStruct.config_mat;
    end
    
    [TPR,FPR,K_1_0,K_0_0, CohensKappa,PPV,NPV, EFF,chi_square,P, Q] = confusion2roc(confusion, sample_size);
    tROC_fig = figure;
    set(tROC_fig,'name','test Receiving Operating Characteristics','visible','off');

    if(handles.user.settings.plot_tpr_vs_fpr)
        roc_h = plot(FPR,TPR,'b.','markersize',10);
        xlabel('False Positive Rate','fontsize',12);
        ylabel('True Positive Rate','fontsize',12);
        tag = 'tpr_vs_fpr';
    else
        roc_h= plot(1-FPR,TPR,'b.','markersize',10);
        xlabel('Specificity','fontsize',12);
        ylabel('Sensitivity','fontsize',12);
        tag = 'se_vs_sp';
    end
    roc_highlight_h = line('parent',gca,'linestyle','none','hittest','off','marker','o',...
        'visible','off','color','r','markersize',12,'tag','highlighter','userdata',tag); %keep the tag to differentiate between the two plot configurations for the tROC
%     title('Test Receiver Operating Characteristics (tROC)','fontsize',12);
    makeTitle(handles);
    %     roc_text_h = text('parent',gca,'visible','off','edgecolor','k','hittest','off','position',[0.05,0.15,0]);
    roc_text_h = annotation(gcf,'textbox','visible','off','edgecolor','k','hittest','off','fitboxtotext','on','fontsize',12); %,position',[0.175,0.25,0.1,0.1]);
    roc_axes_h = gca;
    
    qROC_fig = figure;
    set(gcf,'name','quality Receiver Operating Characteristics');
    
    qroc_h= plot(K_0_0,K_1_0,'r.','markersize',10);
    %     qroc_text_h = text('parent',gca,'visible','off','edgecolor','k','hittest','off','position',[0.05,0.15,0]);
    qroc_text_h = annotation(gcf,'textbox','visible','off','edgecolor','k','hittest','off','fitboxtotext','on','fontsize',12);%'position',[0.175,0.25,0.1,0.1]);
    qroc_axes_h = gca;
    qroc_highlight_h = line('parent',gca,'linestyle','none','hittest','off','marker','o',...
        'visible','off','color','g','markersize',12,'tag','highlighter');
    
    xlabel('\kappa (0,0)');
    ylabel('\kappa (1,0)');
    title('Quality Receiver Operating Characteristics (qROC)','fontsize',12);
    
    userdata.qroc_text_h = qroc_text_h;
    userdata.roc_text_h = roc_text_h;
    userdata.roc_highlight_h = roc_highlight_h;
    userdata.qroc_highlight_h = qroc_highlight_h;
    userdata.param_labels = cell(numel(handles.user.paramStruct.params),1);
    [userdata.param_labels{:}] = handles.user.paramStruct.params.label;
    userdata.config_mat = config_mat;
    userdata.TPR = TPR;
    userdata.FPR = FPR;
    userdata.K_1_0 = K_1_0;
    userdata.K_0_0 = K_0_0;
    userdata.CohensKappa = CohensKappa;
    userdata.PPV = PPV;
    userdata.NPV = NPV;
    userdata.EFF = EFF;
    userdata.chi_square = chi_square;
    userdata.Q = Q;
    userdata.P = P;
    userdata.sample_size = sample_size;
    userdata.decimal_places = handles.user.settings.decimal_places;
    userdata.tROC_fig = tROC_fig;
    
    %apply userdata to both functions so that I can differentiate between
    %the two
    set(roc_h,'buttondownfcn',{@ROC_Callback,userdata},'tag',tag); %use tag here to differentiate between the two tROC views (i.e. sp or fpr)
    set(qroc_h,'buttondownfcn',{@ROC_Callback,userdata},'tag','qROC');
    handles_to_switch = [roc_text_h,qroc_text_h,qroc_highlight_h,roc_highlight_h];
    set(roc_axes_h,'xlim',[0,1],'ylim',[0,1],'buttondownfcn',{@text_off_Callback,handles_to_switch});
    set(qroc_axes_h,'xlim',[0,1],'ylim',[0,1],'buttondownfcn',{@text_off_Callback,handles_to_switch});
    
    %these are added to make it easier to manipulate the titles and such
    %when dealing with cross-validation techniques
    userdata.roc_h = roc_h;
    userdata.qroc_h = qroc_h;
    userdata.roc_axes_h = roc_axes_h;
    userdata.qroc_axes_h = qroc_axes_h;
    userdata.handles_to_switch = handles_to_switch;

    %make it visible so that it is also shown on top
    set(tROC_fig,'userdata',userdata,'visible','on');

    if(nargout>0)
        varargout{1} = tROC_fig;
    end
    if(nargout>1)
        varargout{2} = qROC_fig;
    end

function text_off_Callback(hObject,eventdata,handles2hide)
if(~strcmp(get(gcf,'selectiontype'),'alt'))
    set(handles2hide,'visible','off');
end

function [trainingControlStruct, validationControlStruct, trainingParamsStruct, validationParamsStruct] = getTrainingAndValidationSets(handles)
%produce training and validation sets based on current control structure
%and validation cutoff selected.  The positive and negative diagnostic
%cases are separated and then each is broken into a training and validation
%group corresponding ot the percentange cutoff established by the user.
%The default is 2/3 to be allocated to the training set and 1/3 allocated
%to the validation set.

params = handles.user.paramStruct.params;
controls = handles.user.controlStruct;

diagnosis_ind = (controls.raw == controls.positive_diagnosis_value); 

positive_diagnosis_ind = find(diagnosis_ind);
positive_sample_size = numel(positive_diagnosis_ind);
[~,rand_positive_diagnosis_ind] = sort(rand(positive_sample_size,1));
rand_positive_diagnosis_ind = positive_diagnosis_ind(rand_positive_diagnosis_ind);

negative_diagnosis_ind = find(~diagnosis_ind);
negative_sample_size = numel(negative_diagnosis_ind);
[~,rand_negative_diagnosis_ind] = sort(rand(negative_sample_size,1));
rand_negative_diagnosis_ind = negative_diagnosis_ind(rand_negative_diagnosis_ind);

%split into two groups of two
%first group of two is positive values, the first is for the training the
%second is for validation
numPosTraining = round(positive_sample_size*handles.user.settings.training_validation_split);
numPosValidation = positive_sample_size-numPosTraining;

%second group is negative values, the first is for the training and the
%second is for the validation
numNegTraining = round(negative_sample_size*handles.user.settings.training_validation_split);
numNegValidation = negative_sample_size-numNegTraining;

%create two structures with same fields as handles.user.controlStruct but
%with values separated based on inclusion in either the training or
%validation set
trainingControlStruct.raw = [repmat(controls.positive_diagnosis_value,numPosTraining,1);repmat(~controls.positive_diagnosis_value,numNegTraining,1)];
trainingControlStruct.positive_diagnosis_value = controls.positive_diagnosis_value;

validationControlStruct.raw = [repmat(controls.positive_diagnosis_value,numPosValidation,1);repmat(~controls.positive_diagnosis_value,numNegValidation,1)];
validationControlStruct.positive_diagnosis_value = controls.positive_diagnosis_value;

trainingParamsStruct = params;
validationParamsStruct = params;

%grab the corresponding parameter values here.
for entry=1:numel(params) 
   trainingParamsStruct(entry).raw = [params(entry).raw(rand_positive_diagnosis_ind(1:numPosTraining));params(entry).raw(rand_negative_diagnosis_ind(1:numNegTraining))];
   validationParamsStruct(entry).raw = [params(entry).raw(rand_positive_diagnosis_ind(numPosTraining+1:end));params(entry).raw(rand_negative_diagnosis_ind(numNegTraining+1:end))];
end

function optimal_indices = getOutermostROCIndices(TPR,FPR)
%retrieves the optimal ROC indices based on the following algorithm...
%optimality here is the outer ROC curve when traversing the x-axis values
%(i.e. FPR or (1-SP))
u_FPR = unique(FPR);
final_index = zeros(numel(u_FPR),1);

%obtain unique specificity values, and then search for maximum sensitivity
%where duplicates may exist.  The indices found represent the
%outtermost ROC curve/contour
for k=1:numel(u_FPR)
    cand_ind = find(u_FPR(k)>=FPR); %obtain the candidate indices
    [~,max_ind] = max(TPR(cand_ind));
    final_index(k) = cand_ind(max_ind);
end
optimal_indices = final_index;

function [confusion_training, sample_size_training, confusion_validation, sample_size_validation, config_mat_trainingvalidation] = getConfusionCrossValidation(handles)
%returns the confusion matrices and sample size for the training and
%validation sets;
%the config_mat_validation is also returned so it can be used by the
%interactive plotting
handles = excludeBadRows(handles); %remove unwanted rows based on exclusion criteria established in settings

[trainingControlStruct, validationControlStruct, trainingParamsStruct, validationParamsStruct] =getTrainingAndValidationSets(handles);

training_config_mat = handles.user.paramStruct.config_mat;
[confusion_training,sample_size_training] = getConfusion(handles, trainingControlStruct, trainingParamsStruct, training_config_mat);

[TPR,FPR,~,~, ~,~,~, ~,~,~, ~] = confusion2roc(confusion_training, sample_size_training);

%now want to obtain indices of unique FPR with matching maximum TPR
optimal_indices = getOutermostROCIndices(TPR,FPR);

config_mat_trainingvalidation = training_config_mat(optimal_indices,:);
confusion_training = confusion_training(optimal_indices,:);  %only get the best data points to cross-validate with

[confusion_validation,sample_size_validation] = getConfusion(handles, validationControlStruct, validationParamsStruct, config_mat_trainingvalidation);


function [confusion, sample_size] = getConfusionAll(handles)
%obtain confusion matrix for exhaustive search approach
%
% confusion= N x 4 array of confusion values where N is the number
% of combinations/thresholds applied
% sample_size = the final sample_size used to generate each confusion
% matrix
handles = excludeBadRows(handles); %remove unwanted rows based on exclusion criteria established in settings

params = handles.user.paramStruct.params;
controls = handles.user.controlStruct;
config_mat = handles.user.paramStruct.config_mat;
[confusion,sample_size] = getConfusion(handles, controls, params, config_mat);

function [confusion, sample_size] = getConfusion(handles, controls,params, config_mat)
%returns the confusion matrix and sample_size for the given controls and
%params structures
% controls fields include
%   .raw = a vector whose values are either .positive_diagnosis_value or
%   ~positive_diagnosis_value
%   .positive_diagnosis_value = the value in raw that represents positivity
% params is an array of structures (of size numEntries,1) whose  fields
% include
%    .range = vector of cutoff points
%    .binop = vector of binary operators that determine how .range is
%    evaluated against .raw
%    .raw = vector of symptom/parameter values obtained from the input file
%    .logic = vector of logic operations (AND, OR,etc) that establish how
%    to combine the next entry.
% obtained from the current gui settings/configuration
% confusion(configuration,:) = [TP,FN,FP,TN];

set(handles.fig1,'pointer','watch');
drawnow();

sample_size = size(controls.raw,1);
[num_configurations,num_parameters] = size(config_mat);
confusion = zeros(num_configurations,4); %[TP,FN,FP,TN]

%positive samples have a vaule of two, while negative diagnosis are zero; 
%this is done to help distinguish TP, FN, FP, and TN below when compared to
%the estimated diagnosis
diagnosis = 2*(controls.raw == controls.positive_diagnosis_value); 

group_cell = handles.user.grouping_cell;
if(isempty(group_cell))
    for c = 1:num_configurations       
        config = config_mat(c,:);
        for k=1:num_parameters
            %apply the cutoff criteria according to the inequality.  
            estimate = feval(params(k).binop,params(k).raw,config(k));
            
            if(k==1)
                running_estimate = estimate; %the first estimate; not combined with anything
            else
                %additional estimates are combined using the logic operator 
                %selected by the user
                running_estimate = feval(lower(params(k).logic),running_estimate,estimate);
            end
        end
        
        result = running_estimate+diagnosis;
        TP = sum(result==3)/sample_size; %b11
        FN = sum(result==2)/sample_size; %b10
        FP = sum(result==1)/sample_size; %b01
        TN = sum(result==0)/sample_size; %b00
        confusion(c,:) = [TP,FN,FP,TN];
    end
else
    num_groups = numel(group_cell);
    grouping_mat = false(num_parameters,2);
    for k = 1:num_groups
       grouping_mat(group_cell{k}(1),1) = true;
       grouping_mat(group_cell{k}(2),2) = true;
    end
    for c = 1:num_configurations       
        config = config_mat(c,:);
        group_on = false;
        for k=1:num_parameters
            estimate = feval(params(k).binop,params(k).raw,config(k));
            if(grouping_mat(k,1))
                group_on = true;
                if(k==1) %this is what will be applied after the grouping is done to recombine with the other set
                    group_logic = []; %no combination method on the first entry
                else
                    group_logic = lower(params(k).logic); %saves the logic combination method at the start point
                end
            end
            if(group_on)
                group_estimate = estimate;
            
                if(grouping_mat(k,1)) %if this is the first element in the group then there is not a running_group_estimate yet
                    running_group_estimate = group_estimate;
                else
                    running_group_estimate = feval(lower(params(k).logic),running_group_estimate,group_estimate);
                end
                if(grouping_mat(k,2)) % we got to the end of the grouping
                    group_on = false;
                    if(isempty(group_logic))
                        running_estimate = running_group_estimate;
                    else
                        running_estimate = feval(group_logic,running_estimate, running_group_estimate);
                    end
                end
            else
                if(k==1)
                    running_estimate = estimate;
                else
                    running_estimate = feval(lower(params(k).logic),running_estimate,estimate);
                end
            end
        end
        
        result = running_estimate+diagnosis;
        TP = sum(result==3)/sample_size; %b11
        FN = sum(result==2)/sample_size; %b10
        FP = sum(result==1)/sample_size; %b01
        TN = sum(result==0)/sample_size; %b00
        confusion(c,:) = [TP,FN,FP,TN];
    end
end
    
set(handles.fig1,'pointer','arrow');

function varargout = ROC_Callback(hObject,~,userdata,varargin)
%userdata is a struct with the following format
%   userdata.param_labels = handles.user.paramStruct.params;
%   userdata.config_mat = handles.user.paramStruct.config_mat;
%   userdata.TPR = TPR;
%   userdata.FPR = FPR;
%   userdata.K_0_0;
%   userdata.K_1_0;
%   userdata.CohensKappa = CohensKappa;
%   userdata.PPV = PPV;
%   userdata.NPV = NPV;
%   userdata.config_mat = config_mat;
%
%varargout{1} is config_ind = the index of the configuration that was selected
%varargin{1} is a preassigned configuration index, most likely from optimal
%point - if this is included, then the function does not need to search for
%a mouse click to obtain config_ind, but rather uses this value for it
if(numel(varargin)>0)
    config_ind = varargin{1};
else
    config_ind = [];
end

if(~strcmp(get(gcf,'selectiontype'),'alt') || ~isempty(config_ind))

    tag = get(hObject,'tag'); %can be tROC or qROC
    
    %find the point that the user clicked on
    if(isempty(config_ind))
        point = get(gca,'currentpoint');
        x = point(1,1);
        y = point(1,2);
        delta = 0.0075;
        deltaX = diff(get(gca,'xlim'))*delta;
        deltaY = diff(get(gca,'ylim'))*delta;
        if(strcmpi(tag,'tpr_vs_fpr'))
            hits = x<userdata.FPR+deltaX & x>userdata.FPR-deltaX & y<userdata.TPR+deltaY & y>userdata.TPR-deltaY;
            config_ind = find(hits,1,'first');
        elseif(strcmpi(tag,'se_vs_sp'))
            hits = x<(1-userdata.FPR+deltaX) & x>(1-userdata.FPR-deltaX) & y<userdata.TPR+deltaX & y>userdata.TPR-deltaY;
            config_ind = find(hits,1,'first');
        else
            hits = x<userdata.K_0_0+deltaX & x>userdata.K_0_0-deltaX & y<userdata.K_1_0+deltaY & y>userdata.K_1_0-deltaY;
            config_ind = find(hits,1,'first');
        end
    end
    %update all plots
    tROC_fig = userdata.tROC_fig;
    
    for k = 1:numel(tROC_fig)
        userdata = get(tROC_fig(k),'userdata'); %can have different userdata depending on the fig, specifically when showing training-validation plots
        
        if(isempty(config_ind))
            qroc_text_label = {'Sample not determined.','Move cursor closer to data point','and click again.'};
            roc_text_label = {'Sample not determined.','Move cursor closer to data point','and click again.'};
            set(userdata.roc_highlight_h,'visible','off');
            set(userdata.qroc_highlight_h,'visible','off');
        else
            
            plot_config = get(userdata.roc_highlight_h,'tag');
            decimal_format = ['%0.',num2str(userdata.decimal_places,'%d'),'f'];
            
            if(strcmpi(plot_config,'tpr_vs_fpr'))
                xlabelStr = 'False Positive Rate';
                ylabelStr= 'True Positive Rate';
                xdata = userdata.FPR(config_ind);
            else
                ylabelStr = 'Sensitivity';
                xlabelStr = 'Specificity';
                xdata = 1-userdata.FPR(config_ind);
            end
            roc_text_label = {[ylabelStr,' = ',num2str(userdata.TPR(config_ind),decimal_format)],[xlabelStr,' = ',num2str(xdata,decimal_format)]};
            
            ylabelStr = '\kappa(1,0)';
            xlabelStr = '\kappa(0,0)';
            
            qroc_text_label = {[ylabelStr,' = ',num2str(userdata.K_1_0(config_ind),decimal_format)],[xlabelStr,' = ',num2str(userdata.K_0_0(config_ind),decimal_format)]};
            
            for k = 1:numel(userdata.param_labels)
                qroc_text_label{end+1} = [userdata.param_labels{k},' = ',num2str(userdata.config_mat(config_ind,k),decimal_format)];
                roc_text_label{end+1} = [userdata.param_labels{k},' = ',num2str(userdata.config_mat(config_ind,k),decimal_format)];
            end
            set(userdata.roc_highlight_h,'visible','on','xdata',xdata,'ydata',userdata.TPR(config_ind));
            set(userdata.qroc_highlight_h,'visible','on','xdata',userdata.K_0_0(config_ind),'ydata',userdata.K_1_0(config_ind));
        end
        
        set(userdata.roc_text_h,'string',roc_text_label);
        pos = get(userdata.roc_text_h,'position');
        set(userdata.roc_text_h,'visible','on','position',[0.15, 0.15,pos(3:4)]);
        
        set(userdata.qroc_text_h,'string',qroc_text_label);
        pos = get(userdata.qroc_text_h,'position');
        set(userdata.qroc_text_h,'visible','on','position',[0.15, 0.15,pos(3:4)]);
    end
    
end
if(nargout>0)
    varargout{1} = config_ind;
end

function handles = excludeBadRows(handles)
%remove those rows which are excluded based on the settings set by the user

settings = handles.user.settings;
controls = handles.user.controlStruct;
params = handles.user.paramStruct.params;

if(handles.user.settings.exclude_number_data)

    bad_rows =  feval(settings.binop_number,controls.raw,settings.exclude_number_value);
    
    for k = 1:numel(params)
        bad_rows = bad_rows | feval(settings.binop_number,params(k).raw,settings.exclude_number_value);
    end
    if(isnumeric(settings.exclude_number_value))
        exclude_value = num2str(settings.exclude_number_value);
    else
        exclude_value = settings.exclude_number_value;
    end
    disp([num2str(sum(bad_rows)),' numeric rows removed due to exclusion criteria: remove numbers ',settings.binop_number,' ', exclude_value]);

    for k = 1:numel(params)    
        handles.user.paramStruct.params(k).raw = params(k).raw(~bad_rows,:);
    end
    handles.user.controlStruct.raw = controls.raw(~bad_rows,:);
end


% if(handles.user.settings.exclude_string_data)
% 
%     bad_rows =  strcmp(settings.exclude_string_value,controls.raw);
%     
%     for k = 1:numel(params)
%         bad_rows = bad_rows | strcmp(params(k).raw,settings.exclude_string_value);
%     end
%     
%     %invert answer if we have negative operation (i..e NOT equal)
%     if(strcmp(settings.binop_string,'~='))
%         bad_rows = ~bad_rows;
%     end
%     
%     disp([num2str(sum(bad_rows)),' text rows removed due to exclusion criteria: remove string data ',settings.binop_string,' ', settings.exclude_string_value]);
% 
%     for k = 1:numel(params)    
%         handles.user.paramStruct.params(k).raw = params(k).raw(~bad_rows,:);
%     end
%     handles.user.controlStruct.raw = controls.raw(~bad_rows,:);
%     
%    
% end


% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.user.grouping_cell = getGroupings(handles);

if(~areGroupingsOkay(handles))
    warndlg('Incorrect groupings found - will not continue');
else
    [pathstr, name, ext] = fileparts(handles.user.filename) ;
    suggested_filename = ['result - ',name,'.txt']; 
    [filename, pathname] = uiputfile('*.txt','Pick a text file',fullfile(handles.user.settings.current_pathname,suggested_filename));
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
    else
        %follow line delimiter convention based on operating system in use
        if(ispc)
            delim = '\r\n';
        else
            delim = '\n';
        end;
        
        decimal_format = ['%0.',num2str(handles.user.settings.decimal_places,'%d'),'f'];

        handles.user.paramStruct = getParameterConfigurations(handles);
        handles.user.controlStruct = getControlConfiguration(handles);
        [handles.user.confusion, sample_size] =  getConfusionAll(handles);
        
        [TPR,FPR,K_1_0,K_0_0, CohensKappa,PPV,NPV, EFF,chi_square,P, Q] = confusion2roc(handles.user.confusion, sample_size);
        
        fid = fopen(fullfile(pathname,filename),'w');
        
        fprintf(fid,'Sensitivity\tSpecificity\tCohensKappa\tPPV\tNPV\tTP\tFN\tFP\tTN\tK(1,0)\tK(0,0)\tEFF\tChi_square\tP\tQ');
        
        
        params = handles.user.paramStruct.params;
        config_mat = handles.user.paramStruct.config_mat;
        for k=1:numel(params)
            fprintf(fid,['\t',params(k).label]);
        end
        fprintf(fid,delim);
        fclose(fid);
        
        %group the data to be saved in tab-delimited, column order
        save_data = [TPR,1-FPR,CohensKappa,PPV,NPV,handles.user.confusion,K_1_0,K_0_0,EFF,chi_square,P,Q,config_mat];
        save(fullfile(pathname,filename),'save_data','-ascii','-tabs','-append');
        timestamp = datestr(now);        
        disp(['Data saved to ',fullfile(pathname,filename)]);
        
        log_fid = fopen(handles.user.settings.log_filename,'a');
        fprintf(log_fid,[delim,'%s\t%s',delim],timestamp,fullfile(pathname,filename));

        positive_diagnosis_value = handles.user.controlStruct.positive_diagnosis_value;
        
        %format for saving to file
        if(isnumeric(positive_diagnosis_value))
            positive_diagnosis_value = num2str(positive_diagnosis_value,decimal_format);
        end
        params = handles.user.paramStruct.params;
        group_cell = handles.user.grouping_cell;
        group_mat = zeros(numel(params),2);
        for k = 1:numel(group_cell)
           group_mat(group_cell{k}(1),1) = 1; 
           group_mat(group_cell{k}(2),2) = 1; 
        end
        fprintf(log_fid,['\tGold Standard = %s\tPositivity Value = %s',delim],handles.user.controlStruct.label,positive_diagnosis_value);
        for k=1:numel(params)
            start = params(k).range(1);
            stop = params(k).range(end);
            if(isnumeric(start))
                start = num2str(start,decimal_format);
                stop = num2str(stop,decimal_format);
            end
            group_label = '  ';
            if(group_mat(k,1))
               group_label(1) = '(';
            end
            if(group_mat(k,2))
               group_label(2) = ')';
            end
            
            % [(] Logic Label Binop [Start:Stop] [)]
            fprintf(log_fid,['\t%s %c\t%s\t%s [%s : %s] %c',delim],params(k).logic,group_label(1),params(k).label,params(k).binop,start,stop,group_label(2));
            
        end
        fclose(log_fid);
        
    end
end

% --- Executes on button press in push_settings.
function push_settings_Callback(hObject, eventdata, handles)
% hObject    handle to push_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

settings = softROC_settings(handles.user.settings);
fields = fieldnames(settings);

%this ensures I don't drop any fields from my user.settings that are not
%returned by softROC_settings
for k=1:numel(fields)
    handles.user.settings.(fields{k}) = settings.(fields{k});
end

% settings.pfile = handles.user.settings.pfile; %these fields are not returned by softROC_settings
% settings.current_pathname = handles.settings.current_pathname ;
    
plist.saveXMLPlist(handles.user.settings.pfile,handles.user.settings); %do this to save the most recent change to the pathname...
updateStepsText(handles);
guidata(hObject, handles);

function group_cell = getGroupings(handles)
%returns a cell whose elements contain a start stop vector of the tag
%suffixes/i.e. rows that begin and end the selected groupings
num_entries = handles.user.current_row;
cur_start = 0;
num_groups = 0;
group_cell = cell(num_entries,1);
failed = false;
for k = 1:num_entries
    start = get(handles.(['check_start',num2str(k)]),'value');
    stop = get(handles.(['check_end',num2str(k)]),'value');

    if(start)
        cur_start = k;
        failed = true; %avoid the problem of having a start with no matching stop
    end
    if(stop)
        if(cur_start==0)
            disp('failed!');
            failed = true;
        else
            num_groups = num_groups+1;
            group_cell{num_groups} = [cur_start,k]; 
            failed = false;
        end
    end
    
end

if(failed)
    group_cell = {};
else
    group_cell = group_cell(1:num_groups);
end




function check_start_Callback(hObject,eventdata)
handles = guidata(hObject);

handles.user.grouping_cell = getGroupings(handles);

if(~isempty(handles.user.grouping_cell))

else
    
end
guidata(hObject,handles);

function check_end_Callback(hObject,eventdata)
handles = guidata(hObject);
handles.user.grouping_cell = getGroupings(handles);
if(~isempty(handles.user.grouping_cell))

else
    
end
guidata(hObject,handles);


% --- Executes on button press in push_reset.
function push_reset_Callback(hObject, eventdata, handles)
% hObject    handle to push_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reset_gui(hObject,handles);


% --- Executes on button press in push_delete_entry.
function push_delete_entry_Callback(hObject, eventdata, handles)
% hObject    handle to push_delete_entry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
removeRow(handles,handles.user.current_row);
handles.user.current_row = handles.user.current_row-1;
if(handles.user.current_row == 1)
    set(hObject,'enable','off');
end;

guidata(hObject,handles);


% --- Executes when user attempts to close fig1.
function fig1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fig1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
try
    warning(handles.user.xlsread_warning_state,'MATLAB:xlsread:Mode');
catch    
    warning('on','MATLAB:xlsread:Mode');
end;
delete(hObject);
