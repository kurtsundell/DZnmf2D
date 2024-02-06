%% DZnmf %%

function varargout = DZnmf2D_v001(varargin)
% DZNMF2D_V001 MATLAB code for DZnmf2D_v001.fig
%      DZNMF2D_V001, by itself, creates a new DZNMF2D_V001 or raises the existing
%      singleton*.
%
%      H = DZNMF2D_V001 returns the handle to a new DZNMF2D_V001 or the handle to
%      the existing singleton*.
%
%      DZNMF2D_V001('CALLBACK',hObject,eventData,H,...) calls the local
%      function named CALLBACK in DZNMF2D_V001.M with the given input arguments.
%
%      DZNMF2D_V001('Property','Value',...) creates a new DZNMF2D_V001 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DZnmf2D_v001_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DZnmf2D_v001_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIH

% Edit the above text to modify the response to help DZnmf2D_v001

% Last Modified by GUIDE v2.5 29-Apr-2020 08:27:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DZnmf2D_v001_OpeningFcn, ...
                   'gui_OutputFcn',  @DZnmf2D_v001_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before DZnmf2D_v001 is made visible.
function DZnmf2D_v001_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for DZnmf2D_v001
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function varargout = DZnmf2D_v001_OutputFcn(hObject, eventdata, H) 
varargout{1} = H.output;

set(H.sigma1,'Enable','on')  %disable checkbox
set(H.sigma2,'Enable','on')  %disable checkbox
set(H.sigma1,'Value',1);

set(H.check_nmf,'Value',1)
set(H.check_nmf_opt,'Value',0);
set(H.num_sources_opt,'Enable','off');
set(H.max_num,'Enable','off');
set(H.opt_num,'Enable','off');
set(H.opt_num_result,'Enable','off');
set(H.num_sources,'Enable','on');
set(H.srcs,'Enable','on');
set(H.cancel,'Enable','off');
set(H.replot_res,'Enable','off');



set(H.pdps,'Value',0)
set(H.kdes,'value',1)
set(H.kern_opt,'Enable','off');
set(H.kern_set,'Enable','on');
set(H.kern,'Enable','on');

optimized = 0;
data = [];
N = [];
name = [];
text = [];
pdps_active = [];
density_active = [];
density_active_vector = [];
nsources = [];

global cancel
cancel = 0;

H.optimized = optimized;
H.data = data;
H.N = N;
H.name = name;
H.text = text;
H.pdps_active = pdps_active;
H.density_active = density_active;
H.density_active_vector = density_active_vector;
H.nsources = nsources;


%setting for default 
set(H.one_dimensional, 'Value', 0);
set(H.two_dimensional, 'Value', 0);
set(H.Dist_Opts_1D, 'Visible', 'off')
set(H.Dist_Opts_2D, 'Visible', 'off')
set(H.iterations,'String','1000');
set(H.tof,'string','1.5E-17');
set(H.Final_Residual,'String','0.015');
set(H.input_data, 'Visible', 'off');
set(H.Dist_Opts_1D, 'Visible', 'off')
set(H.Dist_Opts_2D, 'Visible', 'off', 'Position', [0.7 0.83 0.166 .17]);
set(H.Load, 'visible','off');
set(H.run_nmf, 'visible','off');
set(H.input_and_reconstructed_metrics, 'visible','off');
set(H.text61,'visible','off');
set(H.show_sinks, 'value',1);

set(H.stopping_criteria, 'visible','off');
set(H.calc_type,'visible','off');
set(H.compare_sources,'visible','off');
set(H.clear_all,'visible','off');
set(H.Export_NMF_Results,'visible','off');
set(H.export_plots,'visible','off');
set(H.uipanel20,'visible','off');


%initialize variables
H.density_active = [];
H.x_min = [];
H.x_max = [];
H.x_int = [];

guidata(hObject, H);

%% --- Executes during object creation, after setting all properties.
function Loaded_samples_CreateFcn(hObject, eventdata, H)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in show_sinks.
function show_sinks_Callback(hObject, eventdata, H)
if get(H.show_sinks,'value') == 1  
    set(H.hide_sinks,'value',0)
    set(H.stacked_sink_samples_uipanel, 'position', [.25 .05 .35 .75], 'visible','on')
    set(H.stacked_source_samples_uipanel, 'position', [.625 0.05 0.35 0.75], 'visible', 'on') 
    set(H.text19,'ForegroundColor',[0 0 0])
end



% --- Executes on button press in hide_sinks.
function hide_sinks_Callback(hObject, eventdata, H)
if get(H.hide_sinks,'Value') == 1  
    set(H.show_sinks,'Value',0)
    set(H.stacked_sink_samples_uipanel, 'position', [.25 .05 .35 .75], 'visible','off')
    set(H.stacked_source_samples_uipanel, 'position', [.25 0.05 0.725 0.75], 'visible', 'on') 
    set(H.text19,'ForegroundColor',[0.5 0.5 0.5])
end


%% CHECKBOX One Dimensional %%
function one_dimensional_Callback(hObject, eventdata, H)
if get(H.one_dimensional, 'Value') == 1
    set(H.two_dimensional, 'Value', 0);
    set(H.iterations,'String','10000');
    set(H.Final_Residual,'String','1E-8');
    set(H.input_data, 'Visible', 'on');
    set(H.Dist_Opts_2D, 'Visible', 'off')
    set(H.Dist_Opts_1D, 'Position', [0.7 0.83 0.166 .17]);
    set(H.Load,'visible','on');
    set(H.uipanel20,'visible','on');
    set(H.x_variable_name,'visible','on');
    set(H.text63,'visible','on')
    set(H.y_variable_name,'visible','off');
    set(H.text62,'visible','off');
    set(H.text45,'string','Mean / Range Kuiper V value of sink samples:');
end


%% CHECKBOX Two Dimensional %%
function two_dimensional_Callback(hObject, eventdata, H)
if get(H.two_dimensional, 'Value') == 1
    set(H.one_dimensional, 'Value', 0);
    set(H.iterations,'String','1000');
    set(H.Final_Residual,'String','0.015');
    set(H.tof,'string','1.5E-17');
    set(H.input_data, 'Visible', 'off');
    set(H.Dist_Opts_1D, 'Visible', 'off')
    set(H.Dist_Opts_2D, 'Position', [0.7 0.83 0.166 .17]);
    set(H.Load,'visible','on');
    set(H.uipanel20,'visible','on');
    set(H.x_variable_name,'visible','on');
    set(H.text63,'visible','on')
    set(H.y_variable_name,'visible','on');
    set(H.text62,'visible','on');
    set(H.text45,'string','Mean / Range KS D value of sink samples:');
    
end





%% CHECKBOX Load Ages %%
function load_ages_Callback(hObject, eventdata, H)

if get(H.load_ages,'Value') == 1
	set(H.load_dists,'Value',0);
	set(H.sigma1,'Enable','on')  %disable checkbox
	set(H.sigma2,'Enable','on')  %disable checkbox
	set(H.sigma1,'Value',1);
	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.minx,'Enable','on');
	set(H.maxx,'Enable','on');
	set(H.intx,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
end

%% CHECKBOX Load Distributions %%
function load_dists_Callback(hObject, eventdata, H)

if get(H.load_dists,'Value') == 1
	set(H.load_ages,'Value',0);
	set(H.sigma1,'Enable','off');
	set(H.sigma2,'Enable','off');
	set(H.sigma1,'Value',0);
	set(H.sigma2,'Value',0);
	set(H.xaxis_min,'Enable','off');
	set(H.xaxis_max,'Enable','off');
	set(H.xaxis_int,'Enable','off');
	set(H.minx,'Enable','off');
	set(H.maxx,'Enable','off');
	set(H.intx,'Enable','off');
	set(H.pdps,'Enable','off');
	set(H.kdes,'Enable','off');
end

%% CHECKBOX Run NMF %%
function check_nmf_Callback(hObject, eventdata, H)

if get(H.check_nmf,'Value') == 1
	set(H.check_nmf_opt,'Value',0);
	set(H.num_sources_opt,'Enable','off');
	set(H.max_num,'Enable','off');
	set(H.opt_num,'Enable','off');
	set(H.opt_num_result,'Enable','off');
	set(H.num_sources,'Enable','on');
	set(H.srcs,'Enable','on');
	set(H.cancel,'Enable','off');
	set(H.replot_res,'Enable','off');
end


%% CHECKBOX Run NMF and find optimal number of sources %%
function check_nmf_opt_Callback(hObject, eventdata, H)

if get(H.check_nmf_opt,'Value') == 1
	set(H.check_nmf,'Value',0);
	set(H.num_sources_opt,'Enable','on');
	set(H.max_num,'Enable','on');
	set(H.opt_num,'Enable','on');
	set(H.opt_num_result,'Enable','on');
	set(H.num_sources,'Enable','off');
	set(H.srcs,'Enable','off');
	set(H.cancel,'Enable','on');
	set(H.replot_res,'Enable','on');
end

%% CHECKBOX PDPs %%
function pdps_Callback(hObject, eventdata, H)

if get(H.pdps,'Value') == 1
	set(H.kdes,'Value',0)
	set(H.kern_opt,'Enable','off');
	set(H.kern_opt,'Value',0);
	set(H.kern_set,'Value',0);
	set(H.kern_set,'Enable','off');
	set(H.kern,'Enable','off');
	
end
%{
cla(H.stacked_reconstructed_sources,'reset');
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

pdps_active = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

ylabel([])
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])
%}
%active_selected = {1};
%loaded_selected = {1};
%H.active_selected = active_selected;
%H.loaded_selected = loaded_selected;
%H.pdps_active = pdps_active;
guidata(hObject, H);

%% CHECKBOX KDEs %%
function kdes_Callback(hObject, eventdata, H)

if get(H.kdes,'Value') == 1
	set(H.pdps,'Value',0)
	set(H.kern_opt,'Enable','on');
	set(H.kern_opt,'Value',1);
	set(H.kern_set,'Value',0);
	set(H.kern_set,'Enable','on');
	set(H.kern,'Enable','on');
end
%{
cla(H.stacked_reconstructed_sources,'reset');
set(H.stacked_reconstructed_sources,'ytick',[])
set(H.stacked_reconstructed_sources,'xtick',[])

cla(H.stacked_sink_samples,'reset');
set(H.stacked_sink_samples,'ytick',[])
set(H.stacked_sink_samples,'xtick',[])

pdps_active = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

ylabel([])
ytickformat('%.0f')
set(H.stacked_sink_samples,'ytick',[])
%}
%active_selected = {1};
%loaded_selected = {1};
%H.active_selected = active_selected;
%H.loaded_selected = loaded_selected;
%H.pdps_active = pdps_active;
guidata(hObject, H);

%% PUSHBUTTON Reset to Default Parameters
function reset_params_Callback(hObject, eventdata, H)


rad_on_dimensions=get(H.input_dimensions,'selectedobject');
switch rad_on_dimensions
    case H.one_dimensional
        set(H.iterations,'String','10000');
        set(H.Final_Residual,'String','1E-8');
        set(H.tof,'String','2E-16');
    case H.two_dimensional
        set(H.iterations,'String','1000');
        set(H.Final_Residual,'String','0.015');
        set(H.tof,'String','1.5E-17');
end



%% PUSHBUTTON Load Data %%
function Load_Callback(hObject, eventdata, H)

pdps_active = H.pdps_active;
density_active = H.density_active;
density_active_vector = H.density_active_vector;
set(H.stopping_criteria, 'visible','on');
set(H.calc_type,'visible','on');
if get(H.one_dimensional,'value') == 1
    set(H.Dist_Opts_1D,'visible','on');
else
    set(H.Dist_Opts_2D,'visible','on');
    set(H.color_chkbx,'enable','on');
    set(H.contour_chkbx,'enable','on');
    set(H.contour_percentile_bx,'enable','on');
end
set(H.clear_all, 'visible','on');

if isempty(pdps_active) == 0
    %clear sink axes
	cla(H.stacked_sink_samples,'reset');
    set(H.stacked_sink_samples,'ytick',[])
    set(H.stacked_sink_samples,'xtick',[])
    %clear factorized axes
    cla(H.stacked_reconstructed_sources,'reset');
    set(H.stacked_reconstructed_sources,'ytick',[])
    set(H.stacked_reconstructed_sources,'xtick',[])
    pdps_active = [];
    name = H.name;
    set(H.Active_samples, 'String', []);
    set(H.Loaded_samples, 'String', name);
    set(H.Active_samples, 'Value',1);
    ylabel([])
    ytickformat('%.0f')
    set(H.stacked_sink_samples,'ytick',[])
end




data = H.data;
N = H.N;
name = H.name;
text = H.text;

if isempty(data) == 0
	[filename, pathname] = uigetfile({'*'},'File Selector');
	fullpathname = strcat(pathname, filename);
	[numbers, text_new, data_new_tmp] = xlsread(fullpathname);
	if isempty(text_new) == 1
		err_dlg=errordlg('Please add headers (names) to your input data. Press Example Input button to see an example.','Error!');
		waitfor(err_dlg);
	else
		new_path = [get(H.filepath,'String'), '   &   ', fullpathname];
		set(H.filepath,'String',new_path);
		data_new = cell2mat(data_new_tmp(3:end,:));
		data_new_length = max(length(data_new(:,1)),length(data(:,1)));
		data_new_N = length(data_new(1,:))/2;
		data_tmp = data;
		data = zeros(data_new_length,N*2+data_new_N*2);
		data(1:length(data_tmp(:,1)),1:N*2) = data_tmp;
		data(1:length(data_new(:,1)),N*2+1:N*2+data_new_N*2) = data_new;
		N = N + data_new_N;
		for i = 1:data_new_N
			name_new(i,1) = text_new(1,i*2-1);
		end
		name = [name;name_new];
		text = [text,text_new];
	end
end

if isempty(data) == 1
	[filename, pathname] = uigetfile({'*'},'File Selector');
	fullpathname = strcat(pathname, filename);
	[numbers, text, data_tmp] = xlsread(fullpathname);
	if isempty(text) == 1
		err_dlg=errordlg('Please add headers (names) to your input data. Press Example Input button to see an example.','Error!');
		waitfor(err_dlg);
	else
		set(H.filepath, 'String', fullpathname);
		data = cell2mat(data_tmp(2:end,:));
		[dataR,dataC]=size(data);
		N = (dataC/2); % Number of source samples
		clear name
	end
end

%set axes
variable1=[];
variable2=[];
for i=1:N
    variable1 = vertcat(variable1, data(:,2*i-1));
    variable2 = vertcat(variable2, data(:,2*i));
    
end
if get(H.one_dimensional, 'Value') == 1
    set(H.xaxis_min,'string', num2str(floor(min(min(variable1)/500))*500));
    set(H.xaxis_max,'string', num2str(ceil(max(max(variable1)/500))*500));
    
end
if get(H.two_dimensional, 'Value') == 1
    set(H.xaxis_min_2D, 'string', num2str(floor(min(min(variable1)/500))*500));
    set(H.xaxis_max_2D,'string', num2str(ceil(max(max(variable1)/500))*500));
    set(H.yaxis_min_2D,'string', num2str(floor(min(min(variable2)))));
    set(H.yaxis_max_2D,'string', num2str(ceil(max(max(variable2)))));
end
clear variable1, variable2;
%end set axes
	
	for i = 1:N
		name(i,1) = text(1,i*2-1);
        %{
        if size(text, 1)<2;
            err_dlg=errordlg('Please add headers (names & input variables) to your input data. Press Example Input button to see an example.','Error!');
            waitfor(err_dlg);
        end
        input_variables = text(2,1:2);
        %}
	end



if get(H.one_dimensional,'value') == 1
rad_on=get(H.input_data,'selectedobject');
switch rad_on
	case H.sigma2
	for i = 1:N
		data(:,i*2) = data(:,i*2)./2;
	end
end
end

set(H.Loaded_samples, 'String', name);



set(H.Loaded_samples,'Value', 1);


active_selected = {1};
loaded_selected = {1};
stacked_sink   = NaN;            %set variable for 2D uipanel "stacked_sink_2D"
stacked_source = NaN;            %set variable for 2D uipanel "stacked_source_2D"
H.stacked_sink = stacked_sink;
H.stacked_source = stacked_source;
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
H.density_active = density_active;
H.density_active_vector = density_active_vector;
H.name = name;
H.text = text;
H.N = N;
H.data = data;

guidata(hObject, H);

%% PUSHBUTTON Run NMF! %%
function run_nmf_Callback(hObject, eventdata, H)
opengl hardware
%cla(H.stacked_reconstructed_sources,'reset');
%set(H.stacked_reconstructed_sources,'ytick',[])
%set(H.stacked_reconstructed_sources,'xtick',[])
%cla(H.Final_Residual_v_Rank,'reset');
%cla(H.SSR_v_Rank,'reset');

input_variables = H.input_variables;

%clear uipanels with 2D density plots, if they exist
try child = allchild(H.stacked_source_samples_uipanel)
    delete (child); 
catch; 
end;

active_names = cellstr(get(H.Active_samples,'String'));
check = active_names{1,1};

if isempty(check) == 1
	err_dlg=errordlg('There are no active samples!','But wait...');
	waitfor(err_dlg);
else
    
set(H.xaxis_min,'Enable','off');
set(H.xaxis_max,'Enable','off');
set(H.xaxis_int,'Enable','off');
set(H.pdps,'Enable','off');
set(H.kdes,'Enable','off');
set(H.compare_sources,'visible','on');
set(H.Export_NMF_Results, 'visible','on');
set(H.input_and_reconstructed_metrics, 'visible', 'on');
set(H.text61,'visible','on');
set(H.export_plots,'visible','on');
    
X1=H.X1;
Y1=H.Y1;
pdps_active = H.pdps_active;
density_active = H.density_active;
name = H.name;

option.iterations = str2num(get(H.iterations,'String'));
option.residual = str2num(get(H.Final_Residual,'String'));
option.tof = str2num(get(H.tof,'String'));

%Get/set axes parameters
%1D
x_min = H.x_min;
x_max = H.x_max;
x_int = H.x_int;
x = x_min:x_int:x_max;
set(H.xaxis_min,'enable','off');
set(H.xaxis_max,'enable','off');
set(H.xaxis_int,'enable','off');
set(H.pdps,'enable','off');
set(H.kdes,'enable','off');
set(H.kern_opt,'enable','off');
set(H.kern_set,'enable','off');
set(H.kern,'enable','off');
%2D
grid = H.grid;
lower_Xlim = H.lower_Xlim;
lower_Ylim = H.lower_Ylim;
upper_Xlim = H.upper_Xlim;
upper_Ylim = H.upper_Ylim;
MIN_XY     = [lower_Xlim,lower_Ylim];
MAX_XY     = [upper_Xlim,upper_Ylim];
X_bandwidth= H.X_bandwidth;
Y_bandwidth= H.Y_bandwidth;
color_ramp = H.color_ramp;
contour_var = H.contour_var;
contour_percentile = H.contour_percentile;
set(H.xaxis_min_2D,'enable','off');
set(H.xaxis_max_2D,'enable','off');
set(H.bandwidth_x,'enable','off');
set(H.yaxis_min_2D,'enable','off');
set(H.yaxis_max_2D,'enable','off');
set(H.bandwidth_y,'enable','off');
set(H.gridspc,'enable','off');
%end get/set parameters



global iter_num
global cancel

cancel = 0;

if get(H.check_nmf, 'Value') == 1
	nsources = str2num(get(H.num_sources,'String'));
	iter_num = nsources;
    if get(H.one_dimensional, 'Value') == 1
        [source_PDP,Weightings,numIter,tElapsed,finalResidual]=DZnmf(pdps_active,nsources,option);
        factorized_density_2D = []; 
    end;                 %set factorized_density as empty (for exporting)
    if get(H.two_dimensional, 'Value') == 1
        for i=1:length(active_names)
            sink_density_vector(:,i)=reshape(density_active(:,:,i),1,[]); 
        end
        [factorized_density_vector,Weightings,numIter,tElapsed,finalResidual]=...
            DZnmf(sink_density_vector,nsources, option);    
        for i=1:nsources
            factorized_density_2D(:,:,i)=reshape(factorized_density_vector(:,i), [size(X1,1), size(X1,1)]);
        end
        pdps_active = sink_density_vector;      %set pdps_active equal to sink_density_vector for calculation of metrics below
        source_PDP = []                         %set source_PDP as empty (for exporting)
    end 
    
    while numIter==str2num(get(H.iterations,'String'))
        answer=questdlg('Factorization did not converge before maximum number of iterations was reached.  Run multiple iterations?', ...
            'Maximum number of iterations reached','Yes','No','Yes');
        switch answer
            case 'Yes'
               prompt='Input number of iterations';
               title='Input';
               dims=[1 50];
               default_input={'50'};
               repeat_iterations=str2num(cell2mat(inputdlg(prompt,title,dims,default_input)));
               if get(H.one_dimensional, 'Value') == 1
                    for i=1:repeat_iterations
                        [source_PDP(:,:,i),Weightings(:,:,i),numIter(i),tElapsed(i),finalResidual(i)]...
                            = DZnmf(pdps_active,nsources,option);
                    end
                    [minimum_finalResidual, min_index]=min(finalResidual);
                    source_PDP=source_PDP(:,:,min_index);
                    Weightings=Weightings(:,:,min_index);
                    numIter=numIter(min_index);
                    finalResidual=finalResidual(min_index);
               end
               if get(H.two_dimensional, 'Value') == 1
                    for i=1:repeat_iterations %run repeat NMF on 2D densities 
                        [factorized_density_vector(:,:,i),Weightings(:,:,i),numIter(i),tElapsed(i),finalResidual(i)]...
                            = DZnmf(sink_density_vector,nsources,option);
                    end
                    [minimum_finalResidual, min_index]=min(finalResidual); %determine the minimum final residual for multiple NMF runs
                    factorized_density_vector=factorized_density_vector(:,:,min_index);  %select density associated with minimum final residual
                    for i=1:nsources                                % reshape vector to 2D density
                        factorized_density_2D(:,:,i)=reshape(factorized_density_vector(:,i), [size(X1,1), size(X1,1)]);
                    end
                    Weightings=Weightings(:,:,min_index);   %select Weightings for minimum final residual
                    numIter=numIter(min_index);             %select iterations for minimum final residual
                    finalResidual=finalResidual(min_index); %select final residual for minimum final residual
                    
                end
            case 'No'
                opts.WindowStyle='modal';
                opts.Interpreter = 'none';
                f=errordlg('Treat results with caution', 'Error',opts);
                break
        end
    end
    
    %reconstruct sink age distributions and CDFs
    if get(H.one_dimensional, 'Value') == 1
        reconstructed_sinks_PDP = source_PDP*Weightings;
        reconstructed_sinks_CDF = cumsum(reconstructed_sinks_PDP);
        cdfs_active = cumsum(pdps_active);
    end
    if get(H.two_dimensional, 'Value') == 1
        reconstructed_sinks_PDP = factorized_density_vector*Weightings; %set reconstructed_sinks_PDP...
        %...equal to the vectorized densities for the purposes of calculating cross-correlation coefficient
        for i=1:size(active_names)
            reconstructed_density(:,:,i)=reshape(reconstructed_sinks_PDP(:,i), [size(X1,1), size(X1,1)]);
            %CDFx(:,:,i) = cumsum(reconstructed_density(:,:,i),1);    % take the CDF of x along y
            %reconstructed_density_cdf (:,:,i) = cumsum(CDFx(:,:,i),2)';   % take the CDF of y along x                
            %sink_density_CDFx(:,:,i) = cumsum(density_active(:,:,i),1);
            %sink_density_cdf (:,:,i) = cumsum(sink_density_CDFx(:,:,i),2)';
        end
        [trash, Kuiper, mean_V, range_V_sink,trash1,K_S, mean_D, range_D_sink] = ...
            KS_from_density(reconstructed_density, density_active, get(H.one_dimensional,'Value'),hObject);
               
        
        Kuiper = Kuiper';
        K_S = K_S';
        mean_Kuiper = mean(Kuiper,1);
        stdev_Kuiper = std(Kuiper,1);
        mean_K_S = mean(K_S,1);
        stdev_K_S = std(K_S,1);
    end
        
    
    %run metrics comparing input to reconstructed sinks
    %Cross-correlation
    [trash,number_of_sinks]=size(reconstructed_sinks_PDP);
    for i=1:number_of_sinks
        R2(1,i)=r2(reconstructed_sinks_PDP(:,i), pdps_active(:,i));
    end
    mean_R2=mean(R2);
    stdev_R2=std(R2);
        
    %Likeness
    for i=1:number_of_sinks
        Likeness(1,i)=like(reconstructed_sinks_PDP(:,i), pdps_active(:,i));
    end
    mean_Likeness=mean(Likeness);
    stdev_Likeness=std(Likeness);
        
    %Kuiper
    if get(H.one_dimensional, 'Value') == 1
        for i=1:number_of_sinks
            Kuiper(1,i)=max(reconstructed_sinks_CDF(:,i)-cdfs_active(:,i))+max(cdfs_active(:,i)-reconstructed_sinks_CDF(:,i));
        end
    end
    mean_Kuiper=mean(Kuiper);
    stdev_Kuiper=std(Kuiper);
        
    %KS
    if get(H.one_dimensional, 'Value') == 1
    for i=1:number_of_sinks
        K_S(1,i)=max(abs(reconstructed_sinks_CDF(:,i)-cdfs_active(:,i)));
    end
    end
    mean_K_S=mean(K_S);
    stdev_K_S=std(K_S);
    
    %%%Print comparison to GUI
    set(H.mean_R2, 'String', round(mean_R2,3));
    set(H.stdev_R2, 'String', round(stdev_R2,3));
    set(H.mean_Likeness, 'String', round(mean_Likeness,3));
    set(H.stdev_Likeness, 'String', round(stdev_Likeness,3));
    if get(H.one_dimensional, 'Value') == 1 % Use V values for 1 dimensional 
        set(H.mean_Kuiper, 'String', round(mean_Kuiper,3));
        set(H.stdev_Kuiper, 'String', round(stdev_Kuiper,3));
    end
    if get(H.two_dimensional, 'Value') == 1 % Don't Use V values for 2 dimensional
        set(H.mean_Kuiper, 'String', 'N/A');
        set(H.stdev_Kuiper, 'String', 'N/A');
    end
    set(H.mean_K_S, 'String', round(mean_K_S,3));
    set(H.stdev_K_S, 'String', round(stdev_K_S,3));
    set(H.mean_finalResidual, 'String', round(finalResidual,3));
    %%%end Print comparison to GUI
    optimized = 0;
    
end

%Check for Optimum Rank
if get(H.check_nmf_opt, 'Value') == 1	
%Set number of iterations for 1D or 2D approaches to avoid extremely
    %long runs for 2D approach
    if get(H.one_dimensional, 'Value') == 1
        option.iterations=500;
    end
    if get(H.two_dimensional, 'Value') == 1
        option.iterations=200;
    end
%Set other ending criteria for NMF so that the maximum number of iterations
%is always reached.  This ensures that a uniform criterion is always
%applied when calculating the factorization rank
    option.residual = 1E-20;
    option.tof = 1E-100;
    nsources = str2num(get(H.num_sources_opt,'String'));
    
%Run looped NMF "nsources" times
    if get(H.one_dimensional, 'Value') == 1
        
        [source_PDP_loop,Weightings_loop,finalResidual_loop,...
            coefficient_count_loop,R2_loop, numIter_loop]=DZnmf_loop(pdps_active,nsources, option);
        factorized_density_2D = [];
        factorized_density_optimum = [];
    end
    if get(H.two_dimensional, 'Value') == 1
        
        for i=1:length(active_names)
            sink_density_vector(:,i)=reshape(density_active(:,:,i),1,[]); 
        end
        [source_PDP_loop,Weightings_loop,finalResidual_loop,...
            coefficient_count_loop,R2_loop, numIter_loop]=DZnmf_loop(sink_density_vector,nsources, option);    
    end
    
	xdata=transpose(2:1:nsources);
	ydata=transpose(finalResidual_loop);
    %{
    %%test optimum rank using bi-cross-validation (Owen and Perry, 2009)
    BCVsources = nsources;
    %option.iterations=100;
    %option.residual = 1E-20;
    %option.tof = 1E-100;
    
    if get(H.one_dimensional, 'Value') == 1
        input_density = pdps_active; %for 1D samples 
    end
    if get(H.two_dimensional, 'Value') == 1
        input_density = sink_density_vector; %for 2D samples
    end
    %divide the input matrix into four matrices
    col = ceil(0.5*size(input_density, 2));
    row = ceil(0.5*size(input_density, 1));
    A_in = input_density(1:row, 1:col);
    B_in = input_density(1:row, col+1:end);
    C_in = input_density(row+1:end, 1:col);
    D_in = input_density(row+1:end,col+1:end);
    
    [A_factorized,A_Weightings,A_numIter,A_tElapsed,A_finalResidual]=DZnmf_loop(A_in,BCVsources,option);
    [B_factorized,B_Weightings,B_numIter,B_tElapsed,B_finalResidual]=DZnmf_loop(B_in,BCVsources,option);
    [C_factorized,C_Weightings,C_numIter,C_tElapsed,C_finalResidual]=DZnmf_loop(C_in,BCVsources,option);
    
    [D_factorized,D_Weightings,D_numIter,D_tElapsed,D_finalResidual]=DZnmf_loop(D_in,BCVsources,option);
    
    for i=2:BCVsources
        
        A_out_inv = pinv(A_factorized(:,1:i,i-1));
        A_Weightings_inv = pinv(A_Weightings(1:i,:,i-1));
        D_out = (C_in * A_Weightings_inv) * (A_out_inv * B_in);
        Residual_D(i-1,1) = norm(D_in - D_out, 'fro');
        
        B_out_inv = pinv(B_factorized(:,1:i,i-1));
        B_Weightings_inv = pinv(B_Weightings(1:i,:,i-1));
        C_out = (D_in * B_Weightings_inv) * (B_out_inv * A_in);
        Residual_C(i-1,1) = norm(C_in - C_out, 'fro');
        
        C_out_inv = pinv(C_factorized(:,1:i,i-1));
        C_Weightings_inv = pinv(C_Weightings(1:i,:,i-1));
        B_out = (A_in * C_Weightings_inv) * (C_out_inv * D_in);
        Residual_B(i-1,1) = norm(B_in - B_out, 'fro');
        
        D_out_inv = pinv(D_factorized(:,1:i,i-1));
        D_Weightings_inv = pinv(D_Weightings(1:i,:,i-1));
        A_out = (B_in * D_Weightings_inv) * (D_out_inv * C_in);
        Residual_A(i-1,1) = norm(A_in - A_out, 'fro');
    end
    Residual_sum = Residual_A + Residual_B + Residual_C + Residual_D;
    [min_CV_Tot,optimum_rank_CV_Tot] = min(Residual_sum);
    [min_CV_A,optimum_rank_CV_A] = min(Residual_A);
    [min_CV_B,optimum_rank_CV_B] = min(Residual_B);
    [min_CV_C,optimum_rank_CV_C] = min(Residual_C);
    [min_CV_D,optimum_rank_CV_D] = min(Residual_D);
    
    optimum_rank_CV_Tot = optimum_rank_CV_Tot + 1;
    optimum_rank_CV_A   = optimum_rank_CV_A + 1;
    optimum_rank_CV_B   = optimum_rank_CV_B + 1;
    optimum_rank_CV_C   = optimum_rank_CV_C + 1;
    optimum_rank_CV_D   = optimum_rank_CV_D + 1;
    BCV_analysis = horzcat(xdata,Residual_A,Residual_B,Residual_C,Residual_D, Residual_sum);
    BCV_analysis = num2cell(BCV_analysis);
    BCV_text = {'Rank','Residual_A', 'Residual_B', 'Residual_C', 'Residual_D', 'Total_Residual'}
    BCV_analysis = vertcat(BCV_text, BCV_analysis);
    H.BCV_analysis = BCV_analysis;
    %%END "test optimum rank using bi-cross-validation (Owen and Perry,
    %2009)"
    %}
	%test integer breakpoint
	i=1;
	clear breakpoint brkpnt_coef SSR
	for xb=2:max(xdata)
		P0=[2 2 2];
		plusfun = @(x,xb)min(x,xb);
		model = @(P,x) P(1) + P(2)*x + P(3)*(plusfun(x,xb)-xb);
		Pfit = lsqcurvefit(model,P0,xdata,ydata);
		modelpred = model(Pfit,sort(xdata));
		squared_residuals=(modelpred-ydata).^2;
		SSR(i)=max(cumsum(squared_residuals));
		i=i+1;
	end
	[trash,brkpnt_coef]=min(SSR);
	breakpoint=xdata(brkpnt_coef);
        
    [trash,number_of_sinks, number_of_loops]=size(Weightings_loop);
    %reconstruct sink age distributions and CDFs
    if get(H.one_dimensional, 'Value') == 1
        for i=1:number_of_loops
            reconstructed_sinks_PDP(:,:,i)=source_PDP_loop(:,:,i)*Weightings_loop(:,:,i);
            reconstructed_sinks_CDF(:,:,i)=cumsum(reconstructed_sinks_PDP(:,:,i));
        end
        cdfs_active=cumsum(pdps_active);
    end
    if get(H.two_dimensional, 'Value') == 1
        for j=1:number_of_loops
        reconstructed_density = [];
        reconstructed_sinks_PDP(:,:,j)=source_PDP_loop(:,:,j)*Weightings_loop(:,:,j); %set reconstructed_sinks_PDP...
        for i = 1: number_of_sinks;
            %equal to the vectorized densities for the purposes of calculating cross-correlation coefficient
            reconstructed_density(:,:,i)=reshape(reconstructed_sinks_PDP(:,i,j), [size(X1,1), size(X1,1)]);
            
            %{
            CDFx(:,:,i) = cumsum(sum(reconstructed_density(:,:,i),1));    % take the CDF of x at y = max (y)
            reconstructed_density_cdf (:,:,i) = cumsum(CDFx(:,:,i),2)';                   % take the outer product of the two CDFs
            sink_density_CDFx(:,:,i) = cumsum(sum(density_active(:,:,i),1));
            sink_density_cdf (:,:,i) = cumsum(sink_density_CDFx(:,:,i),2)';
            Kuiper(j,i)=max(reconstructed_density_cdf(:,i,j)-sink_density_cdf(:,i,j))+...
                max(sink_density_cdf(:,i, j)-reconstructed_density_cdf(:,i,j));
            K_S(j,i)=max(abs(reconstructed_density_cdf(:,i,j)-sink_density_cdf(:,i, j)));
            %}
        end
        [trash, Kuiper(:,j), mean_V(j), range_V_sink(j),trash1,K_S(:,j), mean_D(j), range_D_sink(j)] = ...
            KS_from_density(reconstructed_density, density_active, get(H.one_dimensional,'Value'),hObject);
        %mean_Kuiper(j)=mean(Kuiper(j,:));
        %stdev_Kuiper(j)=std(Kuiper(j, :));
        %mean_K_S(j)=mean(K_S(j,:));
        %stdev_K_S(j)=std(K_S(j,:));
        
        end
        Kuiper = Kuiper';
        K_S = K_S';
        mean_Kuiper = mean(Kuiper,1);
        stdev_Kuiper = std(Kuiper,1);
        mean_K_S = mean(K_S,1);
        stdev_K_S = std(K_S,1);
    end
    
    
    %run metrics comparing input to reconstructed sinks
    
    %Cross-correlation
    for j=1:number_of_loops
        for i=1:number_of_sinks
            if get(H.one_dimensional, 'Value') == 1
            R2(j,i)=r2(reconstructed_sinks_PDP(:,i,j), pdps_active(:,i));
            end
            if get(H.two_dimensional, 'Value') == 1
                R2(j,i)=r2(reconstructed_sinks_PDP(:,i,j), sink_density_vector(:,i));
            end
        end
        mean_R2(j)=mean(R2(j,:));
        stdev_R2(j)=std(R2(j,:));
    end
    
    %Likeness
    for j=1:number_of_loops
        for i=1:number_of_sinks
            if get(H.one_dimensional, 'Value') == 1
                Likeness(j,i)=like(reconstructed_sinks_PDP(:,i,j), pdps_active(:,i));
            end
            if get(H.two_dimensional, 'Value') == 1
                Likeness(j,i)=like(reconstructed_sinks_PDP(:,i,j), sink_density_vector(:,i));
            end
        end
        mean_Likeness(j)=mean(Likeness(j,:));
        stdev_Likeness(j)=std(Likeness(j,:));
    end
    %Kuiper
    if get(H.one_dimensional, 'Value') == 1
    for j=1:number_of_loops
        for i=1:number_of_sinks
            Kuiper(j,i)=max(reconstructed_sinks_CDF(:,i,j)-cdfs_active(:,i))+max(cdfs_active(:,i)-reconstructed_sinks_CDF(:,i,j));
        end
        mean_Kuiper(j)=mean(Kuiper(j,:));
        stdev_Kuiper(j)=std(Kuiper(j, :));
    end
    end
    %KS
    if get(H.one_dimensional, 'Value') == 1
    for j=1:number_of_loops
        for i=1:number_of_sinks
        K_S(j,i)=max(abs(reconstructed_sinks_CDF(:,i,j)-cdfs_active(:,i)));
        end
        mean_K_S(j)=mean(K_S(j,:));
        stdev_K_S(j)=std(K_S(j,:));
    end
    end
    
        
    source_PDP = source_PDP_loop(:,:,brkpnt_coef);
	nsources = breakpoint;
    
    %%%Print comparison metrics to GUI
    set(H.mean_R2, 'String', round(mean_R2(brkpnt_coef),3));
    set(H.stdev_R2, 'String', round(stdev_R2(brkpnt_coef),3));
    set(H.mean_Likeness, 'String', round(mean_Likeness(brkpnt_coef),3));
    set(H.stdev_Likeness, 'String', round(stdev_Likeness(brkpnt_coef),3));
    set(H.mean_Kuiper, 'String', round(mean_Kuiper(brkpnt_coef),3));
    set(H.stdev_Kuiper, 'String', round(stdev_Kuiper(brkpnt_coef),3));
    set(H.mean_K_S, 'String', round(mean_K_S(brkpnt_coef),3));
    set(H.stdev_K_S, 'String', round(stdev_K_S(brkpnt_coef),3));
    if get(H.one_dimensional, 'Value') == 1 % Use V values for 1 dimensional 
        set(H.mean_Kuiper, 'String', round(mean_Kuiper(brkpnt_coef),3));
        set(H.stdev_Kuiper, 'String', round(stdev_Kuiper(brkpnt_coef),3));
    end
    if get(H.two_dimensional, 'Value') == 1 % Don't Use v values for 2 dimensional
        set(H.mean_Kuiper, 'String', 'N/A');
        set(H.stdev_Kuiper, 'String', 'N/A');
    end
    set(H.mean_finalResidual, 'String', round(finalResidual_loop(brkpnt_coef),3));
    %%%end Print comparison metrics to GUI
    
    if get(H.one_dimensional, 'Value') == 1 
        plot_1D_distribution(source_PDP(:,1:breakpoint), x_min, x_max, x_int, H.stacked_source_samples_uipanel, input_variables);
    end
    
    if get(H.two_dimensional, 'Value') == 1
        
        for j=1:size(source_PDP_loop,3)
            for i=1:size(source_PDP_loop,2)
            factorized_density_2D (:,:,i)=reshape(source_PDP_loop(:,i,j), [size(X1,1), size(X1,1)]);
            rows = size(factorized_density_2D,2);
            name = cell(1,rows);
            name(1,1) = cellstr(strcat('Density ',' ', num2str(i)));
            density_temp = num2cell(factorized_density_2D(:,:,i));
            export_temp(:,:,i) = vertcat(name,density_temp);
            end
            if j == brkpnt_coef
                plot_2D_distribution(factorized_density_2D(:,:,1:breakpoint), X1, Y1, H.stacked_source_samples_uipanel, input_variables, color_ramp, contour_var, contour_percentile);
                factorized_density_optimum = factorized_density_2D(:,:,1:breakpoint);
            end
            
            export_2D (:,:,j)= reshape(export_temp, [], size(export_temp,2),1);
            
        end
    end
    
end




if get(H.one_dimensional, 'Value') == 1 && get(H.check_nmf, 'Value') == 1
    nsamples = length(pdps_active(1,:));
    plot_1D_distribution(source_PDP, x_min, x_max, x_int,...
            H.stacked_source_samples_uipanel, input_variables);
    factorized_density_optimum = [];
end

if get(H.two_dimensional, 'Value') == 1  && get(H.check_nmf, 'Value') == 1
    plot_2D_distribution(factorized_density_2D, X1, Y1, H.stacked_source_samples_uipanel, ...
        input_variables, color_ramp, contour_var, contour_percentile);
    factorized_density_optimum = [];
end

if get(H.check_nmf_opt, 'Value') == 1	
	figure;
	xb=breakpoint;
	set(H.opt_num_result,'String',xb);
	plusfun = @(x,xb)min(x,xb);
	model = @(P,x) P(1) + P(2)*x + P(3)*(plusfun(x,xb)-xb);
	P0=[2 2 2];
	Pfit = lsqcurvefit(model,P0,xdata,ydata);
	modelpred = model(Pfit,sort(xdata));
	plot(xdata,ydata,'o',sort(xdata),modelpred,'r-');
	%title('Final Residual versus Rank')
	xlabel('Rank')
	ylabel('Final Residual')

	%plot SSR vs rank
	%axes(H.SSR_v_Rank);
	figure;
	plot(xdata,SSR,'-o');
	
	%title('SSR versus Rank')
	xlabel('Rank')
	ylabel('SSR1 + SSR2')
	
	optimized = 1;
	%concatenate data for output
	breakpoint_output=horzcat(xdata,ydata,modelpred,transpose(SSR));
    breakpoint_output = num2cell(breakpoint_output);
    breakpoint_text = {'Rank', 'Final Residual', 'Linear fit to final residual', 'Sum of squared residuals of linear fit'};
    breakpoint_output = vertcat(breakpoint_text, breakpoint_output);
    H.breakpoint_output = breakpoint_output;
	H.xdata = xdata;
	H.ydata = ydata;
	H.SSR = SSR;
	H.modelpred = modelpred;
end

if get(H.check_nmf, 'Value') == 1 
	
    export{length(x)+2,nsources+9+number_of_sinks} = [];
    export(1,1) = {'Reconstructed Source Distributions'};
	if get(H.one_dimensional, 'Value') == 1
        export(2,1) = {'Age (Ma)'};
        export(3:end,1) = num2cell(x');
        BaseName='Source ';
        for i=1:nsources
            export{2,i+1}=[BaseName,num2str(i)];
        end
        export(3:length(x)+2,2:nsources+1) = num2cell(source_PDP);
        export_2D = [];
    end
    if get(H.two_dimensional, 'Value') == 1
        export(1,2) = {'See associated sheet'};
        export(2,1) = {'x variable'};
        export(2,2) = input_variables(1,1);
        export(3,1) = {'y variable'};
        export(3,2) = input_variables(1,2);
        export(4,1) = {'x min'};
        export(4,2) = num2cell(str2num(get(H.xaxis_min_2D,'String')));
        export(5,1) = {'x max'};
        export(5,2) = num2cell(str2num(get(H.xaxis_max_2D,'String')));
        export(6,1) = {'y min'};
        export(6,2) = num2cell(str2num(get(H.yaxis_min_2D,'String')));
        export(7,1) = {'y max'};
        export(7,2) = num2cell(str2num(get(H.yaxis_max_2D,'String')));
        export_2D = factorized_density_2D;
    end
	export(1,nsources+4) = {'Compare Input & Reconstructed Sink Samples'};
	export(3,nsources+4) = {'Cross-correlation'};
	export(4,nsources+4) = {'Likeness'};
	export(5,nsources+4) = {'Kuiper V Value'};
	export(6,nsources+4) = {'KS D Value'};
	export(8,nsources+4) = {'Number of sink samples'};
    export(9,nsources+4) = {'Mean Cross-correlation of sink samples'};
    export(10,nsources+4)= {'Range Cross-correlation of sink samples'};
    export(11,nsources+4)= {'Mean Kuiper V value of sink samples'};
    export(12,nsources+4)= {'Range Kuiper V value of sink samples'};
    export(13,nsources+4)= {'Final Residual'};
    export(14,nsources+4) = {'Number of Iterations '};
	export(2,nsources+5) = {'Mean'};
	export(3,nsources+5) = num2cell(mean_R2);
	export(4,nsources+5) = num2cell(mean_Likeness);
	export(5,nsources+5) = num2cell(mean_Kuiper);
	export(6,nsources+5) = num2cell(mean_K_S);
    export(8,nsources+5) = num2cell(str2double(get(H.nsamples,'String')));
    export(9,nsources+5) = num2cell(str2double(get(H.mean_R2_sink,'String')));
    export(10,nsources+5)= num2cell(str2double(get(H.range_R2_sink,'String')));
    export(11,nsources+5)= num2cell(str2double(get(H.mean_V_sink,'String')));
    export(12,nsources+5)= num2cell(str2double(get(H.range_V_sink,'String')));
	export(13,nsources+5) = num2cell(finalResidual);
    export(14,nsources+5) = num2cell(numIter);
	export(3:6,nsources+6) = {'±'};
	export(2,nsources+7) = {'St. Dev.'};
	export(3,nsources+7) = num2cell(stdev_R2);
	export(4,nsources+7) = num2cell(stdev_Likeness);
	export(5,nsources+7) = num2cell(stdev_Kuiper);
	export(6,nsources+7) = num2cell(stdev_K_S);
	BaseName='Weights Source ';
	for i=1:nsources
		export{i+2,nsources+9}=[BaseName,num2str(i)];
	end
	export(nsources+4,nsources+9) = {'Cross-correlation'};
	export(nsources+5,nsources+9) = {'Likeness'};
	export(nsources+6,nsources+9) = {'Kuiper V Value'};
	export(nsources+7,nsources+9) = {'KS D Value'};
	export(1,nsources+10) = {'Sink Sample Weights'};
	export(3:nsources+2,nsources+10:nsources+number_of_sinks+9) = num2cell(Weightings);
	export(nsources+4,nsources+10:nsources+number_of_sinks+9) = num2cell(R2);
	export(nsources+5,nsources+10:nsources+number_of_sinks+9) = num2cell(Likeness);
	export(nsources+6,nsources+10:nsources+number_of_sinks+9) = num2cell(Kuiper);
	export(nsources+7,nsources+10:nsources+number_of_sinks+9) = num2cell(K_S);
	export(2,nsources+10:nsources+number_of_sinks+9) = cellstr(get(H.Active_samples,'String'));
end

if get(H.check_nmf_opt, 'Value') == 1
	tst{1,iter_num} = [];
	tst2{1,iter_num} = [];
	BaseName='Source ';
	BaseName2='Weights Source ';
	for i=1:iter_num
		tst{1,i}=[BaseName,num2str(i)];
		tst2{1,i}=[BaseName2,num2str(i)];
	end
	export{length(x)+2,nsources+9+number_of_sinks,iter_num-1} = [];

	for k = 1:iter_num-1
		export(1,1,k) = {'Reconstructed Source Distributions'};
		export(2,1,k) = {'Age (Ma)'};
		
        
        if get(H.one_dimensional, 'Value') == 1
            export(2,1) = {'Age (Ma)'};
            export(3:end,1,k) = num2cell(x');
            export(2,2:k+2,k) = tst(1,1:k+1);
            export(3:length(x)+2,2:k+2,k) = num2cell(source_PDP_loop(:,1:k+1,k));
            export_2D = [];
        end
        if get(H.two_dimensional, 'Value') == 1
        export(1,2,k) = {'See associated sheets'};
        export(2,1,k) = {'x variable'};
        export(2,2,k) = input_variables(1,1);
        export(3,1,k) = {'y variable'};
        export(3,2,k) = input_variables(1,2);
        export(4,1,k) = {'x min'};
        export(4,2,k) = num2cell(str2num(get(H.xaxis_min_2D,'String')));
        export(5,1,k) = {'x max'};
        export(5,2,k) = num2cell(str2num(get(H.xaxis_max_2D,'String')));
        export(6,1,k) = {'y min'};
        export(6,2,k) = num2cell(str2num(get(H.yaxis_min_2D,'String')));
        export(7,1,k) = {'y max'};
        export(7,2,k) = num2cell(str2num(get(H.yaxis_max_2D,'String')));
        export_2D = export_2D;
        end
    
		export(1,k+5,k) = {'Compare Input & Reconstructed Sink Samples'};
		export(3,k+5,k) = {'Cross-correlation'};
		export(4,k+5,k) = {'Likeness'};
		export(5,k+5,k) = {'Kuiper V Value'};
		export(6,k+5,k) = {'KS D Value'};
		export(8,k+5,k) = {'Number of sink samples'};
        export(9,k+5,k) = {'Mean Cross-correlation of sink samples'};
        export(10,k+5,k)= {'Range Cross-correlation of sink samples'};
        export(11,k+5,k)= {'Mean Kuiper V value of sink samples'};
        export(12,k+5,k)= {'Range Kuiper V value of sink samples'};
        export(13,k+5,k)= {'Final Residual'};
        export(14,k+5,k) = {'Number of Iterations '};
		export(2,k+6,k) = {'Mean'};
		export(3,k+6,k) = num2cell(mean_R2(1,k));
		export(4,k+6,k) = num2cell(mean_Likeness(1,k));
		export(5,k+6,k) = num2cell(mean_Kuiper(1,k));
		export(6,k+6,k) = num2cell(mean_K_S(1,k));
        export(8,k+6,k) = num2cell(str2double(get(H.nsamples,'String')));
        export(9,k+6,k) = num2cell(str2double(get(H.mean_R2_sink,'String')));
        export(10,k+6,k)= num2cell(str2double(get(H.range_R2_sink,'String')));
        export(11,k+6,k)= num2cell(str2double(get(H.mean_V_sink,'String')));
        export(12,k+6,k)= num2cell(str2double(get(H.range_V_sink,'String')));
		export(13,k+6,k) = num2cell(finalResidual_loop(1,k));
		export(14,k+6,k) = num2cell(numIter_loop(k));
		export(3:6,k+7,k) = {'±'};
		export(2,k+8,k) = {'St. Dev.'};
		export(3,k+8,k) = num2cell(stdev_R2(1,k));
		export(4,k+8,k) = num2cell(stdev_Likeness(1,k));
		export(5,k+8,k) = num2cell(stdev_Kuiper(1,k));
		export(6,k+8,k) = num2cell(stdev_K_S(1,k));
		export(3:k+3,k+10,k) = tst2(1,1:k+1);
		export(k+5,k+10,k) = {'Cross-correlation'};
		export(k+6,k+10,k) = {'Likeness'};
		export(k+7,k+10,k) = {'Kuiper V Value'};
		export(k+8,k+10,k) = {'KS D Value'};
		export(1,k+10,k) = {'Sink Sample Weights'};
		export(3:k+3,k+11:k+10+number_of_sinks,k) = num2cell(Weightings_loop(1:k+1,:,k));
		export(k+5,k+11:k+10+number_of_sinks,k) = num2cell(R2(k,:));
		export(k+6,k+11:k+10+number_of_sinks,k) = num2cell(Likeness(k,:));
		export(k+7,k+11:k+10+number_of_sinks,k) = num2cell(Kuiper(k,:));
		export(k+8,k+11:k+10+number_of_sinks,k) = num2cell(K_S(k,:));
		export(2,k+11:k+10+number_of_sinks,k) = cellstr(get(H.Active_samples,'String'));
	end
end






H.iter_num = iter_num;
H.nsources = nsources;
H.source_PDP = source_PDP;
H.factorized_density = factorized_density_2D;
H.factorized_density_optimum = factorized_density_optimum;
H.export = export;
H.export_2D = export_2D;
H.optimized = optimized;
end


guidata(hObject, H);

%% PUSHBUTTON Clear All %%
function clear_all_Callback(hObject, eventdata, H)

set(H.filepath,'String',[]);
set(H.Loaded_samples,'String',[]);
set(H.Active_samples,'String',[]);

pdps_active = [];                   %set pdp data set as empty
density_active = [];                %set density data set as empty
density_active_vector = [];

%clear uipanels with 2D density plots, if they exist
try child = allchild(H.stacked_sink_samples_uipanel)
    delete (child); 
catch; 
end;

try child = allchild(H.stacked_source_samples_uipanel)
    delete (child); 
catch; 
end;


%cla(H.Final_Residual_v_Rank,'reset');
%cla(H.SSR_v_Rank,'reset');
if get(H.one_dimensional, 'Value') == 1
        set(H.iterations,'String','10000');
        set(H.Final_Residual,'String','1E-8');
        set(H.tof,'String','2E-16');
end
if get(H.two_dimensional, 'Value') == 1
        set(H.iterations,'String','1000');
        set(H.Final_Residual,'String','0.015');
        set(H.tof,'String','1.5E-17');
end
set(H.one_dimensional,'value',0);
set(H.two_dimensional,'value',0);
set(H.Load, 'visible', 'off');
set(H.one_dimensional, 'Value', 0);
set(H.Dist_Opts_1D, 'Visible', 'off')

%reset 2D options
set(H.Dist_Opts_2D, 'Visible', 'off')
set(H.xaxis_min_2D,'enable','on');
set(H.xaxis_max_2D,'enable','on');
set(H.bandwidth_x,'enable','on');
set(H.yaxis_min_2D,'enable','on');
set(H.yaxis_max_2D,'enable','on');
set(H.bandwidth_y,'enable','on');
set(H.gridspc,'enable','on');

set(H.input_data, 'Visible', 'off');
set(H.Dist_Opts_1D, 'Visible', 'off')
set(H.Dist_Opts_2D, 'Visible', 'off', 'Position', [0.7 0.83 0.166 .17]);
set(H.Load, 'visible','off');
set(H.run_nmf, 'visible','off');
set(H.input_and_reconstructed_metrics, 'visible','off');
set(H.text61,'visible','off');
set(H.stopping_criteria, 'visible','off');
set(H.calc_type,'visible','off');
set(H.compare_sources,'visible','off');
set(H.clear_all,'visible','off');
set(H.Export_NMF_Results,'visible','off');
set(H.export_plots,'visible','off');
set(H.show_sinks, 'value', 1);
set(H.uipanel20,'visible','off');
set(H.x_variable_name,'enable','on');
set(H.y_variable_name,'enable','on');

set(H.tof,'String','2E-16');
set(H.opt_num_result,'String','?');
set(H.gridspc,'string','9');

set(H.xaxis_min,'string', num2str(0));
set(H.xaxis_max,'string', num2str(4000));
set(H.xaxis_min_2D, 'string', num2str(0));
set(H.xaxis_max_2D,'string', num2str(4000));
set(H.yaxis_min_2D,'string', num2str(-30));
set(H.yaxis_max_2D,'string', num2str(15));
set(H.xaxis_min,'enable','on');
set(H.xaxis_max,'enable','on');
set(H.xaxis_int,'enable','on');
set(H.pdps,'enable','on');
set(H.kdes,'enable','on');
set(H.kern_opt,'enable','on');
set(H.kern_set,'enable','on');
set(H.kern,'enable','on');


active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
H.density_active = density_active;
H.density_active_vector = density_active_vector;

optimized = 0;

data = [];
N = [];
name = [];
text = [];

H.optimized = optimized;
H.data = data;
H.N = N;
H.name = name;
H.text = text;

    set(H.mean_R2, 'String', 'N/A');
    set(H.stdev_R2, 'String', 'N/A');
    set(H.mean_Likeness, 'String', 'N/A');
    set(H.stdev_Likeness, 'String', 'N/A');
    set(H.mean_Kuiper, 'String', 'N/A');
    set(H.stdev_Kuiper, 'String', 'N/A');
    set(H.mean_K_S, 'String', 'N/A');
    set(H.stdev_K_S, 'String', 'N/A');
    set(H.mean_finalResidual, 'String', 'N/A');
    
guidata(hObject, H);

%% LISTBOX Loaded Samples %%
function Loaded_samples_Callback(hObject, eventdata, H)

loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_selected = {get(H.Loaded_samples,'Value')};

H.loaded_names = loaded_names;
H.loaded_selected = loaded_selected;
guidata(hObject, H);

%% LISTBOX Active Samples %%
function Active_samples_Callback(hObject, eventdata, H)

active_names = cellstr(get(H.Active_samples,'String'));
active_selected = {get(H.Active_samples,'Value')};

H.active_names = active_names;
H.active_selected = active_selected;
guidata(hObject, H);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PUSHBUTTON Activate Samples
function Activate_Callback(hObject, eventdata, H)
opengl hardware

pdps_active = H.pdps_active;
density_active = H.density_active;
density_active_vector = H.density_active_vector;
text = H.text;
N = H.N;
data = H.data;
loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_selected = H.loaded_selected;
stacked_sink = H.stacked_sink;

%clear uipanels with 2D density plots, if they exist
try child = allchild(H.stacked_sink_samples_uipanel)
    delete (child); 
catch; 
end;

if isempty(loaded_names) == 1
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(loaded_names)) == 0
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
else
	
%Get/set parameters
input_variables = {get(H.x_variable_name,'string') get(H.y_variable_name,'string')};
set(H.x_variable_name,'enable','off');
set(H.y_variable_name,'enable','off');
if size(density_active,1) == 0
    %1D
    x_min = str2num(get(H.xaxis_min,'String')); 
    x_max = str2num(get(H.xaxis_max,'String'));
    x_int = str2num(get(H.xaxis_int,'String'));
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
%2D
    grid = 2^str2num(get(H.gridspc,'String'));
    lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
    lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));
    upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));
    upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));
    MIN_XY=[lower_Xlim,lower_Ylim];
    MAX_XY=[upper_Xlim,upper_Ylim];
    X_bandwidth = str2num(get(H.bandwidth_x,'String'));
    Y_bandwidth = str2num(get(H.bandwidth_y,'String'));
    color_ramp = get(H.color_chkbx,'Value');
    contour_var = get(H.contour_chkbx,'Value');
    contour_percentile = str2num(get(H.contour_percentile_bx,'String'));
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    set(H.color_chkbx,'enable','off');
    set(H.contour_chkbx,'enable','off');
    set(H.contour_percentile_bx,'enable','off');
else
    %1D
    x_min = H.x_min;
    x_max = H.x_max;
    x_int = H.x_int;
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
    %2D
    grid = H.grid;
    lower_Xlim = H.lower_Xlim;
    lower_Ylim = H.lower_Ylim;
    upper_Xlim = H.upper_Xlim;
    upper_Ylim = H.upper_Ylim;
    MIN_XY     = [lower_Xlim,lower_Ylim];
    MAX_XY     = [upper_Xlim,upper_Ylim];
    X_bandwidth= H.X_bandwidth;
    Y_bandwidth= H.Y_bandwidth;
    contour_var = H.contour_var;
    color_ramp = H.color_ramp;
    contour_percentile = H.contour_percentile;
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    set(H.color_chkbx,'enable','off');
    set(H.contour_chkbx,'enable','off');
    set(H.contour_percentile_bx,'enable','off');
end
%end get/set parameters

set(H.run_nmf,'visible','on');

active_names = cellstr(get(H.Active_samples,'String'));
active_names(all(cellfun('isempty',active_names),2),:) = [];
active_names = [active_names; loaded_names(loaded_selected{1,1},1)];
set(H.Active_samples, 'String', active_names);

contents = get(H.Loaded_samples,'String');
if length(contents)<1; return; end %if already empty, do nothing
Index=get(H.Loaded_samples,'Value');
contents(Index)=[]; %remove the item
Value=Index-1;
if Value<1; Value=1;end %take care of exception
set(H.Loaded_samples,'String',contents,'Value',1);

new_names = loaded_names(loaded_selected{1,1});
%Load active samples' names ("new_names") to the Active samples list
for i = 1:N
	for j = 1:length(new_names)
		if strcmp(text(1,i*2-1), new_names(j,1)) == 1
			new_data_name{:,i*2-1} = new_names(j,1);
			new_data(:,i*2-1:i*2) = data(:,i*2-1:i*2);
		end
	end
end

%new_data(isnan(new_data))=0;
%new_data(isnan(new_data))=[];
nonzero = 1;
for i=1:size(new_data,2)
    if nnz(new_data(:,i))> 0
        temp(:,nonzero)=new_data(:,i);
        nonzero = nonzero + 1;
    end
end
new_data = temp; 

r=ones(1,size(new_data,1));
c=ones(1, size(new_data, 2));
new_data=mat2cell(new_data, r, c);
for i = 1:size(r,2)
	for j = 1:size(c,2)
		if cellfun('isempty', new_data(i,j)) == 0
			if cellfun(@isnan, new_data(i,j)) == 1
				new_data(i,j) = {[]};
			end	
		end
	end
end
%new_data=cell2mat(new_data);
%new_data( :, all(~new_data,1) ) = [];
new_data_name = new_data_name(~cellfun('isempty',new_data_name));

%calculate distributions for new samples ("new_data")
if get(H.one_dimensional, 'Value') == 1
    %Make PDPs          
    if get(H.pdps,'Value') == 1
        for i = 1:length(new_data_name)
            temp = cell2mat(new_data(:,i*2-1:i*2));
            
            m = nonzeros(temp(:,1));
    		m = m(isfinite(m(:,1)),:);
        	s = nonzeros(temp(:,2));
            s = s(isfinite(s(:,1)),:);
            if length(m)~=length(s)
                err_dlg=errordlg('Zeros in PDP bandwidth. Use KDEs instead');
                waitfor(err_dlg);
                return 
            end
            f = zeros(length(m),length(x));
            for j = 1:length(m)
            	f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
            end
            pdps_new(:,i) = ((sum(f, 1))/length(m)).';
        end
    end
    
    %Make KDEs
    if get(H.kdes,'Value') == 1
        rad_on_kernel=get(H.Dist_Opts_1D,'selectedobject');
        switch rad_on_kernel
            case H.kern_opt
                xA = transpose(x);
                for i = 1:length(new_data_name)
                    dist_data = nonzeros(cell2mat(new_data(:,i*2-1)));
                    [bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
                    pdps_new(:,i) = transpose(interp1(xmesh1, kdeA, xA));
                    clear dist_data kdeA
                end
            case H.kern_set
                for i = 1:length(new_data_name)
                    dist_data = nonzeros(cell2mat(new_data(:,i*2-1)));
                    kernel = str2num(get(H.kern,'String'));
                    kernel_dist_data(1:length(nonzeros(cell2mat(new_data(:,i*2-1)))),1) = kernel;
                    pdps_new(:,i) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
                end
        end
    end
    pdps_active = [pdps_active,pdps_new];
    density_active=[];
    input_distributions = pdps_active;
end
if get(H.two_dimensional, 'Value') == 1 %calculate density for 2D KDEs
    for i = 1:length(new_data_name)
        data_now=cell2mat(new_data(:,2*i-1:2*i));
        [bandwidth_not_needed,sink_density(:,:,i),X1,Y1]=kde2d(data_now, grid, MIN_XY, MAX_XY, X_bandwidth, Y_bandwidth);
        sink_density(:,:,i) = sink_density(:,:,i)./sum(sum(sink_density(:,:,i)));
        density_vector_new(:,i)=reshape(sink_density(:,:,i),1,[]);
    end
    density_new=sink_density;
    density_active=cat(3,density_active,density_new);
    
    density_active_vector = [density_active_vector, density_vector_new];
    
    pdps_active=[];
    input_distributions = density_active;
end
    
%plot all active samples ("pdp_active")
rad_on_plotting=get(H.Display_sinks,'selectedobject');
if get(H.show_sinks,'value') ==1
    
    if get(H.one_dimensional, 'Value') == 1
        nsamples = length(pdps_active(1,:));
        plot_1D_distribution(pdps_active, x_min, x_max, x_int,...
            H.stacked_sink_samples_uipanel, input_variables);
        X1=NaN;                     %set 2D elements as NaN
        Y1=NaN;                     %set 2D elements as NaN
              %set 2D elements as NaN
    end
    if get(H.two_dimensional, 'Value') == 1
        delete stacked_sink
        nsamples = size(density_active,3);
        nsamples_new = size(density_new,3);
        
        plot_2D_distribution(density_active, X1, Y1, ...
            H.stacked_sink_samples_uipanel, input_variables, color_ramp, contour_var, contour_percentile);   
        
    end
end
if get(H.hide_sinks, 'value') == 1
    if get(H.one_dimensional, 'Value') == 1
        nsamples = length(pdps_active(1,:));
        X1=NaN;                     %set 2D elements as NaN
        Y1=NaN;                     %set 2D elements as NaN
        stacked_sink_2D = NaN;      %set 2D elements as NaN
    end
    if get(H.two_dimensional, 'Value') == 1
        nsamples = size(density_active,3);
        
    end
end

%Mean Kuiper
H.input_variables = input_variables;
H.X1=X1;
H.Y1=Y1;
guidata(hObject, H);%pass variables to hObject
[trash, V, mean_V_sink, range_V_sink,trash1, D, mean_D_sink, range_D_sink] =...
    KS_from_density(input_distributions, input_distributions, get(H.one_dimensional,'Value'),hObject) 
if get(H.one_dimensional, 'Value') == 1 % Use V values for 1 dimensional 
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
if get(H.two_dimensional, 'Value') == 1 % Use D values for 2 dimensional
    set(H.mean_V_sink, 'String', mean_D_sink);
    set(H.range_V_sink, 'String', range_D_sink);
end


%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            if get(H.one_dimensional, 'Value') == 1
                mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));end
            if get(H.two_dimensional, 'Value') == 1
                mean_R2_sink(i,j)=r2(density_active_vector(:,i),density_active_vector(:,j)); end
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end
    
set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.stacked_sink = stacked_sink;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;            %exports the 1D pdps 
H.density_active = density_active       %exports the 2D densities  
H.density_active_vector = density_active_vector;
H.X1=X1;
H.Y1=Y1;
H.active_names = active_names;          %exports active sample names
H.density_active_vector = density_active_vector;
H.grid = grid;
H.lower_Xlim=lower_Xlim;
H.upper_Xlim=upper_Xlim;
H.lower_Ylim=lower_Ylim;
H.upper_Ylim=upper_Ylim;
H.X_bandwidth= X_bandwidth;
H.Y_bandwidth= Y_bandwidth;
H.contour_percentile = contour_percentile;
H.color_ramp = color_ramp;
H.contour_var = contour_var;
guidata(hObject, H);

end

%% PUSHBUTTON Deactivate Samples
function Deactivate_Callback(hObject, eventdata, H)
input_variables = H.input_variables;
pdps_active = H.pdps_active;
density_active = H.density_active;
text = H.text;
N = H.N;
X1=H.X1;
Y1=H.Y1;
data = H.data;
active_names = cellstr(get(H.Active_samples,'String'));
active_selected = H.active_selected;
stacked_sink = H.stacked_sink;

%clear uipanels with 2D density plots, if they exist
try child = allchild(H.stacked_sink_samples_uipanel)
    delete (child); 
catch; 
end;

%Get/set parameters
if size(density_active,1) == 0
%1D
    x_min = str2num(get(H.xaxis_min,'String')); 
    x_max = str2num(get(H.xaxis_max,'String'));
    x_int = str2num(get(H.xaxis_int,'String'));
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','on');
    set(H.xaxis_max,'enable','on');
    set(H.xaxis_int,'enable','on');
    set(H.pdps,'enable','on');
    set(H.kdes,'enable','on');
    set(H.kern_opt,'enable','on');
    set(H.kern_set,'enable','on');
    set(H.kern,'enable','on');
%2D
    grid = 2^str2num(get(H.gridspc,'String'));
    lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
    lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));
    upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));
    upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));
    MIN_XY=[lower_Xlim,lower_Ylim];
    MAX_XY=[upper_Xlim,upper_Ylim];
    X_bandwidth = str2num(get(H.bandwidth_x,'String'));
    Y_bandwidth = str2num(get(H.bandwidth_y,'String'));
    color_ramp = get(H.color_chkbx,'Value');
    contour_var = get(H.contour_chkbx,'Value');
    contour_percentile = str2num(get(H.contour_percentile_bx,'String'));
    set(H.xaxis_min_2D,'enable','on');
    set(H.xaxis_max_2D,'enable','on');
    set(H.bandwidth_x,'enable','on');
    set(H.yaxis_min_2D,'enable','on');
    set(H.yaxis_max_2D,'enable','on');
    set(H.bandwidth_y,'enable','on');
    set(H.gridspc,'enable','on');
    set(H.color_chkbx,'enable','on');
    set(H.contour_chkbx,'enable','on');
    set(H.contour_percentile_bx,'enable','on');
    
    
else
%1D
    x_min = H.x_min;
    x_max = H.x_max;
    x_int = H.x_int;
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
%2D
    grid = H.grid;
    lower_Xlim = H.lower_Xlim;
    lower_Ylim = H.lower_Ylim;
    upper_Xlim = H.upper_Xlim;
    upper_Ylim = H.upper_Ylim;
    MIN_XY     = [lower_Xlim,lower_Ylim];
    MAX_XY     = [upper_Xlim,upper_Ylim];
    X_bandwidth= H.X_bandwidth;
    Y_bandwidth= H.Y_bandwidth;
    contour_var = H.contour_var;
    color_ramp = H.color_ramp;
    contour_percentile = H.contour_percentile;
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    
end
%end get/set parameters

if isempty(active_names) == 1
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
    set(H.x_variable_name,'enable','on');
    set(H.y_variable_name,'enable','on');
elseif length(char(active_names)) == 0
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
else

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));

loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_names(all(cellfun('isempty',loaded_names),2),:) = [];
loaded_names = [loaded_names; active_names(active_selected{1,1},1)];
for i = 1:length(loaded_names)
	if length((strsplit(loaded_names{i,1},','))) > 1
		loaded_names(length(loaded_names)+1:length(loaded_names)+length((strsplit(loaded_names{i,1},','))),1) = strsplit(loaded_names{i,1},',');
		loaded_names{i,1} = [];
	end
end
loaded_names(all(cellfun('isempty',loaded_names),2),:) = [];

set(H.Loaded_samples, 'String', loaded_names);

contents = get(H.Active_samples,'String');
if length(contents)<1; return; end %if already empty, do nothing
Index=get(H.Active_samples,'Value');
contents(Index)=[]; %remove the item
Value=Index-1;
if Value<1; Value=1;end %take care of exception
set(H.Active_samples,'String',contents,'Value',1);

active_names = cellstr(get(H.Active_samples,'String'));

selected_rmv = active_selected{1,1};
if get(H.show_sinks, 'value') == 1
if get(H.one_dimensional, 'Value') == 1
    idx = 1:1:size(pdps_active, 2);
    idx(1,selected_rmv) = 0;
    pdps_active = pdps_active(:, find(idx));
    plot_1D_distribution(pdps_active, x_min, x_max, x_int,...
            H.stacked_sink_samples_uipanel, input_variables); 
    nsamples = length(pdps_active(1,:));
    input_distributions = pdps_active;
end
if get(H.two_dimensional, 'Value') == 1
    idx = 1:1:size(density_active, 3);
    idx(1,selected_rmv) = 0;
    density_active = density_active(:, :, find(idx));
    
    plot_2D_distribution(density_active, X1, Y1, ...
        H.stacked_sink_samples_uipanel, input_variables, color_ramp, contour_var, contour_percentile);
    nsamples = size(density_active, 3);
    input_distributions = density_active;
end
end
if get(H.hide_sinks, 'value') == 1
    %don't plot if sink plot is hidden
end



%Mean Kuiper
[trash, V, mean_V_sink, range_V_sink,trash1, D, mean_D_sink, range_D_sink] = ...
    KS_from_density(input_distributions, input_distributions, get(H.one_dimensional,'Value'),hObject) %calculate mean Kuiper based on function below
if get(H.one_dimensional, 'Value') == 1 % Use V values for 1 dimensional 
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
if get(H.two_dimensional, 'Value') == 1 % Use D values for 2 dimensional
    set(H.mean_V_sink, 'String', mean_D_sink);
    set(H.range_V_sink, 'String', range_D_sink);
end


%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    if get(H.one_dimensional, 'Value') == 1
        for i=1:nsamples
        for j=1:nsamples
            mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));
        end
        end
    end
    if get(H.two_dimensional, 'Value') == 1
        for i=1:nsamples
            density_active_vector(:,i)=reshape(density_active(:,:,i),1,[]);
        end
        
        for i=1:nsamples
        for j=1:nsamples
            
            mean_R2_sink(i,j)=r2(density_active_vector(:,i),density_active_vector(:,j));
        end
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end

set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
H.density_active=density_active;
guidata(hObject, H);


if isempty(active_names) == 1
	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
elseif length(char(active_names)) == 0
	set(H.xaxis_min,'Enable','on');
	set(H.xaxis_max,'Enable','on');
	set(H.xaxis_int,'Enable','on');
	set(H.pdps,'Enable','on');
	set(H.kdes,'Enable','on');
end

end

%% PUSHBUTTON Merge and Activate
function Merge_Callback(hObject, eventdata, H)
opengl hardware
pdps_active = H.pdps_active;
density_active = H.density_active;
X1=H.X1;
Y1=H.Y1;
text = H.text;
N = H.N;
data = H.data;
loaded_names = cellstr(get(H.Loaded_samples,'String'));
loaded_selected = H.loaded_selected;
density_active_vector = H.density_active_vector;
input_variables = H.input_variables;

set(H.xaxis_min,'Enable','off');
set(H.xaxis_max,'Enable','off');
set(H.xaxis_int,'Enable','off');
set(H.pdps,'Enable','off');
set(H.kdes,'Enable','off');

%Get/set parameters
input_variables = {get(H.x_variable_name,'string') get(H.y_variable_name,'string')};
set(H.x_variable_name,'enable','off');
set(H.y_variable_name,'enable','off');
if size(density_active,1) == 0
%1D variables
    x_min = str2num(get(H.xaxis_min,'String')); 
    x_max = str2num(get(H.xaxis_max,'String'));
    x_int = str2num(get(H.xaxis_int,'String'));
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
%2D variables
    grid = 2^str2num(get(H.gridspc,'String'));
    lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
    lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));
    upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));
    upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));
    MIN_XY=[lower_Xlim,lower_Ylim];
    MAX_XY=[upper_Xlim,upper_Ylim];
    X_bandwidth = str2num(get(H.bandwidth_x,'String'));
    Y_bandwidth = str2num(get(H.bandwidth_y,'String'));
    color_ramp = get(H.color_chkbx,'Value');
    contour_var = get(H.contour_chkbx,'Value');
    contour_percentile = str2num(get(H.contour_percentile_bx,'String'));
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    set(H.color_chkbx,'enable','off');
    set(H.contour_chkbx,'enable','off');
    set(H.contour_percentile_bx,'enable','off');
    
else
    %1D
    x_min = H.x_min;
    x_max = H.x_max;
    x_int = H.x_int;
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
    %2D
    grid = H.grid;
    lower_Xlim = H.lower_Xlim;
    lower_Ylim = H.lower_Ylim;
    upper_Xlim = H.upper_Xlim;
    upper_Ylim = H.upper_Ylim;
    MIN_XY     = [lower_Xlim,lower_Ylim];
    MAX_XY     = [upper_Xlim,upper_Ylim];
    X_bandwidth= H.X_bandwidth;
    Y_bandwidth= H.Y_bandwidth;
    contour_var = H.contour_var;
    color_ramp = H.color_ramp;
    contour_percentile = H.contour_percentile;
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    set(H.color_chkbx,'enable','off');
    set(H.contour_chkbx,'enable','off');
    set(H.contour_percentile_bx,'enable','off');
    
end
%end get/set parameters

if isempty(loaded_names) == 1
	err_dlg=errordlg('There are no samples to merge!','Hang on...');
	waitfor(err_dlg);
elseif length(char(loaded_names)) == 0
	err_dlg=errordlg('There are no samples to merge!','Hang on...');
	waitfor(err_dlg);
else
	
%clear uipanels with 2D density plots, if they exist
try child = allchild(H.stacked_sink_samples_uipanel);
    delete (child); 
catch; 
end;

x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));


tmp = loaded_names(loaded_selected{1,1});
result = strcat(tmp(1:end-1,1), ',');
result(end+1,1) = tmp(end,1);
namecat = result{1,1};
for i = 2:length(result)
	namecat = [namecat, result{i,1}];
end

active_names = cellstr(get(H.Active_samples,'String'));
active_names(all(cellfun('isempty',active_names),2),:) = [];
active_names = [active_names; namecat];
set(H.Active_samples, 'String', active_names);

contents = get(H.Loaded_samples,'String');
if length(contents)<1; return; end %if already empty, do nothing
Index=get(H.Loaded_samples,'Value');
contents(Index)=[]; %remove the item
Value=Index-1;
if Value<1; Value=1;end %take care of exception
set(H.Loaded_samples,'String',contents,'Value',1);

comb = (strsplit(active_names{end,1},','))';

for i = 1:N
	for j = 1:length(comb)
		if strcmp(text(1,i*2-1), comb(j,1)) == 1
			comb_data_name{:,i*2-1} = comb(j,1);
			comb_data(:,i*2-1:i*2) = data(:,i*2-1:i*2);
		end
	end
end

comb_data_name = comb_data_name(~cellfun('isempty',comb_data_name));

%calculate and plot for 1D case
if get(H.one_dimensional, 'Value') == 1
    comb_data(isnan(comb_data))=0;
    comb_data( :, all(~comb_data,1) ) = [];

    for i = 1:length(comb_data_name)
        comb_data_length(i,1) = length(nonzeros(comb_data(:,i*2-1)));
    end

    merged_data = [];

    for i = 1:length(comb_data_name)
        merged_data(1+length(merged_data):length(merged_data)+comb_data_length(i,1),1:2) = comb_data(1:comb_data_length(i,1),i*2-1:i*2);
    end

    x = x_min:x_int:x_max;

    if get(H.pdps,'Value') == 1
        m = nonzeros(merged_data(:,1));
        m = m(isfinite(m(:,1)),:);
        s = nonzeros(merged_data(:,2));
        s = s(isfinite(s(:,1)),:);
        f = zeros(length(m),length(x));
        for j = 1:length(m)
            f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
        end
        pdps_merged(:,1) = ((sum(f, 1))/length(m)).';
    end

    if get(H.kdes,'Value') == 1
        rad_on_kernel=get(H.Dist_Opts_1D,'selectedobject');
        switch rad_on_kernel
            case H.kern_opt
                xA = transpose(x);
                dist_data = nonzeros(merged_data(:,1));
                [bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
                pdps_merged(:,1) = transpose(interp1(xmesh1, kdeA, xA));
            case H.kern_set
                dist_data = nonzeros(merged_data(:,1));
                kernel = str2num(get(H.kern,'String'));
                kernel_dist_data(1:length(nonzeros(merged_data(:,1))),1) = kernel;
                pdps_merged(:,1) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
        end
    end


    pdps_active = [pdps_active,pdps_merged];
    input_distributions = pdps_active; %set variable for calculating Kuiper and Cross-correlation

    nsamples = length(pdps_active(1,:));
    if get(H.show_sinks,'value') == 1 %plot results if sinks are displayed
    plot_1D_distribution(pdps_active, x_min, x_max, x_int,...
            H.stacked_sink_samples_uipanel, input_variables);
    end
end

%calculate and plot for 2D case
if get(H.two_dimensional, 'Value') == 1
    comb_data(isnan(comb_data))=0;
    comb_data( :, all(~comb_data,1) ) = [];

    for i = 1:length(comb_data_name)
        comb_data_length(i,1) = length(nonzeros(comb_data(:,i*2-1)));
    end
    
    merged_data = [];

    for i = 1:length(comb_data_name)
        merged_data(1+length(merged_data):length(merged_data)+comb_data_length(i,1),1:2) = comb_data(1:comb_data_length(i,1),i*2-1:i*2);
    end
    grid = 2^str2num(get(H.gridspc,'String'));
    lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
    lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));;
    upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));;
    upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));;
    MIN_XY=[lower_Xlim,lower_Ylim];
    MAX_XY=[upper_Xlim,upper_Ylim];
       
        [bandwidth_not_needed,sink_density,X1,Y1]=kde2d(merged_data, grid, MIN_XY, MAX_XY, X_bandwidth, Y_bandwidth);
        density_new = sink_density./sum(sum(sink_density));
        density_vector_new=reshape(sink_density(:,:),1, []);
        
    density_active=cat(3,density_active,density_new);
    density_active_vector = [density_active_vector, density_vector_new'];
    
    pdps_active=[];
    input_distributions = density_active;
    
    if get(H.show_sinks, 'value') == 1
     plot_2D_distribution(density_active, X1, Y1, ...
         H.stacked_sink_samples_uipanel, input_variables, color_ramp, contour_var, contour_percentile);
    end
    nsamples = size(density_active,3);

end %end 2D case

%Mean Kuiper
[trash, V, mean_V_sink, range_V_sink,trash1, D, mean_D_sink, range_D_sink] =...
    KS_from_density(input_distributions, input_distributions, get(H.one_dimensional,'Value'),hObject) %calculate mean Kuiper based on function below
if get(H.one_dimensional, 'Value') == 1 % Use V values for 1 dimensional 
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
if get(H.two_dimensional, 'Value') == 1 % Use D values for 2 dimensional
    set(H.mean_V_sink, 'String', mean_D_sink);
    set(H.range_V_sink, 'String', range_D_sink);
end


%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    for i=1:nsamples
        for j=1:nsamples
            if get(H.one_dimensional, 'Value') == 1
                mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));end
            if get(H.two_dimensional, 'Value') == 1
                mean_R2_sink(i,j)=r2(density_active_vector(:,i),density_active_vector(:,j)); end
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end

set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
H.density_active_vector = density_active_vector;
H.density_active = density_active;
H.contour_percentile = contour_percentile;
H.color_ramp = color_ramp;
H.contour_var = contour_var;

guidata(hObject, H);

end

%% PUSHBUTTON Activate All
function Activate_all_Callback(hObject, eventdata, H)
opengl hardware
name = H.name;
N = H.N;
text = H.text;
data = H.data;
loaded_names = cellstr(get(H.Loaded_samples,'String'));

density_active = H.density_active;
set(H.run_nmf, 'visible','on');

%Get/set parameters
input_variables = {get(H.x_variable_name,'string') get(H.y_variable_name,'string')};
set(H.x_variable_name,'enable','off');
set(H.y_variable_name,'enable','off');
if size(density_active,1) == 0
%1D variables
    x_min = str2num(get(H.xaxis_min,'String')); 
    x_max = str2num(get(H.xaxis_max,'String'));
    x_int = str2num(get(H.xaxis_int,'String'));
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
%2D variables
    grid = 2^str2num(get(H.gridspc,'String'));
    lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
    lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));
    upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));
    upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));
    MIN_XY=[lower_Xlim,lower_Ylim];
    MAX_XY=[upper_Xlim,upper_Ylim];
    X_bandwidth = str2num(get(H.bandwidth_x,'String'));
    Y_bandwidth = str2num(get(H.bandwidth_y,'String'));
    color_ramp = get(H.color_chkbx,'Value');
    contour_var = get(H.contour_chkbx,'Value');
    contour_percentile = str2num(get(H.contour_percentile_bx,'String'));
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    set(H.color_chkbx,'enable','off');
    set(H.contour_chkbx,'enable','off');
    set(H.contour_percentile_bx,'enable','off');
else
    %1D
    x_min = H.x_min;
    x_max = H.x_max;
    x_int = H.x_int;
    x = x_min:x_int:x_max;
    set(H.xaxis_min,'enable','off');
    set(H.xaxis_max,'enable','off');
    set(H.xaxis_int,'enable','off');
    set(H.pdps,'enable','off');
    set(H.kdes,'enable','off');
    set(H.kern_opt,'enable','off');
    set(H.kern_set,'enable','off');
    set(H.kern,'enable','off');
    %2D
    grid = H.grid;
    lower_Xlim = H.lower_Xlim;
    lower_Ylim = H.lower_Ylim;
    upper_Xlim = H.upper_Xlim;
    upper_Ylim = H.upper_Ylim;
    MIN_XY     = [lower_Xlim,lower_Ylim];
    MAX_XY     = [upper_Xlim,upper_Ylim];
    X_bandwidth= H.X_bandwidth;
    Y_bandwidth= H.Y_bandwidth;
    contour_var = H.contour_var;
    color_ramp = H.color_ramp;
    contour_percentile = H.contour_percentile;
    set(H.xaxis_min_2D,'enable','off');
    set(H.xaxis_max_2D,'enable','off');
    set(H.bandwidth_x,'enable','off');
    set(H.yaxis_min_2D,'enable','off');
    set(H.yaxis_max_2D,'enable','off');
    set(H.bandwidth_y,'enable','off');
    set(H.gridspc,'enable','off');
    set(H.color_chkbx,'enable','off');
    set(H.contour_chkbx,'enable','off');
    set(H.contour_percentile_bx,'enable','off');
end
%end get/set parameters

if isempty(loaded_names) == 1
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(loaded_names)) == 0
	err_dlg=errordlg('There are no samples to activate!','Hang on...');
	waitfor(err_dlg);
else
    
%clear uipanels with 2D density plots, if they exist
try child = allchild(H.stacked_sink_samples_uipanel);
    delete (child); 
catch; 
end;

set(H.Active_samples, 'String', name);
set(H.Loaded_samples, 'String', []);
set(H.Loaded_samples, 'Value',1);

active_names = cellstr(get(H.Active_samples,'String'));
active_names(all(cellfun('isempty',active_names),2),:) = [];


    for i = 1:N
	for j = 1:length(active_names)
		if strcmp(text(1,i*2-1), active_names(j,1)) == 1
			active_data_name{:,i*2-1} = active_names(j,1);
			active_data(:,i*2-1:i*2) = data(:,i*2-1:i*2);
		end
    end
    end
if get(H.one_dimensional, 'Value') == 1
    active_data(isnan(active_data))=0;
    active_data( :, all(~active_data,1) ) = [];
    active_data_name = active_data_name(~cellfun('isempty',active_data_name));

    if get(H.pdps,'Value') == 1
	for i = 1:length(active_data_name)
		m = nonzeros(active_data(:,i*2-1));
		m = m(isfinite(m(:,1)),:);
		s = nonzeros(active_data(:,i*2));
		s = s(isfinite(s(:,1)),:);
		f = zeros(length(m),length(x));
		for j = 1:length(m)
			f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
		end
		pdps_active(:,i) = ((sum(f, 1))/length(m)).';
    end
    end

    if get(H.kdes,'Value') == 1
	rad_on_kernel=get(H.Dist_Opts_1D,'selectedobject');
	switch rad_on_kernel
		case H.kern_opt
			xA = transpose(x);
			for i = 1:length(active_data_name)
				dist_data = nonzeros(active_data(:,i*2-1));
				[bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
				pdps_active(:,i) = transpose(interp1(xmesh1, kdeA, xA));
				clear dist_data kdeA
			end
		case H.kern_set
			for i = 1:length(active_data_name)
				dist_data = nonzeros(active_data(:,i*2-1));
				kernel = str2num(get(H.kern,'String'));
				kernel_dist_data(1:length(nonzeros(active_data(:,i*2-1))),1) = kernel;
				pdps_active(:,i) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
			end
    end
    end

    nsamples = length(pdps_active(1,:));
    
    %Plot 1D distributions
    if get(H.show_sinks, 'value') == 1
    plot_1D_distribution(pdps_active, x_min, x_max, x_int,...
            H.stacked_sink_samples_uipanel, input_variables);
    end
    
    input_distributions = pdps_active;
    density_active = [];
    X1 = [];
    Y1 = [];
    density_active_vector = [];
end
if get(H.two_dimensional, 'Value') == 1
    active_data_name = active_data_name(~cellfun('isempty',active_data_name));
    r=ones(1,size(active_data,1));
    c=ones(1, size(active_data, 2));
    active_data=mat2cell(active_data, r, c);
    for i = 1:size(r,2)
	for j = 1:size(c,2)
		if cellfun('isempty', active_data(i,j)) == 0
			if cellfun(@isnan, active_data(i,j)) == 1
				active_data(i,j) = {[]};
			end	
		end
    end
    end
    grid = 2^str2num(get(H.gridspc,'String'));
    lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
    lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));;
    upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));;
    upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));;
    MIN_XY=[lower_Xlim,lower_Ylim];
    MAX_XY=[upper_Xlim,upper_Ylim];
    for i = 1:length(active_data_name)
        data_now=cell2mat(active_data(:,2*i-1:2*i));
        [bandwidth_not_needed,sink_density(:,:,i),X1,Y1]=kde2d(data_now, grid, MIN_XY, MAX_XY, X_bandwidth, Y_bandwidth);
        sink_density(:,:,i) = sink_density(:,:,i)./sum(sum(sink_density(:,:,i)));
        density_vector_new(:,i)=reshape(sink_density(:,:,i),1,[]);
    end
    density_active=sink_density;
    density_active_vector = density_vector_new;
    if get(H.show_sinks, 'value') == 1
        
        plot_2D_distribution(density_active, X1, Y1, ...
        H.stacked_sink_samples_uipanel, input_variables, color_ramp, contour_var, contour_percentile);
    end
    input_distributions = density_active;
    pdps_active = [];
    nsamples = size(density_active,3);
end


%%%%%%%%%Mean Kuiper
%update hObject
H.X1 = X1;
H.Y1 = Y1;
H.input_variables = input_variables;
guidata(hObject, H);
%end update hObject
[trash, V, mean_V_sink, range_V_sink,trash1, D, mean_D_sink, range_D_sink] =...
    KS_from_density(input_distributions, input_distributions, get(H.one_dimensional,'Value'), hObject) %calculate mean Kuiper based on function below
if get(H.one_dimensional, 'Value') == 1 % Use V values for 1 dimensional 
    set(H.mean_V_sink, 'String', mean_V_sink);
    set(H.range_V_sink, 'String', range_V_sink);
end
if get(H.two_dimensional, 'Value') == 1 % Use D values for 2 dimensional
    set(H.mean_V_sink, 'String', mean_D_sink);
    set(H.range_V_sink, 'String', range_D_sink);
end


%Mean Cross-correlation
if nsamples==1;
    set(H.mean_R2_sink, 'String', 'N/A');
else
    if get(H.one_dimensional, 'Value') == 1
        for i=1:nsamples
        for j=1:nsamples
            mean_R2_sink(i,j)=r2(pdps_active(:,i),pdps_active(:,j));
        end
        end
    end
    if get(H.two_dimensional, 'Value') == 1
        for i=1:nsamples
            density_active_vector(:,i)=reshape(density_active(:,:,i),1,[]);
        end
        
        for i=1:nsamples
        for j=1:nsamples
            
            mean_R2_sink(i,j)=r2(density_active_vector(:,i),density_active_vector(:,j));
        end
        end
    end
    mean_R2_sink_temp=[];
    for i=1:nsamples
        mean_R2_sink_temp=vertcat(mean_R2_sink_temp,diag(mean_R2_sink,i));
    end
    mean_R2_sink=round(mean(mean_R2_sink_temp),3);
    range_R2_sink=round(max(mean_R2_sink_temp)-min(mean_R2_sink_temp),3);
    set(H.mean_R2_sink, 'String', mean_R2_sink);
    set(H.range_R2_sink, 'String', range_R2_sink);
end

%Pass variables to GUI
set(H.nsamples, 'String', nsamples);
set(H.num_sources_opt,'String',round(nsamples*2/3));
active_selected = {1};
loaded_selected = {1};
H.active_selected = active_selected;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;
H.density_active = density_active;
H.density_active_vector = density_active_vector;
H.X1 = X1;
H.Y1 = Y1;
H.x_min = x_min;
H.x_max = x_max;
H.x_int = x_int;
H.grid = grid;
H.lower_Xlim = lower_Xlim;
H.lower_Ylim = lower_Ylim;
H.upper_Xlim = upper_Xlim;
H.upper_Ylim = upper_Ylim;
H.X_bandwidth= X_bandwidth;
H.Y_bandwidth= Y_bandwidth;
H.contour_percentile = contour_percentile;
H.color_ramp = color_ramp;
H.contour_var = contour_var;
%end pass variables to GUI

guidata(hObject, H);

end

%% PUSHBUTTON Deactivate All.
function Deactivate_all_Callback(hObject, eventdata, H)
set(H.x_variable_name,'enable','on');
set(H.y_variable_name,'enable','on');
%1D
set(H.xaxis_min,'enable','on');
set(H.xaxis_max,'enable','on');
set(H.xaxis_int,'enable','on');
set(H.pdps,'enable','on');
set(H.kdes,'enable','on');
set(H.kern_opt,'enable','on');
set(H.kern_set,'enable','on');
set(H.kern,'enable','on');
%2D
set(H.xaxis_min_2D,'enable','on');
set(H.xaxis_max_2D,'enable','on');
set(H.bandwidth_x,'enable','on');
set(H.yaxis_min_2D,'enable','on');
set(H.yaxis_max_2D,'enable','on');
set(H.bandwidth_y,'enable','on');
set(H.gridspc,'enable','on');
set(H.color_chkbx,'enable','on');
set(H.contour_chkbx,'enable','on');
set(H.contour_percentile_bx,'enable','on');

try child = allchild(H.stacked_sink_samples_uipanel)
    delete (child); 
catch; 
end;

try child = allchild(H.stacked_source_samples_uipanel)
    delete (child); 
catch; 
end;

active_names = cellstr(get(H.Active_samples,'String'));

if isempty(active_names) == 1
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
elseif length(char(active_names)) == 0
	err_dlg=errordlg('There are no samples to deactivate!','Hang on...');
	waitfor(err_dlg);
else

pdps_active = [];                   %set pdp data set as empty
density_active = [];                %set density data set as empty
density_active_vector = [];
name = H.name;

set(H.Active_samples, 'String', []);
set(H.Loaded_samples, 'String', name);
set(H.Active_samples, 'Value',1);

set(H.nsamples, 'String', 'N/A');
set(H.mean_R2_sink, 'String', 'N/A');
set(H.range_R2_sink, 'String', 'N/A');
set(H.mean_V_sink, 'String', 'N/A');
set(H.range_V_sink, 'String', 'N/A');
set(H.num_sources_opt,'String',20);

active_selected = {1};
loaded_selected = {1};
active_names = {};
H.active_selected = active_selected;
H.active_names = active_names;
H.loaded_selected = loaded_selected;
H.pdps_active = pdps_active;        %export empty pdp data set
H.density_active = density_active;  %export empty density data set
H.density_active_vector = density_active_vector;
guidata(hObject, H);

end

%% PUSHBUTTON Replot Results %%
function replot_res_Callback(hObject, eventdata, H)

optimized = H.optimized;

if 	optimized == 0
	err_dlg=errordlg('You need to determine the optimal number of sources or upload previous results.','Wait a sec...');
	waitfor(err_dlg);
elseif optimized == 1 && get(H.check_nmf_opt, 'Value') == 1
	xdata = H.xdata;
	ydata = H.ydata;
	SSR = H.SSR;
	modelpred = H.modelpred;

	figure;
	plot(xdata,ydata,'o',sort(xdata),modelpred,'r-');
	title('Final Residual versus Rank')
	xlabel('Rank')
	ylabel('Final Residual')

	figure;
	plot(xdata,SSR,'-o');
	title('SSR versus Rank')
	xlabel('Rank')
	ylabel('SSR1 + SSR2')
	
	optimized = 1;
	H.optimized = optimized;
	H.xdata = xdata;
	H.ydata = ydata;
	H.SSR = SSR;
	guidata(hObject, H);
end


%% PUSHBUTTON Export Plots %%
function export_plots_Callback(hObject, eventdata, H)
%get/set variables
density_active=H.density_active;
optimized = H.optimized;
nsources = H.nsources;
pdps_active = H.pdps_active;
source_PDP = H.source_PDP;
input_variables = H.input_variables;
factorized_density=H.factorized_density;
factorized_density_optimum = H.factorized_density_optimum;
X1=H.X1;
Y1=H.Y1;
lower_Xlim=str2num(get(H.xaxis_min_2D,'String'));
lower_Ylim=str2num(get(H.yaxis_min_2D,'String'));
upper_Xlim=str2num(get(H.xaxis_max_2D,'String'));
upper_Ylim=str2num(get(H.yaxis_max_2D,'String'));
x_min = str2num(get(H.xaxis_min,'String')); 
x_max = str2num(get(H.xaxis_max,'String'));
x_int = str2num(get(H.xaxis_int,'String'));
x = x_min:x_int:x_max;
contour_var = H.contour_var;
color_ramp = H.color_ramp;
contour_percentile = H.contour_percentile;


%open target figures
finput = figure('units','normalized','position', [0.1 0.1 0.5 0.8], 'MenuBar', 'figure', 'ToolBar', 'auto');

foutput = figure('units','normalized','position', [0.3 0.1 0.5 0.8]);

%plot if 1D distributions
if get(H.one_dimensional, 'Value') == 1
plot_1D_distribution(pdps_active, x_min, x_max, x_int, finput, input_variables);
set(finput, 'MenuBar', 'figure');
plot_1D_distribution(source_PDP, x_min, x_max, x_int, foutput, input_variables);

end

%plot if 2D distributions
if get(H.two_dimensional, 'Value') == 1
    plot_2D_distribution(density_active, X1, Y1, finput, input_variables, color_ramp, contour_var, contour_percentile);
    if get(H.check_nmf, 'Value') == 1
    plot_2D_distribution(factorized_density, X1, Y1, foutput, input_variables, color_ramp, contour_var, contour_percentile);
    end
    if get(H.check_nmf_opt, 'Value') == 1
    plot_2D_distribution(factorized_density_optimum, X1, Y1, foutput, input_variables, color_ramp, contour_var, contour_percentile);    
    end    
end

  

if optimized == 1 

	xdata = H.xdata;
	ydata = H.ydata;
	SSR = H.SSR;
	modelpred = H.modelpred;

	figure;
	plot(xdata,ydata,'o',sort(xdata),modelpred,'r-');
	title('Final Residual versus Rank')
	xlabel('Rank')
	ylabel('Final Residual')

	figure;
	plot(xdata,SSR,'-o');
	title('SSR versus Rank')
	xlabel('Rank')
	ylabel('SSR1 + SSR2')
	
	optimized = 1;

	H.xdata = xdata;
	H.ydata = ydata;
	H.SSR = SSR;
	
	
end

H.optimized = optimized;
guidata(hObject, H);

%% PUHBUTTON Example Input %%
function pushbutton20_Callback(hObject, eventdata, H)

Example_Data_Set_1D
Example_Data_Set_2D


%{
%% Calculate Kuiper from densities
function [V_1_to_1, mean_V_sink, range_V_sink] = KS_from_density(input_distribution1, input_distribution2, dimensions)

%Get/set variables
[X1,Y1]=meshgrid(1:(4000-1)/511:4000,-1:(1-(-1))/511:1);
input_variables = {'Age', 'U/Th'};

if length(size(input_distribution1)) ~= length(size(input_distribution2))
    err_dlg=errordlg('Input distributions must have the same dimensions');
	waitfor(err_dlg);
end

    

if length(size(input_distribution1))==1;
    
    mean_V_sink = {'N/A'};
    range_V_sink = {'N/A'};
end
if length(size(input_distribution1))==2 || dimensions == 1
        pdps_active1=input_distribution1;
        pdps_active2=input_distribution2;
        nsamples = size(pdps_active1, 2);
    
        for i=1:nsamples
        for j=1:nsamples
            V(i,j)=max(cumsum(pdps_active1(:,i))-...
                cumsum(pdps_active2(:,j)))+max(cumsum(pdps_active2(:,j))...
                -cumsum(pdps_active1(:,i)));
        end
        end
        mean_V_sink_temp=[];
        V_1_to_1 = diag(V);
        for i=1:nsamples
             mean_V_sink_temp=vertcat(mean_V_sink_temp,diag(V,i));
        end
        mean_V_sink=round(mean(mean_V_sink_temp),3);
        range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
end
if length(size(input_distribution1))==3 || dimensions == 0
        density_active1 = input_distribution1;
        density_active2 = input_distribution2;
        nsamples = size(input_distribution1, 3);
        for i=1:size(density_active1,3)
            CDFx1(:,:,i) = cumsum(density_active1(:,:,i),1);    % take the CDF of x at y = max (y)
            density_active_CDF1(:,:,i) = cumsum(CDFx1(:,:,i),2);   
            CDFx2(:,:,i) = cumsum(density_active2(:,:,i),1);    % take the CDF of x at y = max (y)
            density_active_CDF2(:,:,i) = cumsum(CDFx2(:,:,i),2);  
        end
        
        
%enable for testing        plot_2D_distribution(density_active_CDF1, X1, Y1, figure, input_variables);
%enable for testing        plot_2D_distribution(density_active_CDF2, X1, Y1, figure, input_variables);
        for i=1:nsamples
        for j=1:nsamples
            V(i,j) = max(max(density_active_CDF1(:,:,i)-density_active_CDF2(:,:,j))...
                    +max(density_active_CDF2(:,:,j)-density_active_CDF1(:,:,i)));
        end
        end

    mean_V_sink_temp=[];
    for i=1:nsamples
        mean_V_sink_temp=vertcat(mean_V_sink_temp,diag(V,i));
    end
    mean_V_sink=round(mean(mean_V_sink_temp),3);
    range_V_sink=round(max(mean_V_sink_temp)-min(mean_V_sink_temp),3);
    
end
V_1_to_1 = diag(V);

if length(size(input_distribution1))>3; 
    err_dlg=errordlg('Error in calculating mean Kuiper');
	waitfor(err_dlg);    
end


%%    PLOT 2D Distributions %%%
function [] = plot_2D_distribution ...
    (input_densities, X1, Y1, input_uipanel, input_variables);
%clear uipanels if they exist
try child = allchild(input_uipanel);
    delete (child); 
catch; 
end;

density_active = input_densities;
nsamples = size(density_active,3);
%parent = strcat('H.',input_uipanel);
xtics = min(X1(1,:)):500:max(X1(1,:));
%set(input_uipanel,'Visible', 'off');
f = waitbar(0, 'Plotting 2D densities', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0)
for i=1:nsamples
            left = 0.075;
            bottom = 0.075 + (((nsamples-i)/nsamples)*0.90);
            width = 0.9;
            height = 0.85/nsamples;
            pos = [left bottom width height];
            
            ax(i)=axes(input_uipanel, 'Position', pos);
            
            colormap(jet);
            surf(ax(i), X1,Y1,density_active(:,:,i), 'EdgeColor','none');
            view(ax(i),[0, 90]);
            
            ax(i).TickDir = 'out';
            ax(i).TickLength = [0.01 0.02];
            ax(i).XTickLabel = {};
            %axis(ax(i),[lower_Xlim  upper_Xlim  lower_Ylim  upper_Ylim]);
            
            
            if getappdata(f,'canceling')
                break
            end
            waitbar(i/nsamples, f);
end
        delete(f);
        f=msgbox('Scaling and labeling','Please wait');
            xticklabels(ax(nsamples),'auto');
            xlabel(ax(nsamples), string(input_variables(1,1)),'Visible','on');
            yticklabels('auto')
            plot_top = ax(1);
            plot_bottom = ax(nsamples);
            p1=get(plot_top,'position');
            p2=get(plot_bottom,'position');
            height = -p2(2)+p1(2)+p1(4);
        
            h3=axes(input_uipanel, 'position',[p2(1) p2(2) p2(3) height],'visible','off');
        
            ylab = ylabel(h3, string(input_variables(1,2)), 'Visible','on');
        
        delete(f);
        set(input_uipanel,'Visible', 'on');
%}
        
        
%%    PLOT 1D Distributions %%%
function [] = plot_1D_distribution ...
    (pdps_active, x_min, x_max, x_int, input_uipanel, input_variables);
%clear uipanels if they exist
try child = allchild(input_uipanel);
    clear map; 
catch; 
end;
x = x_min:x_int:x_max;
pdps_active = pdps_active;
nsamples = size(pdps_active,2);
%parent = strcat('H.',input_uipanel);
xtics = x_min:500:x_max;

f = waitbar(0, 'Plotting 2D densities', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f, 'canceling', 0)
for i=1:nsamples
    
            left = 0.075;
            bottom = 0.075 + (((nsamples-i)/nsamples)*0.90);
            width = 0.9;
            height = 0.85/nsamples;
            pos = [left bottom width height];
            ax(i)=axes(input_uipanel, 'Position', pos);
            area(ax(i),x,pdps_active(:,i), 'Facecolor','k','linewidth',1);
            view(ax(i),[0, 90]);
            ax(i).TickDir = 'out';
            ax(i).TickLength = [0.01 0.02];
            ax(i).XTickLabel = {};
            set(ax(i), 'box','off');
            yticks(ax(i),[]);
            ax(i).YTickLabel = {};
            if getappdata(f,'canceling')
                break
            end
            waitbar(i/nsamples, f);
end
        delete(f);
        f=msgbox('Scaling and labeling','Please wait');
            xticklabels(ax(nsamples),'auto');
            xlabel(ax(nsamples), string(input_variables(1,1)),'Visible','on');
            ax(nsamples).YTickLabel = {};
            plot_top = ax(1);
            plot_bottom = ax(nsamples);
            p1=get(plot_top,'position');
            p2=get(plot_bottom,'position');
            height = -p2(2)+p1(2)+p1(4);
        
            h3=axes(input_uipanel, 'position',[p2(1) p2(2) p2(3) height],'visible','off');
        
            ylab = ylabel(h3, 'Normalized Probability' , 'Visible','on');
        
        delete(f);
        
                
        
% --- Executes on button press in compare_sources.
function compare_sources_Callback(hObject, eventdata, H)
%compare_sources
%Get/set variables
grid = 2^str2num(get(H.gridspc,'String'));
if get(H.one_dimensional, 'Value') == 1 
    source_PDP = H.source_PDP;
end
if get(H.two_dimensional, 'Value') == 1 && get(H.check_nmf, 'Value') == 1
    factorized_density = H.factorized_density;
end
if get(H.two_dimensional, 'Value') == 1 && get(H.check_nmf_opt, 'Value') == 1
    factorized_density = H.factorized_density_optimum;
end
grid = 2^str2num(get(H.gridspc,'String'));
lower_Xlim = str2num(get(H.xaxis_min_2D,'String'));
lower_Ylim = str2num(get(H.yaxis_min_2D,'String'));
upper_Xlim = str2num(get(H.xaxis_max_2D,'String'));
upper_Ylim = str2num(get(H.yaxis_max_2D,'String'));
MIN_XY = [lower_Xlim,lower_Ylim];
MAX_XY = [upper_Xlim,upper_Ylim];
X_bandwidth = str2num(get(H.bandwidth_x,'String'));
Y_bandwidth = str2num(get(H.bandwidth_y,'String'));
X1=H.X1;
Y1=H.Y1;
x_min = H.x_min;
x_max = H.x_max;
x_int = H.x_int;
input_variables = H.input_variables;
color_ramp = H.color_ramp;
contour_var = H.contour_var;
contour_percentile = H.contour_percentile;
x = x_min:x_int:x_max;

[filename, pathname] = uigetfile({'*'},'File Selector');
fullpathname = strcat(pathname, filename);
[numbers, text, data_tmp] = xlsread(fullpathname);
data = cell2mat(data_tmp(2:end,:));
[dataR,dataC]=size(data);
N = (dataC/2); % Number of source samples
for i = 1:N
		name(i,1) = text(1,i*2-1);
        if size(text, 1)<1;
            err_dlg=errordlg('Please add headers to your input data. Press Example Input button to see an example.','Error!');
            waitfor(err_dlg);
        end
        
end
    
    
%%%%Format data and remove zeros and NaN
new_names = name;
data = cell2mat(data_tmp(3:end,:));
new_data = data;

nonzero = 1;
for i=1:size(new_data,2)
    if nnz(new_data(:,i))> 0
        temp(:,nonzero)=new_data(:,i);
        nonzero = nonzero + 1;
    end
end
new_data = temp; 

r=ones(1,size(new_data,1));
c=ones(1, size(new_data, 2));
new_data=mat2cell(new_data, r, c);
for i = 1:size(r,2)
	for j = 1:size(c,2)
		if cellfun('isempty', new_data(i,j)) == 0
			if cellfun(@isnan, new_data(i,j)) == 1
				new_data(i,j) = {[]};
			end	
		end
	end
end
%%%%end Format data and remove zeros and NaN

%%%Make density from empirical sources
if get(H.one_dimensional, 'Value') == 1 % for 1D distributions
    %Make PDPs          
    if get(H.pdps,'Value') == 1
        for i = 1:length(new_names)
            temp = cell2mat(new_data(:,i*2-1:i*2));
            
            m = nonzeros(temp(:,1));
    		m = m(isfinite(m(:,1)),:);
        	s = nonzeros(temp(:,2));
            s = s(isfinite(s(:,1)),:);
            if length(m)~=length(s)
                err_dlg=errordlg('Zeros in PDP bandwidth. Use KDEs instead');
                waitfor(err_dlg);
                return 
            end
            f = zeros(length(m),length(x));
            for j = 1:length(m)
            	f(j,:) = (1./ (s(j)*sqrt(2*pi)) .* exp (  (-((x-m(j)).^2)) ./ (2*((s(j)).^2))  ).*x_int);
            end
            pdps_new(:,i) = ((sum(f, 1))/length(m)).';
        end
    end
    
    %Make KDEs
    if get(H.kdes,'Value') == 1
        rad_on_kernel=get(H.Dist_Opts_1D,'selectedobject');
        switch rad_on_kernel
            case H.kern_opt
                xA = transpose(x);
                for i = 1:length(new_names)
                    dist_data = nonzeros(cell2mat(new_data(:,i*2-1)));
                    [bandwidth,kdeA,xmesh1,cdf]=kde(dist_data,length(x),x_min,x_max);
                    pdps_new(:,i) = transpose(interp1(xmesh1, kdeA, xA));
                    clear dist_data kdeA
                end
            case H.kern_set
                for i = 1:length(new_names)
                    dist_data = nonzeros(cell2mat(new_data(:,i*2-1)));
                    kernel = str2num(get(H.kern,'String'));
                    kernel_dist_data(1:length(nonzeros(cell2mat(new_data(:,i*2-1)))),1) = kernel;
                    pdps_new(:,i) = pdp5(dist_data(:,1),kernel_dist_data,x_min,x_max,x_int); 
                end
        end
    end
    empirical_pdps = pdps_new;
    density_active=[];
    input_distributions = empirical_pdps;
end
if get(H.two_dimensional, 'Value') == 1 % for 2D distributions
    for i = 1:length(name)
        data_now=cell2mat(new_data(:,2*i-1:2*i));
        [bandwidth_not_needed,empirical_sources(:,:,i),X1,Y1]=kde2d(data_now, grid, MIN_XY, MAX_XY, X_bandwidth, Y_bandwidth);
        empirical_sources(:,:,i) = empirical_sources(:,:,i)./sum(sum(empirical_sources(:,:,i)));
        density_vector_new(:,i)=reshape(empirical_sources(:,:,i),1,[]);
    end
end

%plot empirical source densities
empirical_figure = figure;
set(empirical_figure,'WindowStyle','normal', 'MenuBar','figure');
if get(H.one_dimensional, 'Value') == 1 % for 1D distributions
    nsamples = length(empirical_pdps(1,:));
    plot_1D_distribution(empirical_pdps, x_min, x_max, x_int,...
    empirical_figure, input_variables);
    X1=NaN;                     %set 2D elements as NaN
    Y1=NaN;                     %set 2D elements as NaN         
end
if get(H.two_dimensional, 'Value') == 1 % for 2D distributions
    plot_2D_distribution(empirical_sources, X1, Y1, ...
            empirical_figure, input_variables, color_ramp, contour_var, contour_percentile);
end
set(empirical_figure, 'NumberTitle', 'off', 'Name', 'Empirical Potential Sources');   

%plot factorized densities
factorized_figure = figure
if get(H.one_dimensional, 'Value') == 1 % for 1D distributions
    plot_1D_distribution(source_PDP, x_min, x_max, x_int,...
    factorized_figure, input_variables);
end
if get(H.two_dimensional, 'Value') == 1 % for 2D distributions
    plot_2D_distribution(factorized_density, X1, Y1, ...
            factorized_figure, input_variables, color_ramp, contour_var, contour_percentile);
end
set(factorized_figure, 'NumberTitle', 'off', 'Name', 'Factorized Sources', 'MenuBar','figure');   

%compare factorized and empirical sources with cross-correlation
%coefficient
if get(H.one_dimensional, 'Value') == 1 % set up distributions for 1D distributions
    factorized_density_vector = source_PDP;
    sink_density_vector = empirical_pdps;
    number_of_empirical = size(empirical_pdps,2);
    number_of_factorized = size(source_PDP,2);
end
if get(H.two_dimensional, 'Value') == 1 % set up distributions for 2D distributions
    number_of_empirical=size(empirical_sources,3);
    for i=1:number_of_empirical
        sink_density_vector(:,i)=reshape(empirical_sources(:,:,i),[],1);
    end
    number_of_factorized=size(factorized_density,3)
    for i=1:number_of_factorized;
        factorized_density_vector(:,i) = reshape(factorized_density(:,:,i),[],1);
    end
end
for i=1:number_of_empirical % run cross-correlation
    for j = 1:number_of_factorized
        R2(i,j)=r2(factorized_density_vector(:,j), sink_density_vector(:,i));
    end
end

[B, ind] = sort(R2);
compare_sources = [];
for i=1:number_of_factorized
    if get(H.one_dimensional, 'Value') == 1
        compare_sources(:,2*i-1) = source_PDP(:,i);
        compare_sources(:,2*i) = empirical_pdps(:,ind(number_of_empirical, i));
        factorized_density = source_PDP; %set variables for KS test
        empirical_sources = empirical_pdps; %set variables for KS test
    end
    if get(H.two_dimensional, 'Value') == 1
        compare_sources(:,:,2*i-1) = factorized_density(:,:,i);
        compare_sources(:,:,2*i) = empirical_sources(:,:,ind(number_of_empirical, i));
    end
end

for i=1:number_of_factorized
    vars(i) = cellstr(strcat('Sources',' ',num2str(i)));
end
cross_correlation_comparison=figure
%f.Units = 'normalized'
cross_correlation_comparison.Units = 'normalized'
%f.Position = [0.3536 0.5157 0.6 0.37]
cross_correlation_comparison.Position = [0.3536 0.5157 0.6 0.37]

table_data_R2 = horzcat(new_names, num2cell(R2));
header = horzcat('Empirical sources', vars);
table_data_R2  = vertcat(header, table_data_R2);
uit = uitable
uit.Units = 'normalized';
uit.Position = [0 0 1 1];
%uit.ColumnName = vars; 
uit.Data = table_data_R2;
%uit.RowName = new_names;
stop = 1
set(cross_correlation_comparison, 'NumberTitle', 'off', 'Name', 'Cross-correlation of empirical and factorize sources');

%Calculate KS D value comparison
[V, V_1_to_1, mean_V_sink, range_V_sink,D, D_1_to_1, mean_D_sink, range_D_sink] = KS_from_density...
    (empirical_sources, factorized_density, get(H.one_dimensional,'Value'),hObject, Y1, input_variables)

KS_comparison=figure
KS_comparison.Units = 'normalized'
KS_comparison.Position = [0.3536 0.5157 0.6 0.37]

table_data_KS = horzcat(new_names, num2cell(D));
%header = horzcat('Empirical sources', vars);
table_data_KS  = vertcat(header, table_data_KS);
uit = uitable
uit.Units = 'normalized';
uit.Position = [0 0 1 1];
%uit.ColumnName = vars; 
uit.Data = table_data_KS;
%uit.RowName = new_names;
stop = 1
set(KS_comparison, 'NumberTitle', 'off', 'Name', 'Kolmogorov-Smirnov D value of empirical and factorize sources');

comparison_figure=figure;
if get(H.one_dimensional, 'Value') == 1 % for 1D distributions
    plot_1D_distribution(compare_sources, x_min, x_max, x_int,...
    comparison_figure, input_variables);
end
if get(H.two_dimensional, 'Value') == 1 % for 2D distributions
    plot_2D_distribution(compare_sources, X1, Y1, ...
            comparison_figure, input_variables, color_ramp, contour_var, contour_percentile);
end
set(comparison_figure, 'NumberTitle', 'off', 'Name', 'Comparison of factorized sources and their closest empirical matches based on Cross-correlation');
guidata(hObject, H);



























































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








































%% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


%% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)



%% PUSHBUTTON Export NMF Results
%% --- Executes on button press in Export_NMF_Results.
function Export_NMF_Results_Callback(hObject, eventdata, H)
global iter_num;
optimized = H.optimized;
        export = H.export;
    
    
        export_2D = H.export_2D;
    
nsources = H.nsources;
export = H.export


if get(H.check_nmf, 'Value') == 1
	BaseName=' Sources';
	tst3{1,1}=[num2str(nsources),BaseName];
	[file,path] = uiputfile('*.xlsx','Save file');
    f=msgbox('Writing file...This box will close when writing is complete','Please wait');
    warning('off','MATLAB:xlswrite:AddSheet');
	writecell(export,[path file], 'Sheet', char(tst3(1,1)));
    if get(H.two_dimensional, 'Value') == 1 && get(H.check_nmf, 'Value') == 1 
        density_file = strcat(file, 'densities');
        rows = size(export_2D,2);
        export_temp = cell(size(export_2D,3)*(size(export_2D,1)+1), size(export_2D, 2));
        for i=1:size(export_2D,3)
           name = strcat('Density ',' ', num2str(i));
           export_temp(i+rows*(i-1),1) = cellstr(name); 
           export_temp(i+1+rows*(i-1):i*rows+i,:) = num2cell(export_2D(:,:,i)); 
        end
        writecell(export_temp, [path file],'Sheet', '2D Densities'); 
    end
        
	%RemoveSheet123([path file]);
	%warning('off');
end

if get(H.check_nmf_opt, 'Value') == 1
    breakpoint_output = H.breakpoint_output;
    
	tst3{1,iter_num-1} = [];
	BaseName=' Sources';
	for i=2:iter_num
		tst3{1,i-1}=[num2str(i),BaseName];
	end
	[file,path] = uiputfile('*.xlsx','Save file');
    f=msgbox('Writing results...This box will close when writing is complete','Please wait');
	warning('off','MATLAB:xlswrite:AddSheet');
    
    %write summary results sheet
    writecell(breakpoint_output,[path file], 'Sheet', 'BreakpointAnalysis');
    try
        BCV_analysis = H.BCV_analysis;
        
        writecell(BCV_analysis,[path,file],'Sheet','BCVanalysis');
    catch
    end
	for i=1:iter_num-1
		writecell(export(:,:,i),[path file], 'Sheet', char(tst3(1,i)));
    end
        
    %if 2D, write densities to seperate sheet
    g = questdlg('Write distributions to file? This may take a long time.','Warning','Yes','No','Yes')
    if strcmp(g,'Yes')
    if get(H.two_dimensional, 'Value') == 1 && get(H.check_nmf_opt, 'Value') == 1 
        delete(f)
        f=uifigure;
        boxmsg = strcat('Writing ', {' '}, string(1), {' '}, ' of ', {' '}, ...
                string(iter_num-1), {' '}, ' distributions. This message will close when writing is complete. Select cancel to end now')
        g = uiprogressdlg(f, 'Title','Saving', 'Message',boxmsg, 'Indeterminate','on','Cancelable','on')
        for i=2:iter_num
            g.Message = strcat('Writing ', {' '}, string(i-1), {' '}, ' of ', {' '}, ...
                string(iter_num-1), {' '}, ' distributions. This message will close when writing is complete. Select cancel to end now')
            if g.CancelRequested
                break
            end
            tst3{1,i-1}=['Densities for ',num2str(i),' Sources'];
            writecell(export_2D(:,:,i-1),[path file], 'Sheet', char(tst3(1,i-1)));
        end
        close(g);
    end
    end
end
try delete (f);
catch
end
	%RemoveSheet123([path file]);
	%warning('off');








	
	
	
	

function iterations_Callback(hObject, eventdata, H)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterations as text
%        str2double(get(hObject,'String')) returns contents of iterations as a double


%% --- Executes during object creation, after setting all properties.
function iterations_CreateFcn(hObject, eventdata, H)
% hObject    handle to iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



































%% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


%% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)





%% --- Executes during object creation, after setting all properties.
function Active_samples_CreateFcn(hObject, eventdata, H)
% hObject    handle to Active_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









function num_sources_Callback(hObject, eventdata, H)
% hObject    handle to num_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_sources as text
%        str2double(get(hObject,'String')) returns contents of num_sources as a double


%% --- Executes during object creation, after setting all properties.
function num_sources_CreateFcn(hObject, eventdata, H)
% hObject    handle to num_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Final_Residual_Callback(hObject, eventdata, H)
% hObject    handle to Final_Residual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Final_Residual as text
%        str2double(get(hObject,'String')) returns contents of Final_Residual as a double


%% --- Executes during object creation, after setting all properties.
function Final_Residual_CreateFcn(hObject, eventdata, H)
% hObject    handle to Final_Residual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tof_Callback(hObject, eventdata, H)
% hObject    handle to tof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tof as text
%        str2double(get(hObject,'String')) returns contents of tof as a double


%% --- Executes during object creation, after setting all properties.
function tof_CreateFcn(hObject, eventdata, H)
% hObject    handle to tof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, H)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, H)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)












function num_sources_opt_Callback(hObject, eventdata, H)
% hObject    handle to num_sources_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_sources_opt as text
%        str2double(get(hObject,'String')) returns contents of num_sources_opt as a double


% --- Executes during object creation, after setting all properties.
function num_sources_opt_CreateFcn(hObject, eventdata, H)
% hObject    handle to num_sources_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, H)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)






function xaxis_min_2D_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_min_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_min_2D as text
%        str2double(get(hObject,'String')) returns contents of xaxis_min_2D as a double


% --- Executes during object creation, after setting all properties.
function xaxis_min_2D_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_min_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xaxis_max_2D_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_max_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_max_2D as text
%        str2double(get(hObject,'String')) returns contents of xaxis_max_2D as a double


% --- Executes during object creation, after setting all properties.
function xaxis_max_2D_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_max_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xaxis_max_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_max_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_max_2D as text
%        str2double(get(hObject,'String')) returns contents of xaxis_max_2D as a double


% --- Executes during object creation, after setting all properties.
function xaxis_max_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_max_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xaxis_min_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_min_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_min_2D as text
%        str2double(get(hObject,'String')) returns contents of xaxis_min_2D as a double


% --- Executes during object creation, after setting all properties.
function xaxis_min_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_min_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xaxis_int_Callback(hObject, eventdata, H)
% hObject    handle to xaxis_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxis_int as text
%        str2double(get(hObject,'String')) returns contents of xaxis_int as a double


% --- Executes during object creation, after setting all properties.
function xaxis_int_CreateFcn(hObject, eventdata, H)
% hObject    handle to xaxis_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kern_opt.
function kern_opt_Callback(hObject, eventdata, H)
% hObject    handle to kern_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kern_opt


% --- Executes on button press in kern_set.
function kern_set_Callback(hObject, eventdata, H)
% hObject    handle to kern_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kern_set



function kern_Callback(hObject, eventdata, H)
% hObject    handle to kern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    structure with H and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kern as text
%        str2double(get(hObject,'String')) returns contents of kern as a double


% --- Executes during object creation, after setting all properties.
function kern_CreateFcn(hObject, eventdata, H)
% hObject    handle to kern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% H    empty - H not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function opt_num_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opt_num_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





% --- Executes on button press in check_nmf_opt.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to check_nmf_opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_nmf_opt


% --- Executes on button press in check_nmf.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to check_nmf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_nmf


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% PUSHBUTTON Cancel %%
function cancel_Callback(hObject, eventdata, handles)

global cancel
cancel = 1;



function bandwidth_x_Callback(hObject, eventdata, handles)
% hObject    handle to bandwidth_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bandwidth_x as text
%        str2double(get(hObject,'String')) returns contents of bandwidth_x as a double


% --- Executes during object creation, after setting all properties.
function bandwidth_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bandwidth_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yaxis_min_2D_Callback(hObject, eventdata, handles)
% hObject    handle to yaxis_min_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaxis_min_2D as text
%        str2double(get(hObject,'String')) returns contents of yaxis_min_2D as a double


% --- Executes during object creation, after setting all properties.
function yaxis_min_2D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxis_min_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yaxis_max_2D_Callback(hObject, eventdata, handles)
% hObject    handle to yaxis_max_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaxis_max_2D as text
%        str2double(get(hObject,'String')) returns contents of yaxis_max_2D as a double


% --- Executes during object creation, after setting all properties.
function yaxis_max_2D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxis_max_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bandwidth_y_Callback(hObject, eventdata, handles)
% hObject    handle to bandwidth_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bandwidth_y as text
%        str2double(get(hObject,'String')) returns contents of bandwidth_y as a double


% --- Executes during object creation, after setting all properties.
function bandwidth_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bandwidth_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function gridspc_Callback(hObject, eventdata, handles)
% hObject    handle to gridspc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gridspc as text
%        str2double(get(hObject,'String')) returns contents of gridspc as a double


% --- Executes during object creation, after setting all properties.
function gridspc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gridspc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_variable_name_Callback(hObject, eventdata, handles)
% hObject    handle to x_variable_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_variable_name as text
%        str2double(get(hObject,'String')) returns contents of x_variable_name as a double


% --- Executes during object creation, after setting all properties.
function x_variable_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_variable_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_variable_name_Callback(hObject, eventdata, handles)
% hObject    handle to y_variable_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_variable_name as text
%        str2double(get(hObject,'String')) returns contents of y_variable_name as a double


% --- Executes during object creation, after setting all properties.
function y_variable_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_variable_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in color_chkbx.
function color_chkbx_Callback(hObject, eventdata, handles)
% hObject    handle to color_chkbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of color_chkbx


% --- Executes on button press in contour_chkbx.
function contour_chkbx_Callback(hObject, eventdata, handles)
% hObject    handle to contour_chkbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of contour_chkbx



function contour_percentile_bx_Callback(hObject, eventdata, handles)
% hObject    handle to contour_percentile_bx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contour_percentile_bx as text
%        str2double(get(hObject,'String')) returns contents of contour_percentile_bx as a double


% --- Executes during object creation, after setting all properties.
function contour_percentile_bx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contour_percentile_bx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
