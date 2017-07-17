function varargout = NaiveBayesClassifier(varargin)
% NAIVEBAYESCLASSIFIER M-file for NaiveBayesClassifier.fig
%      NAIVEBAYESCLASSIFIER, by itself, creates a new NAIVEBAYESCLASSIFIER or raises the existing
%      singleton*.
%
%      H = NAIVEBAYESCLASSIFIER returns the handle to a new NAIVEBAYESCLASSIFIER or the handle to
%      the existing singleton*.
%
%      NAIVEBAYESCLASSIFIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NAIVEBAYESCLASSIFIER.M with the given input arguments.
%
%      NAIVEBAYESCLASSIFIER('Property','Value',...) creates a new NAIVEBAYESCLASSIFIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NaiveBayesClassifier_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NaiveBayesClassifier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NaiveBayesClassifier

% Last Modified by GUIDE v2.5 09-May-2013 06:54:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NaiveBayesClassifier_OpeningFcn, ...
                   'gui_OutputFcn',  @NaiveBayesClassifier_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
    clc;
    
    % set window potition (get_size_screen/gsl_)
    gsl_ = get(0,'ScreenSize');
    
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NaiveBayesClassifier is made visible.
function NaiveBayesClassifier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NaiveBayesClassifier (see VARARGIN)

% Choose default command line output for NaiveBayesClassifier
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NaiveBayesClassifier wait for user response (see UIRESUME)
% uiwait(handles.NaiveBayesClassifier);

% Set the figure icon by matlabfreecode.wordpress.com
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.NaiveBayesClassifier,'javaframe');
jIcon=javax.swing.ImageIcon('citrus-icon.png');
jframe.setFigureIcon(jIcon);

% --- Outputs from this function are returned to the command line.
function varargout = NaiveBayesClassifier_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in NaiveBayesClassifier.
function trainingdata_Callback(hObject, eventdata, handles)
% hObject    handle to NaiveBayesClassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NaiveBayesClassifierProject=guidata(gcbo);
GetImageTraining=get(NaiveBayesClassifierProject.ImageTraining,'Userdata');

% determine the training data path
path_data_train=strrep(cd,...
    'Matlab_Code_To_Classification_Citrus','CitrusImage\Data Training');

% data training of Citrus nipis (jn)
lots_of_data_train_jn=7;

% data training of Citrus lemon (jl)
lots_of_data_train_jl=3;

% data training of Citrus orange (jm)
lots_of_data_train_jm=5;

lots_of_feature=4;
lots_of_class=3;

% initialization of matrix dataset
dataset=zeros(lots_of_data_train_jn,lots_of_feature);

for i=1:(lots_of_data_train_jn+lots_of_data_train_jl+lots_of_data_train_jm)
    
    if(i<=lots_of_data_train_jn)
        % membaca setiap file citra jeruk nipis
        filename=strcat(path_data_train,'\','CitrusNipis',...
        num2str(i),'.jpg');
        class{i}='Nipis';
    elseif(i<=(lots_of_data_train_jn+lots_of_data_train_jl))
        % membaca setiap file citra jeruk lemon
        filename=strcat(path_data_train,'\','CitrusLemon',...
            num2str(i-lots_of_data_train_jn),'.jpg');
        class{i}='Lemon';
    else
        % membaca setiap file citra jeruk manis
        filename=strcat(path_data_train,'\','CitrusOrange',...
            num2str(i-(lots_of_data_train_jn+lots_of_data_train_jl)),'.jpg');
        class{i}='Orange';
    end
    
    I = imread (filename);
    
    % merisize image size
    % I = imresize(I,[256 256]);

    if(size(I,3)==4)
        I(:,:,1)=[]; % convert to I = [MxNx3]
    end
    
    % count mean Red, Green, Blue
    mean_red=mean(mean(I(:,:,1)));
    mean_green=mean(mean(I(:,:,2)));
    mean_blue=mean(mean(I(:,:,3)));
       
    axes(handles.ImageTraining);
    imshow(I);    
    
    set(NaiveBayesClassifierProject.ImageTraining,'Userdata',I);   
    
    SizeImageTraining=size(I);
    StringSizeImageTraining=sprintf(strcat(num2str(SizeImageTraining(1)),'x',num2str(SizeImageTraining(2)),'x',num2str(SizeImageTraining(3))));    
    StringSizeImageTraining = strrep(StringSizeImageTraining,'x',' x ');
    
    set(handles.SizeImageTraining,'String',StringSizeImageTraining);
    
    % set(handles.SizeImageTraining,'String',StringSizeImageTraining);
    
    set(handles.NameImageTraining,'String',sprintf(strcat('Image Training (Color) ->',num2str(i))));
    
    %% create gray-scale image
    I_gray=Function_ColorToGray(I);
    
    axes(handles.ImageGray);
    imshow(I_gray);  
    
    %%
    
    %% create binary image
%     level=graythresh(I);
%     I_biner=im2bw(I_gray,level);
    
    I_biner=zeros(size(I_gray,1),size(I_gray,2));
    I_biner(find(I_gray<255))=1;

    axes(handles.ImageBinary);
    imshow(I_biner);  
    
    %% create max filter image from binary image
       %windowing_size must be valued an odd number >=3
       windowing_size=5;
       max_filter_I_biner=Function_MaxFilterBiner_(I_biner,windowing_size); 
       
       axes(handles.ImageMaxFilter);
       imshow(max_filter_I_biner);
       
       
    %% count diameter with unit length horizontal each pixel 
       % determine index which contains the value 1
         [idx_,idy_]=find(max_filter_I_biner==1);
         diameter=idy_(numel(idy_))-idy_(1)+1;
    
            
    %% create boundary line in diameter horizontal
        
        %numel((XY2Index(1,idy_(1),size(I_gray,1)):XY2Index(size(I_gray,1),idy_(numel(idy_)),size(I_gray,1)))')
        
        
        % replace pixel value
        I_red=I(:,:,1);
        I_red(:,idy_(1):(idy_(1)+10))=105;
        I_red(:,(idy_(numel(idy_))-10):idy_(numel(idy_)))=105;
        I(:,:,1)=I_red;
        
        I_green=I(:,:,2);
        I_green(:,idy_(1):(idy_(1)+10))=75;
        I_green(:,(idy_(numel(idy_))-10):idy_(numel(idy_)))=75;
        I(:,:,2)=I_green;
        
        I_blue=I(:,:,3);
        I_blue(:,idy_(1):(idy_(1)+10))=245; 
        I_blue(:,(idy_(numel(idy_))-10):idy_(numel(idy_)))=245;
        I(:,:,3)=I_blue;
        
        
        axes(handles.ImageTraining);
        imshow(I); 
        
%     %% count diameter with unit length vertical each pixel
%        % determine index which contains the value 1
%          [idx_v,idy_v]=find(max_filter_I_biner==1);         
%          [value_max,idx_v_max]=max(idx_v);
%          [value_min,idx_v_min]=min(idx_v);         
%          diameter_v=value_max-value_min+1;
%          
%     %% create boundary line in diameter vertical    
%         
%         % replace pixel value
%         I_red=I(:,:,1);
%         I_red(idx_v(idx_v_min):(idx_v(idx_v_min)+10),:)=105;
%         I_red((idx_v(idx_v_max)-10):idx_v(idx_v_max),:)=105;
%         I(:,:,1)=I_red;
%         
%         I_green=I(:,:,2);
%         I_green(idx_v(idx_v_min):(idx_v(idx_v_min)+10),:)=75;
%         I_green((idx_v(idx_v_max)-10):idx_v(idx_v_max),:)=75;
%         I(:,:,2)=I_green;
%         
%         I_blue=I(:,:,3);
%         I_blue(idx_v(idx_v_min):(idx_v(idx_v_min)+10),:)=245; 
%         I_blue((idx_v(idx_v_max)-10):idx_v(idx_v_max),:)=245;
%         I(:,:,3)=I_blue;        
%         
%         axes(handles.ImageTraining);
%         imshow(I); 
        
    % collect the feature value
    dataset(i,:)=[mean_red mean_green mean_blue diameter];
    
end

class';
% 
% dataset_class{1}={mat2cell(dataset) class{1}'}
% 
% dataset_class{1}

dataset;

set(handles.NameImageTraining,'String',sprintf(strcat('Image Training (Color) DONE !')));

merge_data={dataset,class'};

% lots of data multiply with lots of feature
data_multiply_feature=(lots_of_data_train_jn+lots_of_data_train_jl+lots_of_data_train_jm)*lots_of_feature;

% collect data feature R, G, B and Diameter
for i=1:data_multiply_feature
   dat_init{i}=num2str(merge_data{1}(i),'%.3f');
end

% collect data class Citrus : Nipis, Lemon and Orange
for j=1:(lots_of_data_train_jn+lots_of_data_train_jl+lots_of_data_train_jm)
   dat_init{data_multiply_feature+j}=char(merge_data{2}(j));
end

dat=reshape(dat_init,[(lots_of_data_train_jn+lots_of_data_train_jl+lots_of_data_train_jm) (lots_of_feature+1)]);

set(NaiveBayesClassifierProject.DataTraining,'Userdata',dat);

%% insert data into dataset
t=uitable('Data', dat, 'ColumnName',...
        {'R (Red)', 'G (Green)', 'B (Blue)', 'D (Diameter)','Class (Citrus)'},...
        'Position', [20 20 430 150]);
    
   
    
    set(NaiveBayesClassifierProject.dataset_all_feature_class,'Userdata',dataset);
    
% --- Executes during object creation, after setting all properties.
function SizeImageTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeImageTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function SizeImageTesting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeImageTesting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function NameImageTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NameImageTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uitabledataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitabledataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in uitabledataset.
function uitabledataset_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitabledataset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in uitabledataset.
function uitabledataset_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitabledataset (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in testingdata.
function testingdata_Callback(hObject, eventdata, handles)
% hObject    handle to testingdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NaiveBayesClassifierProject=guidata(gcbo);

% data training of Citrus nipis (jn)
lots_of_data_train_jn=7;

% data training of Citrus lemon (jl)
lots_of_data_train_jl=3;

% data training of Citrus orange (jm)
lots_of_data_train_jm=5;

lots_of_feature=4;
lots_of_class=3;

[basefilename,path]= uigetfile({'*.*'},'Open All Image File');
filename= fullfile(path, basefilename);

if sum(strfind(lower(basefilename), '.'))==0
else
    
    I_testing = imread (filename);

    % if I = [MxNx4]
    if(size(I_testing,3)==4)
        I_testing(:,:,1)=[]; % convert to I = [MxNx3]
    end
        
    %set(NaiveBayesClassifierProject.NaiveBayesClassifier,'CurrentAxes',NaiveBayesClassifierProject.ImageTraining);
    %set (imshow(I));
    
    axes(handles.ImageTesting);
    imshow(I_testing);
    
    % count mean Red, Green, Blue
    mean_red_testing=mean(mean(I_testing(:,:,1)));
    mean_green_testing=mean(mean(I_testing(:,:,2)));
    mean_blue_testing=mean(mean(I_testing(:,:,3)));
    
    set(NaiveBayesClassifierProject.var_mean_red_testing,...
        'String',num2str(mean_red_testing,'%.3f'));
    
    set(NaiveBayesClassifierProject.var_mean_green_testing,...
        'String',num2str(mean_green_testing,'%.3f'));
    
    set(NaiveBayesClassifierProject.var_mean_blue_testing,...
        'String',num2str(mean_blue_testing,'%.3f'));
    
    set(NaiveBayesClassifierProject.ImageTesting,'Userdata',I_testing);   
    
    SizeImageTesting=size(I_testing);
    StringSizeImageTesting=sprintf(strcat(num2str(SizeImageTesting(1)),'x',num2str(SizeImageTesting(2)),'x',num2str(SizeImageTesting(3))));    
    StringSizeImageTesting = strrep(StringSizeImageTesting,'x',' x ');
    
    set(handles.SizeImageTesting,'String',StringSizeImageTesting);
    
    %% create gray-scale image
    I_gray_testing=Function_ColorToGray(I_testing);
    
    axes(handles.ImageGrayTesting);
    imshow(I_gray_testing);  
    
    %%
    
    %% create binary image _testing
%     level=graythresh(I_testing);
%     I_biner_testing=im2bw(I_gray_testing,level);
    
    I_biner_testing=zeros(size(I_gray_testing,1),size(I_gray_testing,2));
    I_biner_testing(find(I_gray_testing<255))=1;

    axes(handles.ImageBinaryTesting);
    imshow(I_biner_testing);  
    
    %% create max filter image from biner_testing image
       %windowing_size must be valued an odd number >=3
       windowing_size=5;
       max_filter_I_biner_testing=Function_MaxFilterBiner_(I_biner_testing,windowing_size); 
       
       axes(handles.ImageMaxFilterTesting);
       imshow(max_filter_I_biner_testing);
       
       
    %% count diameter with unit length each pixel
       % determine index which contains the value 1 (white color)
         [idx_,idy_]=find(max_filter_I_biner_testing==1);
         diameter_testing=idy_(numel(idy_))-idy_(1)+1;
         
         set(NaiveBayesClassifierProject.var_diameter_testing,...
        'String',num2str(diameter_testing,'%.2f'));
         
    % collect feature value
    dataset(1,:)=[mean_red_testing mean_green_testing mean_blue_testing diameter_testing];
    
    if(isempty(strfind(basefilename, 'Lemon'))==0)
        class_testing{1}='Lemon';
    elseif(isempty(strfind(basefilename, 'Orange'))==0)
        class_testing{1}='Orange';
    elseif(isempty(strfind(basefilename, 'Nipis'))==0)
        class_testing{1}='Nipis'; 
    else
        class_testing{1}='UnKnown';
    end
    
    class_testing{1}
    
    set(NaiveBayesClassifierProject.var_class_testing,...
        'String',char(class_testing{1}));
    
    % replace pixel value
        I_red=I_testing(:,:,1);
        I_red(:,idy_(1):(idy_(1)+10))=105;
        I_red(:,(idy_(numel(idy_))-10):idy_(numel(idy_)))=105;
        I_testing(:,:,1)=I_red;
        
        I_green=I_testing(:,:,2);
        I_green(:,idy_(1):(idy_(1)+10))=75;
        I_green(:,(idy_(numel(idy_))-10):idy_(numel(idy_)))=75;
        I_testing(:,:,2)=I_green;
        
        I_blue=I_testing(:,:,3);
        I_blue(:,idy_(1):(idy_(1)+10))=245; 
        I_blue(:,(idy_(numel(idy_))-10):idy_(numel(idy_)))=245;
        I_testing(:,:,3)=I_blue;
        
        
        axes(handles.ImageTesting);
        imshow(I_testing); 
        
        
      %% count Probability of Posterior
        %GetDataTraining=get(NaiveBayesClassifierProject.DataTraining,'User
        %data');
        %Getmean_varian=get(NaiveBayesClassifierProject.dataset_all_feature_class,'Userdata');      
        dataset=get(NaiveBayesClassifierProject.dataset_all_feature_class,'Userdata');      
        
        Xrgbd=[mean_red_testing;mean_green_testing;mean_blue_testing;...
            diameter_testing];
        
%         % normalization dataset    
%         dataset
%         Xrgbd=Xrgbd'
%         dataset_Xrgbd=[dataset;Xrgbd];
%         dataset_Xrgbd_norm=dataset_Xrgbd./(padarray((sum(dataset_Xrgbd.*dataset_Xrgbd).^0.5),[(size(dataset_Xrgbd,1)-1) 0],'replicate','post'))
%         
%         % replace dataset and Xrgbd with normalization
%         dataset = dataset_Xrgbd_norm(1:(size(dataset_Xrgbd_norm,1)-1),:)
%         Xrgbd = dataset_Xrgbd_norm(size(dataset_Xrgbd_norm,1),:)
        
        %% count & save result of mean and varian from every feature & class
            % Initialization
            Getmean_varian=zeros(2*lots_of_class,lots_of_feature);

            for i=1:lots_of_feature
                % set feature 1 for R, 2 for G, 3 for B, 4 for D
                feature_rgbd=dataset(:,i);

                % count mean_class_citrus_nipis(jn),_lemon(jl),_orange(jm)
                mean_feature_rgbd_jn=mean(feature_rgbd(1:lots_of_data_train_jn));
                mean_feature_rgbd_jl=mean(feature_rgbd(...
                    (lots_of_data_train_jn+1):(lots_of_data_train_jn+lots_of_data_train_jl)));
                mean_feature_rgbd_jm=mean(feature_rgbd(...
                    (lots_of_data_train_jn+lots_of_data_train_jl+1):...
                    (lots_of_data_train_jn+lots_of_data_train_jl+lots_of_data_train_jm)));

                % count varian_class_citrus_nipis(jn),_lemon(jl),_orange(jm
                varian_feature_rgbd_jn=var(feature_rgbd(1:lots_of_data_train_jn));
                varian_feature_rgbd_jl=var(feature_rgbd(...
                    (lots_of_data_train_jn+1):(lots_of_data_train_jn+lots_of_data_train_jl)));
                varian_feature_rgbd_jm=var(feature_rgbd(...
                    (lots_of_data_train_jn+lots_of_data_train_jl+1):...
                    (lots_of_data_train_jn+lots_of_data_train_jl+lots_of_data_train_jm)));

                Getmean_varian(:,i)=[mean_feature_rgbd_jn,mean_feature_rgbd_jl,...
                    mean_feature_rgbd_jm,varian_feature_rgbd_jn,varian_feature_rgbd_jl,...
                    varian_feature_rgbd_jm];
            end

        %% note of format mean_varian matrix :
           % size : 
             % rows = 2*lots_of_class
             % column = lots_of_feature

        % content mean_varian matrix :
         % -----------------------------------------------------------------------------------
         % |       Red        |     Green         |      Blue         |      Diameter        | 
         % -----------------------------------------------------------------------------------
         % |   mean_jn_red    |   mean_jn_green   |  mean_jn_blue     |   mean_jn_diameter   |
         % |   mean_jl_red    |   mean_jl_green   |  mean_jl_blue     |   mean_jl_diameter   |
         % |   mean_jm_red    |   mean_jm_green   |  mean_jm_blue     |   mean_jm_diameter   |
         % |  varian_jn_red   |  varian_jn_green  |  varian_jn_blue   |  varian_jn_diameter  |
         % |  varian_jl_red   |  varian_jl_green  |  varian_jl_blue   |  varian_jl_diameter  |
         % |  varian_jm_red   |  varian_jm_green  |  varian_jm_blue   |  varian_jm_diameter  |
         % -----------------------------------------------------------------------------------
        %%
        
        
        
        % count Probability of Prior
        P_Prior_jn=lots_of_data_train_jn/(lots_of_data_train_jn+...
        lots_of_data_train_jl+lots_of_data_train_jm);
        P_Prior_jl=lots_of_data_train_jl/(lots_of_data_train_jn+...
            lots_of_data_train_jl+lots_of_data_train_jm);
        P_Prior_jm=lots_of_data_train_jm/(lots_of_data_train_jn+...
            lots_of_data_train_jl+lots_of_data_train_jm);
        
        % initialization Probability value of Posterior
        P_Posterior_jn=1*P_Prior_jn
        P_Posterior_jl=1*P_Prior_jl
        P_Posterior_jm=1*P_Prior_jm
        
%% note format Getmean_varian matrix :
   % size : 
     % rows = 2*lots_of_class
     % column = lots_of_feature

% content of Getmean_varian matrix :
 % -----------------------------------------------------------------------------------
 % |       Red        |     Green         |      Blue         |      Diameter        | 
 % -----------------------------------------------------------------------------------
 % |   mean_jn_red    |   mean_jn_green   |  mean_jn_blue     |   mean_jn_diameter   |
 % |   mean_jl_red    |   mean_jl_green   |  mean_jl_blue     |   mean_jl_diameter   |
 % |   mean_jm_red    |   mean_jm_green   |  mean_jm_blue     |   mean_jm_diameter   |
 % |  varian_jn_red   |  varian_jn_green  |  varian_jn_blue   |  varian_jn_diameter  |
 % |  varian_jl_red   |  varian_jl_green  |  varian_jl_blue   |  varian_jl_diameter  |
 % |  varian_jm_red   |  varian_jm_green  |  varian_jm_blue   |  varian_jm_diameter  |
 % -----------------------------------------------------------------------------------
%%
  
        % count Probability of Likelihood from RBGD feature
        for i=1:lots_of_feature      
            mean_varian_RGBD=Getmean_varian(:,i);
            P_Likelihood_x_RGB_to_jn=...
                (1/sqrt(2*(22/7)*mean_varian_RGBD(4)))...
                *exp(-1*(((Xrgbd(i)-mean_varian_RGBD(1))^2)/(2*mean_varian_RGBD(4))))
            P_Posterior_jn=P_Posterior_jn*P_Likelihood_x_RGB_to_jn
            
            mean_varian_RGBD(5)
            mean_varian_RGBD(2)
            Xrgbd(i)
            P_Likelihood_x_RGB_to_jl=...
                (1/sqrt(2*(22/7)*mean_varian_RGBD(5)))...
                *exp(-1*(((Xrgbd(i)-mean_varian_RGBD(2))^2)/(2*mean_varian_RGBD(5))))
            P_Posterior_jl=P_Posterior_jl*P_Likelihood_x_RGB_to_jl
            
            P_Likelihood_x_RGB_to_jm=...
                (1/sqrt(2*(22/7)*mean_varian_RGBD(6)))...
                *exp(-1*(((Xrgbd(i)-mean_varian_RGBD(3))^2)/(2*mean_varian_RGBD(6))))
            P_Posterior_jm=P_Posterior_jm*P_Likelihood_x_RGB_to_jm
        end
        
        set(NaiveBayesClassifierProject.posterior_class_jn,...
        'String',strcat(num2str(P_Posterior_jn,'%.20f'),{'  ('},num2str(P_Posterior_jn),')'));
        set(NaiveBayesClassifierProject.posterior_class_jl,...
        'String',strcat(num2str(P_Posterior_jl,'%.20f'),{'  ('},num2str(P_Posterior_jl),')'));
        set(NaiveBayesClassifierProject.posterior_class_jm,...
        'String',strcat(num2str(P_Posterior_jm,'%.20f'),{'  ('},num2str(P_Posterior_jm),')'));
        
        All_P_Posterior=[P_Posterior_jn;P_Posterior_jl;P_Posterior_jm];
      
        [vmax_Posterior,idxmax_Posterior]=max(All_P_Posterior);
        
        Decision_Of_Classification='';
        if(idxmax_Posterior==1)
            Decision_Of_Classification='Nipis';
        elseif(idxmax_Posterior==2)
            Decision_Of_Classification='Lemon';
        else
            Decision_Of_Classification='Orange';
        end
        
        set(NaiveBayesClassifierProject.classification_result,...
        'String',Decision_Of_Classification);
        

end



function var_mean_red_testing_Callback(hObject, eventdata, handles)
% hObject    handle to var_mean_red_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_mean_red_testing as text
%        str2double(get(hObject,'String')) returns contents of var_mean_red_testing as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_mean_red_testing = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_mean_red_testing = var_mean_red_testing;

% --- Executes during object creation, after setting all properties.
function var_mean_red_testing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_mean_red_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function var_mean_green_testing_Callback(hObject, eventdata, handles)
% hObject    handle to var_mean_green_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_mean_green_testing as text
%        str2double(get(hObject,'String')) returns contents of var_mean_green_testing as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_mean_green_testing = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_mean_green_testing = var_mean_green_testing;


% --- Executes during object creation, after setting all properties.
function var_mean_green_testing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_mean_green_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function var_mean_blue_testing_Callback(hObject, eventdata, handles)
% hObject    handle to var_mean_blue_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_mean_blue_testing as text
%        str2double(get(hObject,'String')) returns contents of var_mean_blue_testing as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_mean_blur_testing = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_mean_blue_testing = var_mean_blue_testing;


% --- Executes during object creation, after setting all properties.
function var_mean_blue_testing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_mean_blue_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function var_diameter_testing_Callback(hObject, eventdata, handles)
% hObject    handle to var_diameter_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_diameter_testing as text
%        str2double(get(hObject,'String')) returns contents of var_diameter_testing as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_diameter_testing = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_diameter_testing = var_diameter_testing;


% --- Executes during object creation, after setting all properties.
function var_diameter_testing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_diameter_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function var_class_testing_Callback(hObject, eventdata, handles)
% hObject    handle to var_class_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_class_testing as text
%        str2double(get(hObject,'String')) returns contents of var_class_testing as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_class_testing = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_class_testing = var_class_testing;



% --- Executes during object creation, after setting all properties.
function var_class_testing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_class_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posterior_class_jn_Callback(hObject, eventdata, handles)
% hObject    handle to posterior_class_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posterior_class_jn as text
%        str2double(get(hObject,'String')) returns contents of posterior_class_jn as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_class_jn = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_class_jn = var_class_jn;

% --- Executes during object creation, after setting all properties.
function posterior_class_jn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posterior_class_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posterior_class_jl_Callback(hObject, eventdata, handles)
% hObject    handle to posterior_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posterior_class_jl as text
%        str2double(get(hObject,'String')) returns contents of posterior_class_jl as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_class_jl = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_class_jl = var_class_jl;

% --- Executes during object creation, after setting all properties.
function posterior_class_jl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posterior_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posterior_class_jm_Callback(hObject, eventdata, handles)
% hObject    handle to posterior_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posterior_class_jm as text
%        str2double(get(hObject,'String')) returns contents of posterior_class_jm as a double
NaiveBayesClassifierProject = guidata(gcbo);
var_class_jm = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.var_class_jm = var_class_jm;


% --- Executes during object creation, after setting all properties.
function posterior_class_jm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posterior_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_class_jn_Callback(hObject, eventdata, handles)
% hObject    handle to txt_mean_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_mean_var as text
%        str2double(get(hObject,'String')) returns contents of txt_mean_var as a double
NaiveBayesClassifierProject = guidata(gcbo);
mean_class_jn = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.mean_class_jn = mean_class_jn;

% --- Executes during object creation, after setting all properties.
function mean_class_jn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_mean_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function varian_class_jn_Callback(hObject, eventdata, handles)
% hObject    handle to varian_class_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varian_class_jn as text
%        str2double(get(hObject,'String')) returns contents of varian_class_jn as a double
NaiveBayesClassifierProject = guidata(gcbo);
varian_class_jn = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.varian_class_jn = varian_class_jn;

% --- Executes during object creation, after setting all properties.
function varian_class_jn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varian_class_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_class_jl_Callback(hObject, eventdata, handles)
% hObject    handle to mean_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_class_jl as text
%        str2double(get(hObject,'String')) returns contents of mean_class_jl as a double
NaiveBayesClassifierProject = guidata(gcbo);
mean_class_jl = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.mean_class_jl = mean_class_jl;

% --- Executes during object creation, after setting all properties.
function mean_class_jl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function varian_class_jl_Callback(hObject, eventdata, handles)
% hObject    handle to varian_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varian_class_jl as text
%        str2double(get(hObject,'String')) returns contents of varian_class_jl as a double
NaiveBayesClassifierProject = guidata(gcbo);
varian_class_jl = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.varian_class_jl = varian_class_jl;

% --- Executes during object creation, after setting all properties.
function varian_class_jl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varian_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_class_jm_Callback(hObject, eventdata, handles)
% hObject    handle to mean_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_class_jm as text
%        str2double(get(hObject,'String')) returns contents of mean_class_jm as a double
NaiveBayesClassifierProject = guidata(gcbo);
mean_class_jm = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.mean_class_jm = mean_class_jm;

% --- Executes during object creation, after setting all properties.
function mean_class_jm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function varian_class_jm_Callback(hObject, eventdata, handles)
% hObject    handle to varian_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varian_class_jm as text
%        str2double(get(hObject,'String')) returns contents of varian_class_jm as a double
NaiveBayesClassifierProject = guidata(gcbo);
varian_class_jm = str2double(get(hObject, 'String'));
NaiveBayesClassifierProject.varian_class_jm = varian_class_jm;

% --- Executes during object creation, after setting all properties.
function varian_class_jm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varian_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classification_result_Callback(hObject, eventdata, handles)
% hObject    handle to classification_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classification_result as text
%        str2double(get(hObject,'String')) returns contents of classification_result as a double


% --- Executes during object creation, after setting all properties.
function classification_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classification_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ImageTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ImageTraining


% --- Executes during object creation, after setting all properties.
function ImageTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ImageTraining


% --- Executes during object creation, after setting all properties.
function SizeImageTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeImageTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function NameImageTraining_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NameImageTraining (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function var_class_testing_Callback(hObject, eventdata, handles)
% hObject    handle to var_class_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_class_testing as text
%        str2double(get(hObject,'String')) returns contents of var_class_testing as a double


% --- Executes during object creation, after setting all properties.
function var_class_testing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_class_testing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posterior_class_jn_Callback(hObject, eventdata, handles)
% hObject    handle to posterior_class_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posterior_class_jn as text
%        str2double(get(hObject,'String')) returns contents of posterior_class_jn as a double


% --- Executes during object creation, after setting all properties.
function posterior_class_jn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posterior_class_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posterior_class_jl_Callback(hObject, eventdata, handles)
% hObject    handle to posterior_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posterior_class_jl as text
%        str2double(get(hObject,'String')) returns contents of posterior_class_jl as a double


% --- Executes during object creation, after setting all properties.
function posterior_class_jl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posterior_class_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posterior_class_jm_Callback(hObject, eventdata, handles)
% hObject    handle to posterior_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posterior_class_jm as text
%        str2double(get(hObject,'String')) returns contents of posterior_class_jm as a double


% --- Executes during object creation, after setting all properties.
function posterior_class_jm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posterior_class_jm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
