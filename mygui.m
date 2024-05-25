function varargout = mygui(varargin)
% MYGUI MATLAB code for mygui.fig
%      MYGUI, by itself, creates a new MYGUI or raises the existing
%      singleton*.
%
%      H = MYGUI returns the handle to a new MYGUI or the handle to
%      the existing singleton*.
%
%      MYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYGUI.M with the given input arguments.
%
%      MYGUI('Property','Value',...) creates a new MYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mygui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mygui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mygui

% Last Modified by GUIDE v2.5 02-May-2024 21:20:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mygui_OpeningFcn, ...
                   'gui_OutputFcn',  @mygui_OutputFcn, ...
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


% --- Executes just before mygui is made visible.
function mygui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mygui (see VARARGIN)

% Choose default command line output for mygui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mygui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mygui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im   %定义一个全局变量im
global im2
[filename,pathname]=...
    uigetfile({'*.*';'*.bmp';'*.tif';'*.png'},'select picture');  %选择图片路径
str=[pathname filename];  %合成路径+文件名
im=imread(str);   %读取图片
im2=im;
axes(handles.axes1);  %使用第一个axes
imshow(im);  %显示图片


% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global BW 
set(handles.axes2,'HandleVisibility','ON');
axes(handles.axes2);
[filename,pathname]=uiputfile({'*.jpg';'*.bmp';'*.tif';'*.*'},'save image as');
file=strcat(pathname,filename);
BW=getimage(gca);
imwrite(BW,file);
set(handles.axes2,'HandleVisibility','Off');

% --------------------------------------------------------------------
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)  %关闭当前Figure窗口句柄


% --------------------------------------------------------------------
function transformation_Callback(hObject, eventdata, handles)
% hObject    handle to transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function enhancement_Callback(hObject, eventdata, handles)
% hObject    handle to enhancement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Denoising_Callback(hObject, eventdata, handles)
% hObject    handle to Denoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Canny',0.06);
axes(handles.axes2); 
imshow(out);title('Canny算子边缘检测');

% --------------------------------------------------------------------
function xuanzhuan_Callback(hObject, eventdata, handles)
% hObject    handle to xuanzhuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

% A=getimage(handles.axes1);
A  = im;
axes(handles.axes2); 
prompt = {'度数：'};
def={'90'};
answer = inputdlg(prompt,'请输入：',1,def);
if ~isempty(answer)
a = str2num(answer{1});
J = imrotate(A,360-a);
out = J;
axes(handles.axes2); 
imshow(out,[]);title('旋转后图像'); 
end

% --------------------------------------------------------------------
function DFT_Callback(hObject, eventdata, handles)
% hObject    handle to DFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im = rgb;
end

I1 = im2double(im);
I2 = fft2(I1);
I3 = fftshift(I2);
out = log(abs(I3)+1); 
axes(handles.axes2); 
imshow(out,[]);title('离散傅里叶变换'); 


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global BW  %定义全局变量  
axes(handles.axes2); 
BW=rgb2gray(im);
im=BW;
imshow(BW);

% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Sobel',0.06);
axes(handles.axes2); 
imshow(out);title('Sobel算子边缘检测');

% --------------------------------------------------------------------
function Log_Callback(hObject, eventdata, handles)
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Log',0.012);
axes(handles.axes2); 
imshow(out);title('Log算子边缘检测');

% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Prewitt',0.06);
axes(handles.axes2); 
imshow(out);title('Prewitt算子边缘检测');

% --------------------------------------------------------------------
function Roberts_Callback(hObject, eventdata, handles)
% hObject    handle to Roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = edge(im1,'Roberts',0.06);
axes(handles.axes2); 
imshow(out);title('Roberts算子边缘检测');


% --------------------------------------------------------------------
function gaosi_noise_Callback(hObject, eventdata, handles)
% hObject    handle to gaosi_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'gaussian');
out = J;
axes(handles.axes2); 
imshow(out);title('添加高斯噪声后的图像'); 

% --------------------------------------------------------------------
function posong_noise_Callback(hObject, eventdata, handles)
% hObject    handle to posong_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'poisson');
out = J;
axes(handles.axes2); 
imshow(out);title('添加泊松噪声后的图像'); 

% --------------------------------------------------------------------
function salt_Callback(hObject, eventdata, handles)
% hObject    handle to salt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'salt');
out = J;
axes(handles.axes2); 
imshow(out);title('添加盐噪声后的图像'); 

% --------------------------------------------------------------------
function speckle_Callback(hObject, eventdata, handles)
% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'speckle');
out = J;
axes(handles.axes2); 
imshow(out);title('添加胡椒噪声后的图像'); 


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function zhongzhilvbo_Callback(hObject, eventdata, handles)
% hObject    handle to zhongzhilvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im = rgb;
end
fn=imnoise(im,'salt & pepper',0.2);%用函数imnosie产生椒盐噪声，0.2代表图中白点黑点出现的概率为0.2
out = medfilt2(fn);%中值滤波
axes(handles.axes2); 
imshow(out);title('中值滤波处理后的图像'); 


% --------------------------------------------------------------------
function batewosi_Callback(hObject, eventdata, handles)
% hObject    handle to batewosi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  I1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
I1 = rgb;
end
m = double(I1);
f = fft2(m);
f = fftshift(f);
[N1,N2] = size(f);
n1 = round(N1/2);
n2 = round(N2/2);
n = 2;
d0 = 10;
for i = 1:N1
    for j = 1:N2
        d = sqrt((i-n1)^2+(j-n2)^2);
        h = (1/(1+(d0/d)^(2*n)))+0.5;
        y(i,j) = h*f(i,j);
    end
end
y = ifftshift(y);
A = ifft2(y);
out = uint8(real(A));
axes(handles.axes2); 
imshow(out);title('巴特沃斯高通滤波后的图像'); 


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(gca,'XColor',get(gca,'Color')) ;% 这两行代码功能：将坐标轴和坐标刻度转为白色
set(gca,'YColor',get(gca,'Color'));
 
set(gca,'XTickLabel',[]); % 这两行代码功能：去除坐标刻度
set(gca,'YTickLabel',[]);
% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(gca,'XColor',get(gca,'Color')) ;% 这两行代码功能：将坐标轴和坐标刻度转为白色
set(gca,'YColor',get(gca,'Color'));
 
set(gca,'XTickLabel',[]); % 这两行代码功能：去除坐标刻度
set(gca,'YTickLabel',[]);
% Hint: place code in OpeningFcn to populate axes2


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

Inoised = imnoise(im,'gaussian',0.1,0.005);%对图像进行高斯噪声加噪
%制定卷积核
h=ones(3,3)/5;
h(1,1) = 0;
h(1,3) = 0;
h(3,1) = 0;
h(1,3) = 0;
%平滑运算
out = imfilter(Inoised,h);
axes(handles.axes2); 
imshow(out);title('平滑滤波处理后的图像'); 


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
global im2
global BW  %定义全局变量 
BW=im2;
im=im2;
axes(handles.axes2); 
imshow(BW);


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  I1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
I1 = rgb;
end
model=[-1,0,1;
       -2,0,2;
       -1,0,1];
[m,n]=size(I1);
I2=double(I1);

for i=2:m-1
    for j=2:n-1
        I2(i,j)=I1(i+1,j+1)+2*I1(i+1,j)+I1(i+1,j-1)-I1(i-1,j+1)-2*I1(i-1,j)-I1(i-1,j-1);
    end
end
I2 = I2 + double(I1);
out = I2;
axes(handles.axes2); 
imshow(uint8(out));title('锐化后的图像'); 
