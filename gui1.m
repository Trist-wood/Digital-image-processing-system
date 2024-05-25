function varargout = gui1(varargin)
% GUI1 MATLAB code for gui1.fig
%      GUI1, by itself, creates a new GUI1 or raises the existing
%      singleton*.
%
%      H = GUI1 returns the handle to a new GUI1 or the handle to
%      the existing singleton*.
%
%      GUI1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI1.M with the given input arguments.
%
%      GUI1('Property','Value',...) creates a new GUI1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui1

% Last Modified by GUIDE v2.5 07-Dec-2021 18:53:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui1_OpeningFcn, ...
                   'gui_OutputFcn',  @gui1_OutputFcn, ...
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


function gui1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui1 (see VARARGIN)

% Choose default command line output for gui1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton(选择).
function pushbutton5_Callback(hObject, eventdata, handles)
%axis off  %%关闭坐标轴显示  
[filename,pathname] =uigetfile({'*.jpg';'*.bmp';'*.*'},'打开图片');
global im
str=[pathname,filename];  
%%打开图像  
im=imread(str);  
%%打开axes1的句柄 进行axes1的操作  
axes(handles.axes1);  
%%在axes1中显示 图像  
imshow(im); 
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton2(灰度图).
function pushbutton2_Callback(hObject, eventdata, handles)
%I = imread('D:\matlabruanjian\bin\123.jpg');   %读取路径下的图片
%axes(handles.axes1);              %在第一个轴中显示
%imshow(I);title('原图');
global im
global out
out = rgb2gray(im);
axes(handles.axes2); 
imshow(out);title('灰度图');
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4(直方图).
function pushbutton4_Callback(hObject, eventdata, handles)
global im
global out

im1 = rgb2gray(im);
out = imhist(im1);
axes(handles.axes2); 
imshow(out);title('直方图');
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8（膨胀）.
function pushbutton8_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
out = imdilate(im,x);
axes(handles.axes2); 
imshow(out);title('膨胀处理');

% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10（腐蚀）.
function pushbutton10_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
out = imerode(im,x);
axes(handles.axes2); 
imshow(out);title('腐蚀处理');

% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11（开运算）.
function pushbutton11_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
im1 = rgb2gray(im);
out = imopen(im1,x);
axes(handles.axes2); 
imshow(out);title('开运算');
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12（闭运算）.
function pushbutton12_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
im1 = rgb2gray(im);
out = imclose(im1,x);
axes(handles.axes2); 
imshow(out);title('闭运算');
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13（线检测）.
function pushbutton13_Callback(hObject, eventdata, handles)
global im;
global out;

w = [2 -1 -1;-1 2 -1;-1 -1 2];          % -45°方向检测线
g = imfilter(double(im),w);
axes(handles.axes2); 
imshow(out);title('线检测');

% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton14(Sobel算子检测).
function pushbutton14_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Sobel',0.06);
axes(handles.axes2); 
imshow(out);title('Sobel算子边缘检测');

% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15(Canny算子检测).
function pushbutton15_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Canny',0.06);
axes(handles.axes2); 
imshow(out);title('Canny算子边缘检测');

% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16(Log算子检测).
function pushbutton16_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Log',0.012);
axes(handles.axes2); 
imshow(out);title('Log算子边缘检测');

% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17(Prewitt算子检测).
function pushbutton17_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Prewitt',0.06);
axes(handles.axes2); 
imshow(out);title('Prewitt算子边缘检测');

% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18(Roberts算子检测).
function pushbutton18_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Roberts',0.06);
axes(handles.axes2); 
imshow(out);title('Roberts算子边缘检测');
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function jibengongneng_Callback(hObject, eventdata, handles)
% hObject    handle to jibengongneng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function huiduhua_Callback(hObject, eventdata, handles)
global im
global out
out = rgb2gray(im);
axes(handles.axes2); 
imshow(out);title('灰度图');
% hObject    handle to huiduhua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function zhifangtu_Callback(hObject, eventdata, handles)
global im
global out

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  f = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
f = rgb;
end

out = f;
axes(handles.axes2); 
imhist(f);title('直方图');


% hObject    handle to zhifangtu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function erzhituxiangchuli_Callback(hObject, eventdata, handles)
% hObject    handle to erzhituxiangchuli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function pengzhang_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
out = imdilate(im,x);
axes(handles.axes2); 
imshow(out);title('膨胀处理');
% hObject    handle to pengzhang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fushi_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
out = imerode(im,x);
axes(handles.axes2); 
imshow(out);title('腐蚀处理');

% hObject    handle to fushi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function kaiyunsuan_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

out = imopen(im1,x);
axes(handles.axes2); 
imshow(out);title('开运算');
% hObject    handle to kaiyunsuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function biyuansuan_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end
out = imclose(im1,x);
axes(handles.axes2); 
imshow(out);title('闭运算');
% hObject    handle to biyuansuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fenge_Callback(hObject, eventdata, handles)
% hObject    handle to fenge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dian_Callback(hObject, eventdata, handles)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

w = [-1 -1 -1;-1 8 -1;-1 -1 -1];%给定模板
out = ordfilt2(im1,5*5,ones(5,5))-ordfilt2(im1,1,ones(5,5));%在这里采用了5*5模板进行差值处理
T2 = max(out(:));%同理
out = out>=T2/2;%同理
axes(handles.axes2); 
imshow(out);title('点检测');

% hObject    handle to dian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function xian_Callback(hObject, eventdata, handles)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

w = [2 -1 -1;-1 2 -1;-1 -1 2];          % -45°方向检测线
out = imfilter(double(im1),w);
axes(handles.axes2); 
imshow(out);title('线检测');
% hObject    handle to xian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
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
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
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

% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Log_Callback(hObject, eventdata, handles)
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
% hObject    handle to Log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
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
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Roberts_Callback(hObject, eventdata, handles)
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
% hObject    handle to Roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Otsu_Callback(hObject, eventdata, handles)
global im;
global out;

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im1 = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im1 = rgb;
end

I = im2double(im);
k = graythresh(im1);              %得到最优阈值
out = im2bw(im1,k);                 %转换成二值图，k为分割阈值
axes(handles.axes2); 
imshow(out);title('Otsu算子边缘检测');  
% hObject    handle to Otsu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tuxiangdezengqiang_Callback(hObject, eventdata, handles)
% hObject    handle to tuxiangdezengqiang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function jubupinghua_Callback(hObject, eventdata, handles)
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

% hObject    handle to jubupinghua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function zhongzhilvbo_Callback(hObject, eventdata, handles)
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

% hObject    handle to zhongzhilvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tiduruihua_Callback(hObject, eventdata, handles)
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

% hObject    handle to tiduruihua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaosiditonglvbo_Callback(hObject, eventdata, handles)
global im;
global out;

d0=50;  %阈值
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  image = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
image = rgb;
end
[M ,N]=size(image);

img_f = fft2(double(image));%傅里叶变换得到频谱
img_f = fftshift(img_f);  %移到中间

m_mid = floor(M/2);%中心点坐标
n_mid = floor(N/2);  

h = zeros(M,N);%高斯低通滤波器构造
for i = 1:M
    for j = 1:N
        d = ((i-m_mid)^2+(j-n_mid)^2);
        h(i,j) = exp(-(d)/(2*(d0^2)));      
    end
end

img_lpf = h.*img_f;
img_lpf=ifftshift(img_lpf);    %中心平移回原来状态
out = uint8(real(ifft2(img_lpf)));  %反傅里叶变换,取实数部分
axes(handles.axes2); 
imshow(out,[]);title('高斯低通滤波器后的图像'); 

% hObject    handle to gaosiditonglvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaosigaotonglvbo_Callback(hObject, eventdata, handles)
global im;
global out;

d0=50;  %阈值
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  image = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
image = rgb;
end
[M ,N]=size(image);

img_f = fft2(double(image));%傅里叶变换得到频谱
img_f = fftshift(img_f);  %移到中间

m_mid = floor(M/2);%中心点坐标
n_mid = floor(N/2);  

h = zeros(M,N);%高斯低通滤波器构造
for i = 1:M
    for j = 1:N
        d = ((i-m_mid)^2+(j-n_mid)^2);
        h(i,j) = 1-exp(-(d)/(2*(d0^2)));      
    end
end

img_lpf = h.*img_f;
img_lpf=ifftshift(img_lpf);    %中心平移回原来状态
out = uint8(real(ifft2(img_lpf)));  %反傅里叶变换,取实数部分
axes(handles.axes2); 
imshow(out,[]);title('高斯高通滤波器后的图像'); 
% hObject    handle to gaosigaotonglvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function batewosi_Callback(hObject, eventdata, handles)
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
imshow(out);title('巴特沃斯高通滤波器后的图像'); 
% hObject    handle to batewosi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function shuzituxiangbianhaun_Callback(hObject, eventdata, handles)
% hObject    handle to shuzituxiangbianhaun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DFT_Callback(hObject, eventdata, handles)
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

% hObject    handle to DFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DCT_Callback(hObject, eventdata, handles)
global im
global out  

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im = rgb;
end

I1 = double(im);
I2 = dct2(I1);
out = log(abs(I2));
axes(handles.axes2); 
imshow(out);title('离散余弦变换'); 

% hObject    handle to DCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FTT_Callback(hObject, eventdata, handles)
global im
global out  

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  im = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
im = rgb;
end

F = fft2(im2double(im));  %FFT
F = fftshift(F);  %FFT频谱平移
F = real(F);
out= log(F+1);  %频谱对数变换
axes(handles.axes2); 
imshow(out);title('频谱图'); 

% hObject    handle to FTT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DWT_Callback(hObject, eventdata, handles)
global im
global out  
 
%-------------------小波变换一级重构，小波基函数选db4-----------------------

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  gray_pic = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
gray_pic = rgb;
end

figure('name','小波变换一级重构');
[c,s] = wavedec2(gray_pic,1,'db4'); %小波一级分解,小波基函数采用db4
re_ca1 = wrcoef2('a',c,s,'db4',1); %重建第一层低频分量系数
re_ch1 = wrcoef2('h',c,s,'db4',1); %重建第一层高频水平分量系数
re_cv1 = wrcoef2('v',c,s,'db4',1); %重建第一层高频垂直分量系数
re_cd1 = wrcoef2('d',c,s,'db4',1); %重建第一层高频对角分量系数
re_set1 = [re_ca1,re_ch1;re_cv1,re_cd1];  %将各个分量图像拼接在一张图像
subplot(1,2,1);imshow(re_set1,[]);title('第一层小波系数的重构');
out = re_ca1+re_ch1+re_cv1+re_cd1;%将各个分量合并复原
axes(handles.axes2); 
imshow(out,[]);title('一级重构图像'); 



% hObject    handle to DWT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tuxiangxuanzhuan_Callback(hObject, eventdata, handles)

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


% hObject    handle to tuxiangxuanzhuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function zhifangtujunhenghua_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  f = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
f = rgb;
end

out = histeq(f,256);
axes(handles.axes1); 
imshow(out);title('均衡化后的图像'); 
axes(handles.axes2); 
imhist(out);title('均衡化后的直方图'); 

% hObject    handle to zhifangtujunhenghua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function noise_Callback(hObject, eventdata, handles)

% hObject    handle to noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaosi_noise_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'gaussian');
out = J;
axes(handles.axes2); 
imshow(out);title('添加高斯噪声后的图像'); 

% hObject    handle to gaosi_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function posong_noise_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'poisson');
out = J;
axes(handles.axes2); 
imshow(out);title('添加泊松噪声后的图像'); 

% hObject    handle to posong_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function salt_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'salt');
out = J;
axes(handles.axes2); 
imshow(out);title('添加盐噪声后的图像'); 

% hObject    handle to salt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function speckle_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

I = im2double(im);
J = imnoise(I,'speckle');
out = J;
axes(handles.axes2); 
imshow(out);title('添加胡椒噪声后的图像'); 

% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function quyushengzhang_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

I = im;
if isinteger(I)
    I = im2double(I);
end 
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  I = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
I = rgb;
end
[M,N] = size(I);
[y,x] = getpts; %单击取点后，按enter结束
x1 = round(x);
y1 = round(y);
seed = I(x1,y1); %获取中心像素灰度值

J = zeros(M,N);
J(x1,y1) = 1;

count = 1; %待处理点个数
threshold = 0.15;
while count>0
    count = 0;
    for i = 1:M %遍历整幅图像
    for j = 1:N
        if J(i,j) == 1 %点在“栈”内
        if (i-1)>1&(i+1)<M&(j-1)>1&(j+1)<N %3*3邻域在图像范围内
            for u = -1:1 %8-邻域生长
            for v = -1:1
                if J(i+u,j+v) == 0&abs(I(i+u,j+v)-seed) <= threshold
                    J(i+u,j+v) = 1;
                    count = count+1;  %记录此次新生长的点个数
                end
            end
            end
        end
        end
    end
    end
end

out = J;
axes(handles.axes2);
imshow(out);title("segmented image")

% hObject    handle to quyushengzhang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function fenshuiling_Callback(hObject, eventdata, handles)
global im
global out  %定义全局变量 

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  g = rgb2gray(rgb); %将彩色图像转换为灰度图像
else
g = rgb;
end
gc=~g;   %对图像求补
D=bwdist(gc);   %计算其距离变换
L=watershed(-D);   %计算负距离变换的分水岭变换
w=L==0;    %L 的零值即分水岭的脊线像素
out =g & ~w;   %原始二值图像和图像 w 的 “补” 的逻辑 “与” 操作可完成分割
axes(handles.axes2);
imshow(out);title("分水岭分割图像")
% hObject    handle to fenshuiling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function baocun_Callback(hObject, eventdata, handles)
global out;

set(handles.axes2,'HandleVisibility','ON');
axes(handles.axes2);
[filename,pathname]=uiputfile({'*.jpg';'*.bmp';'*.tif';'*.*'},'save image as');
file=strcat(pathname,filename);
out = getimage(gca);
imwrite(out,file);
set(handles.axes2,'HandleVisibility','Off');

% hObject    handle to baocun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tuichu_Callback(hObject, eventdata, handles)
close(gcf)  %关闭当前Figure窗口句柄

% hObject    handle to tuichu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
