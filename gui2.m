function varargout = gui(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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

function gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = gui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

% --------------------------------------------------------------------
%文件
function open_Callback(hObject, eventdata, handles)%打开图片
global im   %定义一个全局变量im
global im2
[filename,pathname]=...
    uigetfile({'*.*';'*.bmp';'*.tif';'*.png'},'select picture');  %选择图片路径
str=[pathname filename];  %合成路径+文件名
im=imread(str);   %读取图片
im2=im;
axes(handles.axes1);  %使用第一个axes
imshow(im);  %显示图片

function save_Callback(hObject, eventdata, handles)%保存图片
global BW 
set(handles.axes2,'HandleVisibility','ON');
axes(handles.axes2);
[filename,pathname]=uiputfile({'*.jpg';'*.bmp';'*.tif';'*.*'},'save image as');
file=strcat(pathname,filename);
BW=getimage(gca);
imwrite(BW,file);
set(handles.axes2,'HandleVisibility','Off');
function quit_Callback(hObject, ~, handles)%退出操作
close(gcf)  %关闭当前Figure窗口句柄

% --------------------------------------------------------------------
%菜单栏的调回函数，实际不使用
function t1_Callback(hObject, eventdata, handles)

function t2_Callback(hObject, eventdata, handles)

function t3_Callback(hObject, eventdata, handles)

function t4_Callback(hObject, eventdata, handles)

function t5_Callback(hObject, eventdata, handles)

function t6_Callback(hObject, eventdata, handles)

function t7_Callback(hObject, eventdata, handles)

function t8_Callback(hObject, eventdata, handles)

function t9_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
% 图像类型变换
function rgb2gray_Callback(hObject, eventdata, handles)%原图-灰度
global im
global BW  %定义全局变量  
axes(handles.axes2); 
BW=rgb2gray(im);
im=BW;
imshow(BW);

function im2bw_Callback(hObject, eventdata, handles)%原图-二值
global im
global BW  %定义全局变量  
axes(handles.axes2); 
BW=im2bw(im);
im=BW;
imshow(BW);

function dither_Callback(hObject, eventdata, handles)%灰度-二值
global im
global BW  %定义全局变量  
axes(handles.axes2); 
BW=dither(im);
im=BW;
imshow(BW);

% --------------------------------------------------------------------
% 边缘检测
function roberts_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
BW=edge(im,'roberts',0.04);
imshow(BW);

function sobel_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
BW=edge(im,'sobel',0.04);
imshow(BW);

function prewitt_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
BW=edge(im,'prewitt',0.04);
imshow(BW);

function log_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
BW=edge(im,'log',0.003);
imshow(BW);

function canny_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
BW=edge(im,'canny',0.2);
imshow(BW);


% --------------------------------------------------------------------
%图像变换
function DFT_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2); 
I1=double(im);
I2=fft2(I1);
I3=fftshift(I2);
I3=log(abs(I3));
BW=I3;
imshow(BW,[]);

function DCT_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2); 
I1=double(im);
I2=dct2(I1);
I3=log(abs(I2));
BW=I3;
imshow(BW);

% --------------------------------------------------------------------
% 图像旋转
function rotate_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量 
% A=getimage(handles.axes1);
A=im;
axes(handles.axes2); 
prompt={'度数：'};
def={'90'};
answer=inputdlg(prompt,'请输入：',1,def);
if ~isempty(answer)
a = str2num(answer{1});
J=imrotate(A,360-a);
BW=J;
imshow(BW);
end

function Initial_Callback(hObject, eventdata, handles)%初始化
global im
global im2
global BW  %定义全局变量 
BW=im2;
im=im2;
axes(handles.axes2); 
imshow(BW);

% --------------------------------------------------------------------
% 图像噪声添加
function gaussian_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量 
I=im2double(im);
J=imnoise(I,'gaussian');
BW=J;
axes(handles.axes2); 
imshow(BW);

function salt_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量 
I=im2double(im);
J=imnoise(I,'salt & pepper');
BW=J;
axes(handles.axes2); 
imshow(BW);

function speckle_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量 
I=im2double(im);
J=imnoise(I,'speckle');
BW=J;
axes(handles.axes2); 
imshow(BW);

function poisson_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量 
I=im2double(im);
J=imnoise(I,'poisson');
BW=J;
axes(handles.axes2); 
imshow(BW);

% --------------------------------------------------------------------
% 图像滤波
function medilt_Callback(hObject, eventdata, handles)%中值滤波
global BW  %定义全局变量 
J=medfilt2(BW, [3,3]);
BW=J;
axes(handles.axes2); 
imshow(BW);

function wiener_Callback(hObject, eventdata, handles)%自适应滤波
global BW  %定义全局变量 
J=wiener2(BW,[3,3]);
BW=J;
axes(handles.axes2); 
imshow(BW);

function filter2_Callback(hObject, eventdata, handles)%均值滤波
global BW  %定义全局变量 
M1=ones(3);
M1=M1/9;
J=filter2(M1,BW);
BW=J;
axes(handles.axes2); 
imshow(BW);

% --------------------------------------------------------------------
% 形态学图像处理
function bwmorph_Callback(hObject, eventdata, handles)%骨骼化
global im
global BW  %定义全局变量 
I=im2double(im);
I=im2bw(I);
J=bwmorph(I,'remove');
G=bwmorph(J,'skel',inf);
BW=G;
axes(handles.axes2); 
imshow(BW);

function imfill_Callback(hObject, eventdata, handles)%区域填充
global im
global BW  %定义全局变量 
axes(handles.axes2); 
I1=im2bw(im);

I2=1-I1;
se=ones(5);
I3=imerode(I2,se);
I4=1-I3;
I5=imerode(I4,se);
I6=imerode(I5,se);
I7=imdilate(I6,se);
BW=I7;

imshow(BW);

function diagonal_Callback(hObject, eventdata, handles)%对角线特征提取
global im
global BW  %定义全局变量 
axes(handles.axes2);
I1=im2bw(im);
v=[1,1,1,1,1,1,1,1,1,1];
se=diag(v);
I2=imerode(I1,se);
I3=imdilate(I2,se);
BW=I3;
imshow(BW);

% --------------------------------------------------------------------
%图像灰度变化
function plotchange_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
A=im2double(im);
    a=0.3;%0.3 0.7 0.5 0.9
    b=0.7;
    c=0.1;
    d=0.9;
    %0.3 0.7 0.1 0.9
    B=A;
    [m,n]=size(B);
    Mg=max(max(B));
    Mf=max(max(A));
    for (i=1:m)
      for (j=1:n)
        if(A(i,j)>=0&&A(i,j)<=a)
             B(i,j)=(c/a)*A(i,j);
        end
        if(A(i,j)>=a&&A(i,j)<=b)
            B(i,j)=(((d-c)/(b-a))*(A(i,j)-a))+c;
        end
        if(A(i,j)>=b&&A(i,j)<=1)
             B(i,j)=(((Mg-d)/(Mf-b))*(A(i,j)-b))+d;
        end
      end
    end
   BW=B;
   imshow(BW);

function imhist_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
 BW=im;
 imhist(BW);
 
function histeq_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
BW=histeq(im);
imhist(BW);


% --------------------------------------------------------------------
function histeqafter_Callback(hObject, eventdata, handles)
global im
global BW  %定义全局变量
axes(handles.axes2);   %使用第二个axes
imshow(BW);
