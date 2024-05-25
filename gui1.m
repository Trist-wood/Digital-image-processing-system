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

% --- Executes on button press in pushbutton(ѡ��).
function pushbutton5_Callback(hObject, eventdata, handles)
%axis off  %%�ر���������ʾ  
[filename,pathname] =uigetfile({'*.jpg';'*.bmp';'*.*'},'��ͼƬ');
global im
str=[pathname,filename];  
%%��ͼ��  
im=imread(str);  
%%��axes1�ľ�� ����axes1�Ĳ���  
axes(handles.axes1);  
%%��axes1����ʾ ͼ��  
imshow(im); 
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton2(�Ҷ�ͼ).
function pushbutton2_Callback(hObject, eventdata, handles)
%I = imread('D:\matlabruanjian\bin\123.jpg');   %��ȡ·���µ�ͼƬ
%axes(handles.axes1);              %�ڵ�һ��������ʾ
%imshow(I);title('ԭͼ');
global im
global out
out = rgb2gray(im);
axes(handles.axes2); 
imshow(out);title('�Ҷ�ͼ');
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4(ֱ��ͼ).
function pushbutton4_Callback(hObject, eventdata, handles)
global im
global out

im1 = rgb2gray(im);
out = imhist(im1);
axes(handles.axes2); 
imshow(out);title('ֱ��ͼ');
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8�����ͣ�.
function pushbutton8_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
out = imdilate(im,x);
axes(handles.axes2); 
imshow(out);title('���ʹ���');

% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10����ʴ��.
function pushbutton10_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
out = imerode(im,x);
axes(handles.axes2); 
imshow(out);title('��ʴ����');

% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11�������㣩.
function pushbutton11_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
im1 = rgb2gray(im);
out = imopen(im1,x);
axes(handles.axes2); 
imshow(out);title('������');
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12�������㣩.
function pushbutton12_Callback(hObject, eventdata, handles)
global im
global out

x = [1,1,1;1,1,1;1,1,1];
im1 = rgb2gray(im);
out = imclose(im1,x);
axes(handles.axes2); 
imshow(out);title('������');
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13���߼�⣩.
function pushbutton13_Callback(hObject, eventdata, handles)
global im;
global out;

w = [2 -1 -1;-1 2 -1;-1 -1 2];          % -45�㷽������
g = imfilter(double(im),w);
axes(handles.axes2); 
imshow(out);title('�߼��');

% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton14(Sobel���Ӽ��).
function pushbutton14_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Sobel',0.06);
axes(handles.axes2); 
imshow(out);title('Sobel���ӱ�Ե���');

% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15(Canny���Ӽ��).
function pushbutton15_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Canny',0.06);
axes(handles.axes2); 
imshow(out);title('Canny���ӱ�Ե���');

% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16(Log���Ӽ��).
function pushbutton16_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Log',0.012);
axes(handles.axes2); 
imshow(out);title('Log���ӱ�Ե���');

% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17(Prewitt���Ӽ��).
function pushbutton17_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Prewitt',0.06);
axes(handles.axes2); 
imshow(out);title('Prewitt���ӱ�Ե���');

% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18(Roberts���Ӽ��).
function pushbutton18_Callback(hObject, eventdata, handles)
global im;
global out;

im1 = rgb2gray(im);
out = edge(im1,'Roberts',0.06);
axes(handles.axes2); 
imshow(out);title('Roberts���ӱ�Ե���');
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
imshow(out);title('�Ҷ�ͼ');
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
  f = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
f = rgb;
end

out = f;
axes(handles.axes2); 
imhist(f);title('ֱ��ͼ');


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
imshow(out);title('���ʹ���');
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
imshow(out);title('��ʴ����');

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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

out = imopen(im1,x);
axes(handles.axes2); 
imshow(out);title('������');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end
out = imclose(im1,x);
axes(handles.axes2); 
imshow(out);title('������');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

w = [-1 -1 -1;-1 8 -1;-1 -1 -1];%����ģ��
out = ordfilt2(im1,5*5,ones(5,5))-ordfilt2(im1,1,ones(5,5));%�����������5*5ģ����в�ֵ����
T2 = max(out(:));%ͬ��
out = out>=T2/2;%ͬ��
axes(handles.axes2); 
imshow(out);title('����');

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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

w = [2 -1 -1;-1 2 -1;-1 -1 2];          % -45�㷽������
out = imfilter(double(im1),w);
axes(handles.axes2); 
imshow(out);title('�߼��');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

out = edge(im1,'Sobel',0.06);
axes(handles.axes2); 
imshow(out);title('Sobel���ӱ�Ե���');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

out = edge(im1,'Canny',0.06);
axes(handles.axes2); 
imshow(out);title('Canny���ӱ�Ե���');

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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

out = edge(im1,'Log',0.012);
axes(handles.axes2); 
imshow(out);title('Log���ӱ�Ե���');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

out = edge(im1,'Prewitt',0.06);
axes(handles.axes2); 
imshow(out);title('Prewitt���ӱ�Ե���');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

out = edge(im1,'Roberts',0.06);
axes(handles.axes2); 
imshow(out);title('Roberts���ӱ�Ե���');
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
  im1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im1 = rgb;
end

I = im2double(im);
k = graythresh(im1);              %�õ�������ֵ
out = im2bw(im1,k);                 %ת���ɶ�ֵͼ��kΪ�ָ���ֵ
axes(handles.axes2); 
imshow(out);title('Otsu���ӱ�Ե���');  
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

Inoised = imnoise(im,'gaussian',0.1,0.005);%��ͼ����и�˹��������
%�ƶ������
h=ones(3,3)/5;
h(1,1) = 0;
h(1,3) = 0;
h(3,1) = 0;
h(1,3) = 0;
%ƽ������
out = imfilter(Inoised,h);
axes(handles.axes2); 
imshow(out);title('ƽ���˲�������ͼ��'); 

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
  im = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im = rgb;
end
fn=imnoise(im,'salt & pepper',0.2);%�ú���imnosie��������������0.2����ͼ�а׵�ڵ���ֵĸ���Ϊ0.2
out = medfilt2(fn);%��ֵ�˲�
axes(handles.axes2); 
imshow(out);title('��ֵ�˲�������ͼ��'); 

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
  I1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
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
imshow(uint8(out));title('�񻯺��ͼ��'); 

% hObject    handle to tiduruihua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaosiditonglvbo_Callback(hObject, eventdata, handles)
global im;
global out;

d0=50;  %��ֵ
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  image = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
image = rgb;
end
[M ,N]=size(image);

img_f = fft2(double(image));%����Ҷ�任�õ�Ƶ��
img_f = fftshift(img_f);  %�Ƶ��м�

m_mid = floor(M/2);%���ĵ�����
n_mid = floor(N/2);  

h = zeros(M,N);%��˹��ͨ�˲�������
for i = 1:M
    for j = 1:N
        d = ((i-m_mid)^2+(j-n_mid)^2);
        h(i,j) = exp(-(d)/(2*(d0^2)));      
    end
end

img_lpf = h.*img_f;
img_lpf=ifftshift(img_lpf);    %����ƽ�ƻ�ԭ��״̬
out = uint8(real(ifft2(img_lpf)));  %������Ҷ�任,ȡʵ������
axes(handles.axes2); 
imshow(out,[]);title('��˹��ͨ�˲������ͼ��'); 

% hObject    handle to gaosiditonglvbo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaosigaotonglvbo_Callback(hObject, eventdata, handles)
global im;
global out;

d0=50;  %��ֵ
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  image = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
image = rgb;
end
[M ,N]=size(image);

img_f = fft2(double(image));%����Ҷ�任�õ�Ƶ��
img_f = fftshift(img_f);  %�Ƶ��м�

m_mid = floor(M/2);%���ĵ�����
n_mid = floor(N/2);  

h = zeros(M,N);%��˹��ͨ�˲�������
for i = 1:M
    for j = 1:N
        d = ((i-m_mid)^2+(j-n_mid)^2);
        h(i,j) = 1-exp(-(d)/(2*(d0^2)));      
    end
end

img_lpf = h.*img_f;
img_lpf=ifftshift(img_lpf);    %����ƽ�ƻ�ԭ��״̬
out = uint8(real(ifft2(img_lpf)));  %������Ҷ�任,ȡʵ������
axes(handles.axes2); 
imshow(out,[]);title('��˹��ͨ�˲������ͼ��'); 
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
  I1 = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
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
imshow(out);title('������˹��ͨ�˲������ͼ��'); 
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
  im = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im = rgb;
end

I1 = im2double(im);
I2 = fft2(I1);
I3 = fftshift(I2);
out = log(abs(I3)+1); 
axes(handles.axes2); 
imshow(out,[]);title('��ɢ����Ҷ�任'); 

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
  im = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im = rgb;
end

I1 = double(im);
I2 = dct2(I1);
out = log(abs(I2));
axes(handles.axes2); 
imshow(out);title('��ɢ���ұ任'); 

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
  im = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
im = rgb;
end

F = fft2(im2double(im));  %FFT
F = fftshift(F);  %FFTƵ��ƽ��
F = real(F);
out= log(F+1);  %Ƶ�׶����任
axes(handles.axes2); 
imshow(out);title('Ƶ��ͼ'); 

% hObject    handle to FTT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DWT_Callback(hObject, eventdata, handles)
global im
global out  
 
%-------------------С���任һ���ع���С��������ѡdb4-----------------------

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  gray_pic = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
gray_pic = rgb;
end

figure('name','С���任һ���ع�');
[c,s] = wavedec2(gray_pic,1,'db4'); %С��һ���ֽ�,С������������db4
re_ca1 = wrcoef2('a',c,s,'db4',1); %�ؽ���һ���Ƶ����ϵ��
re_ch1 = wrcoef2('h',c,s,'db4',1); %�ؽ���һ���Ƶˮƽ����ϵ��
re_cv1 = wrcoef2('v',c,s,'db4',1); %�ؽ���һ���Ƶ��ֱ����ϵ��
re_cd1 = wrcoef2('d',c,s,'db4',1); %�ؽ���һ���Ƶ�ԽǷ���ϵ��
re_set1 = [re_ca1,re_ch1;re_cv1,re_cd1];  %����������ͼ��ƴ����һ��ͼ��
subplot(1,2,1);imshow(re_set1,[]);title('��һ��С��ϵ�����ع�');
out = re_ca1+re_ch1+re_cv1+re_cd1;%�����������ϲ���ԭ
axes(handles.axes2); 
imshow(out,[]);title('һ���ع�ͼ��'); 



% hObject    handle to DWT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function tuxiangxuanzhuan_Callback(hObject, eventdata, handles)

global im
global out  %����ȫ�ֱ��� 

% A=getimage(handles.axes1);
A  = im;
axes(handles.axes2); 
prompt = {'������'};
def={'90'};
answer = inputdlg(prompt,'�����룺',1,def);
if ~isempty(answer)
a = str2num(answer{1});
J = imrotate(A,360-a);
out = J;
axes(handles.axes2); 
imshow(out,[]);title('��ת��ͼ��'); 
end


% hObject    handle to tuxiangxuanzhuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function zhifangtujunhenghua_Callback(hObject, eventdata, handles)
global im
global out  %����ȫ�ֱ��� 

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  f = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
f = rgb;
end

out = histeq(f,256);
axes(handles.axes1); 
imshow(out);title('���⻯���ͼ��'); 
axes(handles.axes2); 
imhist(out);title('���⻯���ֱ��ͼ'); 

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
global out  %����ȫ�ֱ��� 

I = im2double(im);
J = imnoise(I,'gaussian');
out = J;
axes(handles.axes2); 
imshow(out);title('��Ӹ�˹�������ͼ��'); 

% hObject    handle to gaosi_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function posong_noise_Callback(hObject, eventdata, handles)
global im
global out  %����ȫ�ֱ��� 

I = im2double(im);
J = imnoise(I,'poisson');
out = J;
axes(handles.axes2); 
imshow(out);title('��Ӳ����������ͼ��'); 

% hObject    handle to posong_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function salt_Callback(hObject, eventdata, handles)
global im
global out  %����ȫ�ֱ��� 

I = im2double(im);
J = imnoise(I,'salt');
out = J;
axes(handles.axes2); 
imshow(out);title('������������ͼ��'); 

% hObject    handle to salt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function speckle_Callback(hObject, eventdata, handles)
global im
global out  %����ȫ�ֱ��� 

I = im2double(im);
J = imnoise(I,'speckle');
out = J;
axes(handles.axes2); 
imshow(out);title('��Ӻ����������ͼ��'); 

% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function quyushengzhang_Callback(hObject, eventdata, handles)
global im
global out  %����ȫ�ֱ��� 

I = im;
if isinteger(I)
    I = im2double(I);
end 
rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  I = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
I = rgb;
end
[M,N] = size(I);
[y,x] = getpts; %����ȡ��󣬰�enter����
x1 = round(x);
y1 = round(y);
seed = I(x1,y1); %��ȡ�������ػҶ�ֵ

J = zeros(M,N);
J(x1,y1) = 1;

count = 1; %����������
threshold = 0.15;
while count>0
    count = 0;
    for i = 1:M %��������ͼ��
    for j = 1:N
        if J(i,j) == 1 %���ڡ�ջ����
        if (i-1)>1&(i+1)<M&(j-1)>1&(j+1)<N %3*3������ͼ��Χ��
            for u = -1:1 %8-��������
            for v = -1:1
                if J(i+u,j+v) == 0&abs(I(i+u,j+v)-seed) <= threshold
                    J(i+u,j+v) = 1;
                    count = count+1;  %��¼�˴��������ĵ����
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
global out  %����ȫ�ֱ��� 

rgb = im;
mysize = size(rgb);
if numel(mysize)>2
  g = rgb2gray(rgb); %����ɫͼ��ת��Ϊ�Ҷ�ͼ��
else
g = rgb;
end
gc=~g;   %��ͼ����
D=bwdist(gc);   %���������任
L=watershed(-D);   %���㸺����任�ķ�ˮ��任
w=L==0;    %L ����ֵ����ˮ��ļ�������
out =g & ~w;   %ԭʼ��ֵͼ���ͼ�� w �� ������ ���߼� ���롱 ��������ɷָ�
axes(handles.axes2);
imshow(out);title("��ˮ��ָ�ͼ��")
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
close(gcf)  %�رյ�ǰFigure���ھ��

% hObject    handle to tuichu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
