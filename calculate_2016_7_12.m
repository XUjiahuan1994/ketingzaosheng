%%变量声明
% dc_component(j)  不同电压等级直流分量均值，其中j对应不同的文件
% dc_component0(j) 背景噪声直流分量均值
%用于计算三分之一倍频程的中心频率和带宽从20Hz-20KHz
clear;
clc; 
close all;
%读取excel文件
[A_third_octave_info]=xlsread('20_150k的三分之一倍频程');
%  [A_third_octave_info]=xlsread('20_150k的等宽频程5K');
center_freq=A_third_octave_info(1,:)';
NN=length(center_freq);

%% ================ 确定文件路径 ====================
% 修改文件路径，输入背景干扰数据与电晕电流测量数据

% path0='E:\特高压直流输电\电晕笼数据\0\';     %背景干扰
% path='E:\特高压直流输电\电晕笼数据\1100\'; 
% path0='C:\Users\m1881\Desktop\电晕电流可听噪声数据\数据\';     %背景干扰
% path0='E:\特高压直流输电\实验线段数据\20120712_20m+\0\';
% path0='D:\开题\shuju\20120531+\0\';
 path0='D:\开题\shuju\20120712_20m+\0\';

% path='C:\Users\m1881\Desktop\电晕电流可听噪声数据\数据\';
% path='E:\特高压直流输电\实验线段数据\20120518+\1000\'; 
% path='E:\特高压直流输电\实验线段数据\20130531xiawu\1000\';
% path='E:\特高压直流输电\实验线段数据\20120712_20m+\1000\';
path='D:\开题\shuju\20130531xiawu\1000\';
% path='D:\开题\shuju\20120712_20m-\-700\';

size_n=1024*1024*2;                   % 读取长度
% 修改文件名
filenames=dir([path,'*.txt']);    %列出当前目录下符合正则表达式的文件夹和文件
NN1=length(filenames);
% NN1=2;

sum3=zeros(NN,NN1);
current_T=[];

% for j=1:5
  for j=6
% for j=2:length(filenames) 
     filename=[path,filenames(j).name];
     [fid]=fopen(filename,'r+');
     data_7=fread(fid,7,'schar');
    %add(j,10)=data_7(1);
    disp(['接下来运行第',num2str(j),'组数据。还有',num2str(length(filenames)-j),'组数据'])      %在command 窗口里显示的情况

    switch data_7(1)
        case 0;        Fs8=1024;
        case 1;        Fs8=1024*2;  
        case 2;        Fs8=1024*4;
        case 3;        Fs8=1024*8;
        case 4;        Fs8=1024*16;
        case 5;        Fs8=1024*32;
        case 6;        Fs8=1024*64;
        case 7;        Fs8=1024*128;
        case 8;        Fs8=1024*256;
        case 9;        Fs8=1024*512;
        case 10;       Fs8=1024*1024;
        case 18;       Fs8=1024*1024*62.5;   data_voltage1=fread(fid,size_n+7,'schar');  data_voltage=data_voltage1(8:size_n+7);
        case 21;       Fs8=1024*1024*500;    data_voltage1=fread(fid,size_n+1,'schar',7);data_voltage=data_voltage1(2:size_n+1);  %这些数据都没有包含前7个数
    end

fclose(fid);
x(1)=data_7(1);x(2)=data_7(3);
clear data_voltage1;
%% 读取背景数据
% filename=strcat(path0,'2014-03-08-21_10_11.txt');      % 背景干扰62.5MHz的电晕笼
%  filename=strcat(path0,'20120518084929.txt');        %20120518+  背景干扰62.5MHz的
%  filename=strcat(path0,'20120518085138.txt');        %20120518-  背景干扰62.5MHz的
% filename=strcat(path0,'2012-05-31-08_56_01.txt');     %20120531 背景干扰62.5MHz的 531下午没有背景数据
filename=strcat(path0,'2012-07-12-12_25_56.txt');   %20120712_20m+背景干扰62.5MHz的
% filename=strcat(path0,'2012-07-12-12_28_46.txt');   %20120712_20m-背景干扰62.5MHz的
% filename=strcat(path0,'2012-07-12-12_25_56.txt');    %20120712+
%  filename=strcat(path0,'2012-07-12-12_26_38.txt');    %20120712-
% filename=strcat(path0,'2012-07-10-17_15_24.txt');    %20120710+
% filename=strcat(path0,'2012-07-10-17_13_37.txt');    %20120710-
[fid]=fopen(filename,'r+');
data_background_7=fread(fid,7,'schar');   %先读取了前7个
disp('数据量较大，计算机运行速度有限，可爱的你耐心等候哦~么么哒')  
if data_background_7(1)==18
    data_background1=fread(fid,size_n+7,'schar');   
    data_background=data_background1(8:size_n+7);

else if data_background_7(1)==21   %此处还是进行了判断，
    data_background1=fread(fid,size_n+1,'schar',7); 
    data_background=data_background1(2:size_n+1);
    end
end

fclose(fid);
clear data_background1;

%% ============== 提取采样率、存储深度 ===================
    switch data_7(1)
        case 0;        Fs8=1024;
        case 1;        Fs8=1024*2;  
        case 2;        Fs8=1024*4;
        case 3;        Fs8=1024*8;
        case 4;        Fs8=1024*16;
        case 5;        Fs8=1024*32;
        case 6;        Fs8=1024*64;
        case 7;        Fs8=1024*128;
        case 8;        Fs8=1024*256;
        case 9;        Fs8=1024*512;
        case 10;       Fs8=1024*1024;
        case 18;       Fs8=1024*1024*62.5;
        case 21;       Fs8=1024*1024*500;
    end

    switch data_background_7(1)
        case 0;        Fs0=1024;
        case 1;        Fs0=1024*2;
        case 2;        Fs0=1024*4;
        case 3;        Fs0=1024*8;
        case 4;        Fs0=1024*16;
        case 5;        Fs0=1024*32;
        case 6;        Fs0=1024*64;
        case 7;        Fs0=1024*128;
        case 8;        Fs0=1024*256;
        case 9;        Fs0=1024*512;
        case 10;       Fs0=1024*1024;
        case 18;       Fs0=1024*1024*62.5;
        case 21;       Fs0=1024*1024*500;
    end

%     if Fs8~=Fs0
        Fs=1024*1024*62.5;
%     else if Fs8==Fs0
%         Fs=Fs8;
%         end
%     end

    clear Fs0; clear Fs8;

%量程提取
    switch data_7(3)
        case 0;        voltage8=0.1;
        case 1;        voltage8=0.2;
        case 2;        voltage8=0.5;
        case 3;        voltage8=1;
        case 4;        voltage8=2;
        case 5;        voltage8=5;
        case 6;        voltage8=10;
        case 7;        voltage8=20;
    end

    switch data_background_7(3)
        case 0;        voltage0=0.1;
        case 1;        voltage0=0.2;
        case 2;        voltage0=0.5;
        case 3;        voltage0=1;
        case 4;        voltage0=2;
        case 5;        voltage0=5;
        case 6;        voltage0=10;
        case 7;        voltage0=20;
    end

if length(data_background)~=length(data_voltage)
    error('错误！数据长度不一致！')
end
clear data_7;
clear data_background_7;

%% =============== 计算数据 ==================
%将电压数据转变为电流
N=size_n;
T=(N-1)/Fs;
t=0:1/Fs:T;

X8=10*voltage8*data_voltage(1:size_n)/127; 
X0=10*voltage0*data_background(1:size_n)/127; 

clear data_voltage; clear data_background; 
clear voltage8; clear voltage0;

% figure(1); plot(t*1000,X0); grid on;
% xlabel('时间(ms)');  ylabel('电流(mA)');
% xlim([0,1000*T])
% title('背景测量数据时域波形')
% figure(2); plot(t*1000,X8); grid on;
% xlabel('时间(ms)'); ylabel('电晕电流(mA)');
% xlim([0,1000*T])
% title('电晕电流测量数据时域波形')

%傅里叶变换
Nfft=2^nextpow2(N);
nfft=1:1:Nfft/2+1;
Y8=fft(X8,Nfft); 
Y0=fft(X0,Nfft);


dc_component(j)=Y8(1)/N;
dc_component0(j)=Y0(1)/N;    %第一个点就是直流分量


Y8(1)=0; Y0(1)=0;    % 去除直流分量
f=Fs/2*linspace(0,1,Nfft/2+1);
abs_Y8=2*abs(Y8(nfft))/N;
abs_Y0=2*abs(Y0(nfft))/N;

clear X0;

% figure(3); subplot(211);        % 背景信号的单边幅值频谱
% plot(f(1:1300),(abs_Y0(1:1300))); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('背景信号的单边幅值频谱')
% xlim([10,13200]);
% 
% figure(3); subplot(212);        % 加压信号的单边幅值频谱
%% 绘图
% figure(3)
% plot(f(1:1048577),(abs_Y8(1:1048577))); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('加压信号的单边幅值频谱')
% xlim([10,1048577]);



%% ============= 工频谐波处理 =================
N50=ceil(1200*Nfft/Fs);     % 截止频率1200Hz
w=zeros(20,1);
i=1;
for n=2:1:N50-1
    if abs_Y8(n)>abs_Y8(n-1)&&abs_Y8(n)>abs_Y8(n+1)   %若第n点的幅值大于两端的幅值
        w(i)=n-1;                                     %用于
        i=i+1;
    end
end
w(w(:)==0)=[];  % 把为0的初始值都剔除了
n_w=length(w);

clear abs_Y8; clear abs_Y0; clear i;

fw=w(1)*Fs/Nfft;      
ww=40/Fs*pi;                  %角度 因为零点是配置在了50Hz处，所以极点不能配置在50Hz处，就配置在50Hz附近，
r=sqrt((1-sin(ww)/cos(ww)));  %几乎等于1
wp=pi*fw/Fs;                  %角度值？
[B,A]=Notch03(wp,r);   % 计算陷波滤波器    r应该是极点的半径   得到的B A是 滤波器的系数
[H]=freqz(B,A,Nfft);   %然后就求了滤波器的频率特性？

 clear B; clear A; clear fw; clear N50;

W=ones(Nfft,1);
for n=1:1:n_w
    for i=2:1:Nfft/10
        W(i+w(n)-w(1))=H(i);
        W(Nfft-i+2)=W(i);
    end
    Y8=Y8.*W;    %作用陷波滤波器
end

clear W; clear H; clear w;

sum_08=sum(X8);
ave=sum_08/N;  
y=real(ifft(Y8,Nfft))+ave;

%% 绘制工频谐波处理后的信号
% figure(4); plot(t*1000,y); grid on;
% xlabel('时间(ms)');  ylabel('电流(mA)');
% axis([0,1000*T,-200,200])
% title('去除工频谐波后时域波形')
%% 对陷波以后的结果进行插值运算
for i=1:Nfft
    ff(i)=(i-1)*Fs/(Nfft+1);
end            %计算所有频率点的频率值

%计算每一点的幅值
nfft_1=1:1:Nfft;
abs_Y8_all=2*abs(Y8(nfft_1))/length(Y8);
abs_Y8=2*abs(Y8(nfft))/length(Y8);
point_num=length(abs_Y8_all);
abs_Y8_all_1=abs_Y8_all;
for i=2:point_num-1
    if abs_Y8_all_1(i)<1e-4
%        abs_Y8_all_1(i)=((ff(i)-ff(i-1))/ (ff(i+1)-ff(i-1)))*( abs_Y8_all(i+1)- abs_Y8_all(i-1))  +abs_Y8_all(i-1);
         Y8(i)=((ff(i)-ff(i-1))/ (ff(i+1)-ff(i-1)))*( Y8(i+1)- Y8(i-1))  +Y8(i-1);
    end
end

% abs_Y8_1=2*abs(Y8(nfft))/length(Y8);
%% 频域绘图
% figure(5)
% plot(f(1:1300),abs(abs_Y8(1:1300)));
% hold on 
% plot(f(1:1300),abs(abs_Y8_1(1:1300)),'r');
% grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('加压信号的单边幅值频谱')
% xlim([10,1300]);
% clear abs_Y8;
% clear y;
%% ============= 互相关函数运算 =================
Z=(Y8.*conj(Y0))/Nfft;
Thn=ceil(200000*Nfft/Fs);
thn=Thn:1:Nfft/2+1;  %低于200k的频谱不做处理

clear Y0;

% 动态调节阈值
r=1;
ncorr=zeros(20,1);   %记录被去除的频率点
% =============处理次数=================
for n=1:1:7                
    abs_Z=2*abs(Z(nfft))/Nfft;
    ths=1/2*max(abs_Z(thn));
    n_corr=find(abs_Z(thn)>=ths)+Thn-1;
    for i=1:1:length(n_corr)
        ncorr(r)=n_corr(i);
        r=r+1;
        n0=ceil(5000/n*Nfft/Fs);
        Z(n_corr(i)-n0:n_corr(i)+n0)=0;       %它是把某个范围内的数均置零了
%         Y8(n_corr(i)-n0:n_corr(i)+n0)=0;
%         Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=0;  %互相关运算去除对称部分的分量
%         Y8(n_corr(i)-n0:n_corr(i)+n0)=1/2*(Y8(n_corr(i)-n0-1)+Y8(n_corr(i)+n0+1));  %用两边的值的平均值代替
% %% 用平均值代替       
%         bx=f(n_corr(i)+n0+1);ax=f(n_corr(i)-n0-1);cx=f(n_corr(i));
%         by=Y8(n_corr(i)+n0+1);ay=Y8(n_corr(i)-n0-1);
%         Y8(n_corr(i)-n0:n_corr(i)+n0)=(cx-ax)*(by-ay)/(bx-ax)+ay;  %用两边的值的平均值代替
%         Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=fliplr(Y8(n_corr(i)-n0:n_corr(i)+n0));
        
%% 用最小二乘拟合
        kk=80;
        x_round_L=f((n_corr(i)-n0-kk):(n_corr(i)-n0-1));
        x_round_R=f((n_corr(i)+n0+1):(n_corr(i)+n0+kk));
        x_round=[x_round_L,x_round_R]';   %这是行向量
        Y8_round_L=Y8((n_corr(i)-n0-kk):(n_corr(i)-n0-1));
        Y8_round_R=Y8((n_corr(i)+n0+1):(n_corr(i)+n0+kk));
        Y8_round=[ Y8_round_L;Y8_round_R];
        XS0=polyfit(x_round,Y8_round,2);  %  用三阶去拟合
        for p=(n_corr(i)-n0):(n_corr(i)+n0)
            Y8(p)=1*polyval(XS0,f(p)); 
        end
        Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=fliplr(Y8(n_corr(i)-n0:n_corr(i)+n0));  
        
    end
    
end

clear Z; clear abs_Z; clear n0; clear Thn; 
clear n_corr; clear thn; clear i; clear r;

% ncorr*Fs/Nfft       %即为去除的窄带干扰模拟频率值
y8=real(ifft(Y8,Nfft))+ave;
abs_Y8=2*abs(Y8(nfft))/Nfft;   %单边幅值
clear Y8; clear ave; 
clear abs_Y8;

%% ================ 小波降噪 ======================
lev=7;
denoise_result=wden(y8,'rigrsure','s','one',lev,'db4');


figure(8); plot(t*1000,y8); grid on;
xlabel('时间(ms)'); ylabel('电晕电流(mA)');
xlim([0,1000*T])
title('电晕电流测量数据时域波形')

%% ================ 数据分析 ====================== 
%求平均值
denoise_result_mean=mean(denoise_result);
current_T=[current_T,denoise_result_mean];
% Y8=fft(y8,Nfft);
% Y8(1)=0;    % 去除直流分量
% abs_Y8=2*abs(Y8(nfft))/N;            %1024*1024个点
% abs_Y8(1)=1e-6; 


% [ Yc(:,j) ] = fs_same_th( A_third_octave_info,denoise_result,Fs,1000);
 [ Yc(:,j) ] = fs_same_th( A_third_octave_info,denoise_result,Fs,2500);   %这里用的也是没有经过小波降噪的信号
%  [ Yc(:,j) ] = one_3th( A_third_octave_info,denoise_result,Fs);   %这里用的也是没有经过小波降噪的信号
end

current_T=current_T';

figure (6)
plot(center_freq,Yc(:,2:j),'*-')
xlabel('中心频率点')
title('700KV各个曲线')
grid on
% 
% for i=1:NN
%     AVE(i)=mean(Yc(i,:));
% end
% figure (7)
% plot(center_freq,AVE,'*-')
% xlabel('中心频率点')
% title('0KV各个曲线平均')
% grid on
% AVE=AVE';
% 
dc_component=dc_component';
