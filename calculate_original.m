
clear; clc; 
close all;

gmax=[20.32 24.38 28.45 32.51 36.8 44.7];
gmax712=[20.2 23.1 26 28.87];%700-1000
gmax531xia=[20.2 23.1 26 28.9];%700-1000
gmax531=[20.1 23 25.87 28.75];%700-1000
gmax518=[22.04 27.55];%700-1000
r=1.81;%单位为cm
D=30.5;
n=6;
E1=38+1.6*(gmax712(4)-24)+46*log10(r)+5*log10(n)+33*log10(20/D);

f_f=0.05:0.05:5;
for i=1:length(f_f)
EF(i)=5*[1-2*(log10(10*f_f(i)))^2]+1.6*(gmax(6)-24)+46*log10(r)+5*log10(n)+4;
end

%% ================ 确定文件路径 ====================
% 修改文件路径，输入背景干扰数据与电晕电流测量数据
path0='E:\特高压直流输电\电晕笼数据\0\';     %背景干扰
path='E:\特高压直流输电\电晕笼数据\800\';     %电晕电流
% path0='E:\特高压直流输电\起晕电压判断\数据处理\example\exampleb\';
% path0='E:\特高压直流输电\实验线段数据\20120518+\0\'; 
% path0='E:\特高压直流输电\实验线段数据\20120712+\0\';     %背景干扰
% path0='E:\特高压直流输电\实验线段数据\20120712_20m+\0\';
% path0='E:\特高压直流输电\实验线段数据\20120531+\0\';
% path0='E:\特高压直流输电\西藏数据\0\';   %西藏背景
% path0='E:\特高压直流输电\实验线段数据\20120518+\0\';       %6*300的背景



% path='E:\特高压直流输电\起晕电压判断\数据处理\example\example\';
% path='E:\特高压直流输电\实验线段数据\20120712+\1000\';   %与20120712_20m+ 的背景信号一致
% path='E:\特高压直流输电\实验线段数据\20120518+\1000\'; 
% path='E:\特高压直流输电\实验线段数据\20130531xiawu\900\';
% path='E:\特高压直流输电\实验线段数据\20120712_20m+\900\';
% path='E:\特高压直流输电\实验线段数据\20120531+\1000\';
% path='E:\特高压直流输电\西藏数据\800\';
% path='E:\特高压直流输电\导线6_300数据\1000\';

size=1024*1024*8+7;      % 读取长度

% 修改文件名
%filename=strcat(path,'11.txt');  %电晕电流
filenames=dir([path,'*.txt']);
% for j=1:length(filenames) 
for j= 3
% for j= 2
 filename=[path,filenames(j).name];
 [fid]=fopen(filename,'r+');
data_7=fread(fid,7,'schar');
% add(j,10)=data_7(1);
switch data_7(1)
    case 0
        Fs8=1024;
	case 1
        Fs8=1024*2;  
	case 2
        Fs8=1024*4;
	case 3
        Fs8=1024*8;
	case 4
        Fs8=1024*16;
	case 5
        Fs8=1024*32;
	case 6
        Fs8=1024*64;
	case 7
        Fs8=1024*128;
	case 8
        Fs8=1024*256;
	case 9
        Fs8=1024*512;
    case 10
        Fs8=1024*1024;
    case 18
        Fs8=1024*1024*62.5;data_voltage1=fread(fid,size-7,'schar');
    case 21
        Fs8=1024*1024*500;data_voltage1=fread(fid,size-7,'schar',7);
end

fclose(fid);
data_voltage=cat(1,data_7,data_voltage1);
x(1)=data_7(1);x(2)=data_7(3);
clear data_7;clear data_voltage1;

% filename=strcat(path0,'20120518084929.txt');        %20120518+  背景干扰62.5MHz的
% filename=strcat(path0,'2014-03-08-21_10_11.txt');   %exampleb里面的
% filename=strcat(path0,'2012-05-31-08_51_32.txt');     %20120531 背景干扰62.5MHz的 531下午没有背景数据
% filename=strcat(path0,'2012-07-12-12_25_56.txt');   %20120712_20m+背景干扰62.5MHz的
% filename=strcat(path0,'2012-07-12-12_25_56.txt');    %20120712+
%  filename=strcat(path0,'2013-05-05-09_30_18.txt');     %西藏背景数据
filename=strcat(path0,'2014-03-08-21_10_11.txt');    %电晕笼的背景数据

[fid]=fopen(filename,'r+');
data_background=fread(fid,7,'schar');x(3)=data_background(1);x(4)=data_background(3);

data_background1=data_background(1);
data_background3=data_background(3);
if data_background(1)==18
    data_background=fread(fid,size,'schar'); 

else if data_background(1)==21
    data_background=fread(fid,size,'schar',7);
    end
end

     
    fclose(fid);

%% ============== 提取采样率、存储深度 ===================

switch data_voltage(1)
    case 0
        Fs8=1024;
	case 1
        Fs8=1024*2;  
	case 2
        Fs8=1024*4;
	case 3
        Fs8=1024*8;
	case 4
        Fs8=1024*16;
	case 5
        Fs8=1024*32;
	case 6
        Fs8=1024*64;
	case 7
        Fs8=1024*128;
	case 8
        Fs8=1024*256;
	case 9
        Fs8=1024*512;
    case 10
        Fs8=1024*1024;
    case 18
        Fs8=1024*1024*62.5;
    case 21
        Fs8=1024*1024*500;
end

switch data_background1
    case 0
        Fs0=1024;
	case 1
        Fs0=1024*2;
	case 2
        Fs0=1024*4;
	case 3
        Fs0=1024*8;
	case 4
        Fs0=1024*16;
	case 5
        Fs0=1024*32;
	case 6
        Fs0=1024*64;
	case 7
        Fs0=1024*128;
	case 8
        Fs0=1024*256;
	case 9
        Fs0=1024*512;
    case 10
        Fs0=1024*1024;
    case 18
        Fs0=1024*1024*62.5;
    case 21
        Fs0=1024*1024*500;
end

% if Fs8~=Fs0
    Fs=1024*1024*62.5;
% else if Fs8==Fs0
%     Fs=Fs8;
%     end
% end

clear Fs0; clear Fs8;

%量程提取
switch data_voltage(3)
    case 0 
        voltage8=0.1;
	case 1
        voltage8=0.2;
	case 2
        voltage8=0.5;
	case 3
        voltage8=1;
	case 4
        voltage8=2;
	case 5
        voltage8=5;
	case 6
        voltage8=10;
	case 7
        voltage8=20;
end

switch data_background3
    case 0 
        voltage0=0.1;
	case 1
        voltage0=0.2;
	case 2
        voltage0=0.5;
	case 3
        voltage0=1;
	case 4
        voltage0=2;
	case 5
        voltage0=5;
	case 6
        voltage0=10;
	case 7
        voltage0=20;
end

if length(data_background)~=length(data_voltage)
    error('错误！数据长度不一致！')
end

%% =============== 计算数据 ==================
%将电压数据转变为电流
N=size-7;
T=(N-1)/Fs;
t=0:1/Fs:T;

X8=10*voltage8*data_voltage(8:size)/127; 
X0=10*voltage0*data_background(8:size)/127; 

clear data_voltage; clear data_background; 
clear voltage8; clear voltage0;

%绘制去噪前的电晕电流波形 
% figure(1); 
% plot(t*1000,X0); 
% % grid on;
% xlabel('时间(ms)');  ylabel('电流(mA)');
% xlim([0,1000*T])
% title('背景测量数据时域波形')
figure(1); 
h = gca; % 获取当前绘图坐标的指针
size_font = 14;
set(h,'FontSize',size_font);
plot(t*1000,X8,'b');
xlabel('Time(ms)','FontSize',size_font); 
ylabel('Corona Current(mA)','FontSize',size_font);
% ylim([-10,10])
% xlim([0,30])
axis([0,30,-10,10])
% title('电晕电流测量数据时域波形')

%傅里叶变换
Nfft=2^nextpow2(N);
nfft=1:1:Nfft/2+1;
Y8=fft(X8,Nfft); Y0=fft(X0,Nfft);
Y8(1)=0; Y0(1)=0;    % 去除直流分量
f=Fs/2*linspace(0,1,Nfft/2+1);
abs_Y8=2*abs(Y8(nfft))/N;
abs_Y0=2*abs(Y0(nfft))/N;

% clear X0;
% 
% figure(3); subplot(211);        % 背景信号的单边幅值频谱
% plot(f,abs_Y0); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('背景信号的单边幅值频谱')
% xlim([10,10^7]);
% figure(3); subplot(212);        % 加压信号的单边幅值频谱
% plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('加压信号的单边幅值频谱')
% xlim([10,10^7]);



%% ============= 工频谐波处理 =================
N50=ceil(1200*Nfft/Fs);     % 截止频率1200Hz
w=zeros(20,1);
i=1;
for n=2:1:N50-1
    if abs_Y8(n)>abs_Y8(n-1)&&abs_Y8(n)>abs_Y8(n+1)
        w(i)=n-1;
        i=i+1;
    end
end
w(w(:)==0)=[];  % 中心频率点
n_w=length(w);

clear abs_Y8; clear abs_Y0; clear i;

fw=w(1)*Fs/Nfft;      
ww=40/Fs*pi;             % 带宽
r=sqrt((1-sin(ww)/cos(ww)));
wp=pi*fw/Fs;  
[B,A]=Notch03(wp,r);   % 计算陷波滤波器 
[H]=freqz(B,A,Nfft);

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

% figure(4); plot(t*1000,y); grid on;
% xlabel('时间(ms)');  ylabel('电流(mA)');
% axis([0,1000*T,-150,150])
% title('去除工频谐波后时域波形')

clear y;


%% ============= 互相关函数运算 =================
Z=(Y8.*conj(Y0))/Nfft;
Thn=ceil(200000*Nfft/Fs);
thn=Thn:1:Nfft/2+1;  %低于200k的频谱不做处理

clear Y0;

% 动态调节阈值
r=1;
ncorr=zeros(20,1);   %记录被去除的频率点
for n=1:1:7                % !!!!!!!!!!!!!!!处理次数!!!!!!!!!!!!!!!!!
    abs_Z=2*abs(Z(nfft))/Nfft;
    ths=1/2*max(abs_Z(thn));
    n_corr=find(abs_Z(thn)>=ths)+Thn-1;
    for i=1:1:length(n_corr)
        ncorr(r)=n_corr(i);
        r=r+1;
        n0=ceil(5000/n*Nfft/Fs);
        Z(n_corr(i)-n0:n_corr(i)+n0)=0;
        Y8(n_corr(i)-n0:n_corr(i)+n0)=0;
    end
end

clear Z; clear abs_Z; clear n0; clear Thn; 
clear n_corr; clear thn; clear i; clear r;

% ncorr*Fs/Nfft       %即为去除的窄带干扰模拟频率值
y8=real(ifft(Y8,Nfft))+ave;
abs_Y8=2*abs(Y8(nfft))/Nfft;

clear Y8; %clear ave; 



clear abs_Y8;


%% ================ 小波降噪 ======================
lev=7;
denoise_result=wden(y8,'rigrsure','s','one',lev,'db4'); 

figure(2); 
% hold on 
h = gca; % 获取当前绘图坐标的指针
set(h,'FontSize',size_font);
size_font = 14;
plot(t*1000,denoise_result,'r','linewidth',0.1);
xlabel('Time(ms)','FontSize',size_font);
ylabel('Corona Current(mA)','FontSize',size_font);
% ylim([-10,10])
% xlim([0,30])
axis([0,30,-10,10])
legend('降噪前','降噪后')
% title('去噪后的时域波形')



%% ================ 数据分析 ======================

% Y8=fft(y8,Nfft);
% Y8(1)=0;    % 去除直流分量
% abs_Y8=2*abs(Y8(nfft))/N;
% abs_Y8(1)=1e-6;  
% plot_Pxx=20*log10(abs_Y8/2e-5);
% 
% 
% for nfft=1:1:Nfft/2+1
%     if plot_Pxx(nfft)<0
%        plot_Pxx_abs(nfft)=0;
%     else plot_Pxx_abs(nfft)=plot_Pxx(nfft);
%     end
%     
%     
% end


% figure(3);        % 加压信号的单边幅值频谱
% plot(f,plot_Pxx_abs); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('加压信号的单边幅值频谱')
% xlim([1,1000*20]);


%正常幅值
% ff=0.15:0.05:30;
% for ii=1:length(ff)
%  pingjunzhi(ii)=mean(abs_Y8(round((140+50*ii)*1000/Fs*Nfft):round((160+50*ii)*1000)/Fs*Nfft));
% end
% plot(ff,pingjunzhi,'r');
%  xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
%  title('小波降噪后的单边幅值频谱')
 
%对数幅值

% ffff=[31.5 40 50 63	80 100 125 160 200 250 315	400	500	630	800	1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000];
% 
% for ii=1:length(ffff)
%  pingjunzhi1abs(ii)=mean(plot_Pxx_abs(round(ffff(ii)/2^(1/6)/Fs*Nfft):round(ffff(ii)*2^(1/6)/Fs*Nfft)));
% end




% fff=0.05:0.05:5;
% for ii=1:length(fff)
%  pingjunzhi(ii)=mean(plot_Pxx(round((1+50*(ii-1))*1000/Fs*Nfft):round((100+50*(ii-1))*1000)/Fs*Nfft));
% end


% figure(1);
% plot(fff,pingjunzhi,'b');
% xlabel('频率(MHz)'); ylabel('对数幅值|Y(f)|（dB）');
% title('电晕笼数据 900kv 对数频域幅值预测值与实际值对比图 绝对误差dB');
% grid on;hold on;
% 
% plot(f_f,EF,'r');
% 
% legend('实际值','预测值');
% 
% 
% wucha=abs(EF-pingjunzhi);
% wucha1=abs(EF-pingjunzhi)/abs(pingjunzhi);
% xxx(j)=mean(wucha);









% f1=[ 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000];
% 
% shengji1100=[61.28345712	62.68491524	64.49003947	65.99801444	66.89919123	67.59959467	68.09978233	68.59992625	68.39993576	68.19995947	68.09997558	67.69998231	66.79998804	66.1999921	63.29999005	59.99999191];
% shengji1000=[59.37435126	61.07817899	62.48420303	63.89677934	64.69865769	65.49934262	65.99964698	66.29987475	66.39989819	65.99993274	66.09996129	65.59997131	64.59998015	63.89998658	60.9999831	54.6999726];
% shengji900=[	56.34867273	57.95533001	59.46842363	60.69326836	61.19699438	62.99883092	62.49920964	62.59970637	62.49975008	62.3998459	62.39990926	61.79993117	60.79995238	60.19996854	57.39996129	51.09993722];
% shengji800=[	50.2968091	52.23557046	53.26990128	54.06914496	53.58267592	53.79026646	54.69523544	54.69818918	54.79852817	54.89913338	54.69946567	54.1996039	53.29973221	52.49981473	49.39975577	42.69956568];
% shengji700=[	47.52285795	48.95738822	49.495984	49.30830173	47.93676724	48.4669298	48.78143426	48.89311144	48.89427111	48.99662749	49.09805964	48.99868825	47.99909254	47.09935758	44.19919123	37.59859442];
% shengji600=[42.23988945	44.00818236	44.99169983	45.27135246	43.31937228	42.77870725	43.63963306	43.7776694	43.88185775	43.88907729	44.19400093	44.69646849	44.19782283	43.19842289	39.69772019	32.59555362];
% shengji500=[	29.52783738	31.72235043	34.90433774	39.49169983	36.85059633	34.78725481	35.84907659	35.24310635	34.75375671	34.19935157	34.13963306	34.56372589	34.47963887	33.98686467	30.38055732	23.66536502];
% 
% alpha=[	3.7299	5.7219	8.5764	12.8254	19.7889	28.7331	40.6027	55.9366
% 74.2386	92.2748	110.5456	128.2559	144.0938	160.6596	182.7201	209.92];
% 
% 
% shengji=[33.8 30.7 35.9 40.8 45.6 49.9 53.6 56.5 59.1 60.4 61.3 62.7 64.5 66 66.9 67.6 68.1 68.6 68.4 68.2 68.1 67.7 66.8 66.2 63.3 60];
% plot(f1,shengji1100+alpha*0.005,'r');
% 
% legend('对数幅值','可听噪声分贝数');



% 
% junzhi500db(j)=mean(plot_Pxx_abs(round(490*1000/Fs*Nfft):round(510*1000/Fs*Nfft)));
% junzhi500(j)=10e5*mean(abs_Y8(round(490*1000/Fs*Nfft):round(510*1000/Fs*Nfft)));































% pingjunzhi(1)=mean(abs_Y8(round(1/Fs*Nfft+1):round(10*1000)/Fs*Nfft));
% pingjunzhi(2)=mean(abs_Y8(round(10*1000/Fs*Nfft+1):round(60*1000)/Fs*Nfft));
% pingjunzhi(3)=mean(abs_Y8(round(60*1000/Fs*Nfft+1):round(1000*1000)/Fs*Nfft));
% pingjunzhi(5)=mean(abs_Y8(round(1000*1000/Fs*Nfft+1):round(1500*1000)/Fs*Nfft));
% pingjunzhi(6)=1000*mean(abs_Y8(round(1500*1000/Fs*Nfft+1):round(10000*1000)/Fs*Nfft));
% pingjunzhi(7)=1000*mean(abs_Y8(round(6000*1000/Fs*Nfft+1):round(7000*1000)/Fs*Nfft));
% pingjunzhi(8)=1000*mean(abs_Y8(round(7500*1000/Fs*Nfft+1):round(8000*1000)/Fs*Nfft));
% pingjunzhi(4)=mean(abs_Y8(round(500*1000/Fs*Nfft+1):round(1000*1000)/Fs*Nfft));
% pingjunzhi
% 
% pingfanghe(1)=sum(abs_Y8(round(1/Fs*Nfft+1):round(150*1000)/Fs*Nfft).^2);
% pingfanghe(2)=sum(abs_Y8(round(1/Fs*Nfft+1):round(500*1000)/Fs*Nfft).^2);
% pingfanghe(3)=sum(abs_Y8(round(150*1000/Fs*Nfft):round(1000*10000)/Fs*Nfft).^2);
% pingfanghe(4)=sum(abs_Y8(round(500*1000/Fs*Nfft+1):round(1000*2000)/Fs*Nfft).^2);
% pingfanghe(5)=sum(abs_Y8(round(1/Fs*Nfft+1):round(1000*10000)/Fs*Nfft).^2);
% pingfanghe
% junzhi=mean(abs_Y8(round(490*1000/Fs*Nfft):round(500*1000)/Fs*Nfft))
% f=0.5;
% g=44.7;
% r=1.81;%单位为cm
% D=30.5;
% n=6;
% 
% e=38+1.6*(g-24)+46*log10(r)+5*log10(n)+5*[1-2*(log10(10*f)^2)]+33*log10(20/D)



% for i=1:5
% figure(i); plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('小波降噪后的单边幅值频谱')
%  xlim([100*(i-1)*1000,100*i*1000]);
% end
 
%  ylim([0,0.3])
%ylim([0,0.025]);
% for i=11:20
% figure(i); plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('小波降噪后的单边幅值频谱')
% xlim([1000*(i-10)*1000,1000*(i-9)*1000]);
% end
% 
% figure(21); plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('小波降噪后的单边幅值频谱')
% xlim([10,1000*1000]);
% %  ylim([0,0.002]);
% 
% figure(22); plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('小波降噪后的单边幅值频谱')
% xlim([500*1000,1000*1000]);
% 
% figure(23); plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('小波降噪后的单边幅值频谱')
% xlim([500*1000,2000*1000]);
% 
% figure(24); plot(f,abs_Y8); grid on;
% xlabel('频率(Hz)'); ylabel('幅值|Y(f)|');
% title('小波降噪后的单边幅值频谱')
% xlim([2000*1000,10000*1000]);



















% totalmean(1)=10^5*mean(abs_Y8(round(1000*480/Fs*Nfft):round(1000*520)/Fs*Nfft));
% totalmean(2)=10^5*mean(abs_Y8(round(1000*450/Fs*Nfft):round(1000*550)/Fs*Nfft));
% totalmean(3)=10^5*mean(abs_Y8(round(1000*400/Fs*Nfft):round(1000*600)/Fs*Nfft));
% totalmean(4)=10^5*mean(abs_Y8(round(1000*350/Fs*Nfft):round(1000*650)/Fs*Nfft));
% totalmean(5)=10^5*mean(abs_Y8(round(1000*300/Fs*Nfft):round(1000*700)/Fs*Nfft));
% totalmean(6)=10^5*mean(abs_Y8(round(1000*450/Fs*Nfft):round(1000*500)/Fs*Nfft));
% totalmean(7)=10^5*mean(abs_Y8(round(1000*500/Fs*Nfft):round(1000*550)/Fs*Nfft));
% totalmean(8)=10^5*mean(abs_Y8(round(1000*400/Fs*Nfft):round(1000*500)/Fs*Nfft));
% totalmean(9)=10^5*mean(abs_Y8(round(1000*500/Fs*Nfft):round(1000*600)/Fs*Nfft));
% totalmean(10)=10^5*mean(abs_Y8(round(1000*350/Fs*Nfft):round(1000*500)/Fs*Nfft));
% totalmean(11)=10^5*mean(abs_Y8(round(1000*500/Fs*Nfft):round(1000*650)/Fs*Nfft));
% totalmean(12)=10^5*mean(abs_Y8(round(1000*500/Fs*Nfft):round(1000*1000)/Fs*Nfft));
% totalmean(14)=10^5*mean(abs_Y8(round(1000*3000/Fs*Nfft):round(1000*8000)/Fs*Nfft));
% totalmean
% j=1;
% 
% mn=round(1000*500/Fs*Nfft)
% for mn=round(1000*500/Fs*Nfft):round(1000*1000)/Fs*Nfft
%     if abs_Y8(mn)>3*10^-5*totalmean(12)
%         yao(j)=abs_Y8(mn);
%         j=j+1;
%     end
% end
% 
% totalmean(13)=10^5*mean(yao);
% 
% 
% 
% ij=1
% for mn=round(1000*3000/Fs*Nfft):round(1000*8000)/Fs*Nfft
%     if abs_Y8(mn)>3*10^-5*totalmean(13)
%         yaoyao(ij)=abs_Y8(mn);
%         ij=ij+1;
%     end
% end
% 
% totalmean(15)=10^5*mean(yaoyao);
% totalmean
% j-1
% ij-1










% 
% for i=1:20
%     k(i)=round(i*1000/Fs*Nfft);
% end
% for i=1:8
%     add(j,i)=mean(plot_Pxx(k(i):k(20)));
% end
% add(j,9)=mean(plot_Pxx(k(10):k(20)));
% average(j)=mean(abs_Y8(round(1000*400/Fs*Nfft):round(1000*600)/Fs*Nfft));
% 
% 
% for i=1:501
%     ave1(i)=mean(abs_Y8(round(1000*20*(i-1)/Fs*Nfft+1):round(1000*20*i/Fs*Nfft)));
% end
% figure(4);
% plot((1:10000)*1000,ave);grid on; 
% %xlim([0,1000000]);
% [a,b]=max(ave)



end





