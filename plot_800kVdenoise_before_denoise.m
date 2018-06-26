%--------------------------------------------------------------------
% ��������
% dc_component(j)  ��ͬ��ѹ�ȼ�ֱ��������ֵ������j��Ӧ��ͬ���ļ�
% dc_component0(j) ��������ֱ��������ֵ
%���������ڻ���ȥ��ǰ��ĵ��ε���ͼ
%--------------------------------------------------------------------

clear; clc; 
close all;


width=440;%��ȣ�������
height=330;%�߶�
left=200;%����Ļ���½�ˮƽ����
bottem=100;%����Ļ���½Ǵ�ֱ����

%% ====���߱��泡ǿ=======  ������û�����ص�����
gmax518=[22.04 27.55];                     %800��1000
gmax531=[20.1 23 25.87 28.75];             %700-800-900-1000
gmax_long=[20.32 24.38 28.45 32.51 36.8 44.7];  %�������ĵ糡ǿ��
gmax712=[20.2 23.1 26 28.87];              %700-1000
gmax531xia=[20.2 23.1 26 28.9];            %700-800-900-1000
%% ���߲���
r=1.81;                           %��λΪcm  �Ե��߰뾶
D=30.5;                           %���������Ե��ߵľ���D
n=6;
Y8_ave=zeros(1048577,1);
%% ��Ҫ�޸ĵĲ���
Gmax=gmax518(2);
% d_Hz=0.025;
% f_f=0.15:d_Hz:4;            %���㷶Χ��0.15MHz-4Hz
%% ================ ȷ���ļ�·�� ====================
% �޸��ļ�·�������뱳��������������ε�����������
% path0='E:\�ظ�ѹֱ�����\����������\0\';     %��������
% path='E:\�ظ�ѹֱ�����\����������\1100\'; 
% path0='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120518+\0\';     %��������
path0='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120712_20m+\0\';
% path0='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120531+\0\';

% path='E:\�ظ�ѹֱ�����\����������\1100\'; 
% path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120518+\1000\'; 
% path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20130531xiawu\1000\';
path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120712_20m+\1000\';
% path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120531+\700\';

voltage_list = 1100;
if voltage_list == 500
    sheet = 1;
end
if voltage_list == 600
    sheet = 2;
end
if voltage_list == 700
    sheet = 3;
end
if voltage_list == 800
    sheet = 4;
end
if voltage_list == 900
    sheet = 5;
end
if voltage_list == 1100
    sheet = 6;
end
    
    

size_n=1024*1024*8;      % ��ȡ����
% �޸��ļ���
%filename=strcat(path,'11.txt');  %���ε���
filenames=dir([path,'*.txt']);
% for j=1:length(filenames) 
for j=1:1
 filename=[path,filenames(j).name];
 [fid]=fopen(filename,'r+');
 data_7=fread(fid,7,'schar');

disp(['���������е�',num2str(j),'�����ݡ�����',num2str(length(filenames)-j),'������'])
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
    case 21;       Fs8=1024*1024*500;    data_voltage1=fread(fid,size_n+1,'schar',7);data_voltage=data_voltage1(2:size_n+1);  %��Щ���ݶ�û�а���ǰ7����
end

fclose(fid);
clear data_voltage1;
%% ��ȡ��������
%ע�ⱳ�����ݵ��޸�
% filename=strcat(path0,'2014-03-08-21_10_11.txt');      % ��������62.5MHz�ĵ�����
% filename=strcat(path0,'20120518084929.txt');        %20120518+  ��������62.5MHz��
% filename=strcat(path0,'2012-05-31-08_51_32.txt');     %20120531 ��������62.5MHz�� 531����û�б�������
filename=strcat(path0,'2012-07-12-12_25_56.txt');   %20120712_20m+��������62.5MHz��
[fid]=fopen(filename,'r+');
data_background_7=fread(fid,7,'schar');   %�ȶ�ȡ��ǰ7��
disp('�������ϴ󣬼���������ٶ����ޣ��ɰ��������ĵȺ�Ŷ~ôô��')  
if data_background_7(1)==18
    data_background1=fread(fid,size_n+7,'schar');   
    data_background=data_background1(8:size_n+7);

else if data_background_7(1)==21   %�˴����ǽ������жϣ�
    data_background1=fread(fid,size_n+1,'schar',7); 
    data_background=data_background1(2:size_n+1);
    end
end

fclose(fid);
clear data_background1;
%% ============== ��ȡ�����ʡ��洢��� ===================

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

% if Fs8~=Fs0
    Fs=1024*1024*62.5;
% else if Fs8==Fs0
%     Fs=Fs8;
%     end
% end

clear Fs0; clear Fs8;

%% ������ȡ
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
    error('�������ݳ��Ȳ�һ�£�')
end
clear data_7;
clear data_background_7;
%% =============== �������� ==================
%����ѹ����ת��Ϊ����
N=size_n;
T=(N-1)/Fs;
t=0:1/Fs:T;

X8=10*voltage8*data_voltage(1:size_n)/127; 
X0=10*voltage0*data_background(1:size_n)/127; 

clear data_voltage; clear data_background; 
clear voltage8; clear voltage0;


figure(1); 
size_font = 16;
h = gca; % ��ȡ��ǰ��ͼ�����ָ��
set(h,'FontSize',size_font);
plot(t*1000,X8,'k'); 
xlabel('ʱ��(ms)'); 
ylabel('���ε���(mA)');
% xlim([0,40])
axis([0,1000*0.04,-20,20])
% % title('���ε�����������ʱ����')
% % title('Time domain waveform of corona current data before denoising')
% set(gcf,'position',[left,bottem,width,height])

%% ����Ҷ�任
Nfft=2^nextpow2(N);
nfft=1:1:Nfft/2+1;
Y8=fft(X8,Nfft); 
Y0=fft(X0,Nfft);
zero_num1=find(abs(Y8)==0);

dc_component(j)=Y8(1)/N;
dc_component0(j)=Y0(1)/N;


Y8(1)=0; Y0(1)=0;    % ȥ��ֱ������
f=Fs/2*linspace(0,1,Nfft/2+1);
abs_Y8=2*abs(Y8(nfft))/N;
abs_Y0=2*abs(Y0(nfft))/N;

% clear X0;
% figure(3); subplot(211);        % �����źŵĵ��߷�ֵƵ��
% plot(f,abs_Y0); grid on;
% xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
% title('�����źŵĵ��߷�ֵƵ��')
% xlim([10,10^7]);

% figure(2);   % ��ѹ�źŵĵ��߷�ֵƵ��
% plot(f,abs_Y8); grid on;
% xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
% title('��ѹ�źŵ�ԭʼ���߷�ֵƵ��')




%% ============= ��Ƶг������ =================
N50=ceil(1200*Nfft/Fs);     % ��ֹƵ��1200Hz
w=zeros(20,1);
i=1;
for n=2:1:N50-1
    if abs_Y8(n)>abs_Y8(n-1)&&abs_Y8(n)>abs_Y8(n+1)
        w(i)=n-1;
        i=i+1;
    end
end
w(w(:)==0)=[];  % ����Ƶ�ʵ�
n_w=length(w);

clear abs_Y8; clear abs_Y0; clear i;

fw=w(1)*Fs/Nfft;      
ww=40/Fs*pi;             % ����
r=sqrt((1-sin(ww)/cos(ww)));
wp=pi*fw/Fs;  
[B,A]=Notch03(wp,r);   % �����ݲ��˲��� 
[H]=freqz(B,A,Nfft);

clear B; clear A; clear fw; clear N50;

W=ones(Nfft,1);
for n=1:1:n_w
    for i=2:1:Nfft/10
        W(i+w(n)-w(1))=H(i);
        W(Nfft-i+2)=W(i);
    end
    Y8=Y8.*W;    %�����ݲ��˲���
end
zero_num2=find(abs(Y8)==0);
clear W; clear H; clear w;

sum_08=sum(X8);
ave=sum_08/N;  
y=real(ifft(Y8,Nfft))+ave;

% figure(3); plot(t*1000,y); grid on;
% xlabel('ʱ��(ms)');  ylabel('����(mA)');
% axis([0,1000*T,-20,20])
% title('ȥ����Ƶг����ʱ����')

clear y;


%% ============= ����غ������� =================
Z=(Y8.*conj(Y0))/Nfft;      %��ʵ�ʵ����źŵ�Ƶ��ֵ��˱��������Ĺ���
Thn=ceil(200000*Nfft/Fs);  %��0.2MHz���ڵ��źŽ��д���
thn=Thn:1:Nfft/2+1;        

clear Y0;

% ��̬������ֵ
r=1;
ncorr=zeros(20,1);               %��¼��ȥ����Ƶ�ʵ�
for n=1:1:7                      %�������
    abs_Z=2*abs(Z(nfft))/Nfft;
    ths=1/2*max(abs_Z(thn));
    n_corr=find(abs_Z(thn)>=ths)+Thn-1;
    for i=1:1:length(n_corr)
        ncorr(r)=n_corr(i);
        r=r+1;
        n0=ceil(5000/n*Nfft/Fs);
        Z(n_corr(i)-n0:n_corr(i)+n0)=0;       %���ǰ�ĳ����Χ�ڵ�����������
%         Y8(n_corr(i)-n0:n_corr(i)+n0)=0;
%         Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=0;  %���������ȥ���ԳƲ��ֵķ���
%         Y8(n_corr(i)-n0:n_corr(i)+n0)=1/2*(Y8(n_corr(i)-n0-1)+Y8(n_corr(i)+n0+1));  %�����ߵ�ֵ��ƽ��ֵ����
% %% ��ƽ��ֵ����       
%         bx=f(n_corr(i)+n0+1);ax=f(n_corr(i)-n0-1);cx=f(n_corr(i));
%         by=Y8(n_corr(i)+n0+1);ay=Y8(n_corr(i)-n0-1);
%         Y8(n_corr(i)-n0:n_corr(i)+n0)=(cx-ax)*(by-ay)/(bx-ax)+ay;  %�����ߵ�ֵ��ƽ��ֵ����
%         Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=fliplr(Y8(n_corr(i)-n0:n_corr(i)+n0));
        
%% ����С�������
    kk=80;
    x_round_L=f((n_corr(i)-n0-kk):(n_corr(i)-n0-1));
    x_round_R=f((n_corr(i)+n0+1):(n_corr(i)+n0+kk));
    x_round=[x_round_L,x_round_R]';   %����������
    Y8_round_L=Y8((n_corr(i)-n0-kk):(n_corr(i)-n0-1));
    Y8_round_R=Y8((n_corr(i)+n0+1):(n_corr(i)+n0+kk));
    Y8_round=[ Y8_round_L;Y8_round_R];
    XS0=polyfit(x_round,Y8_round,1);  %  ������ȥ���
    for p=(n_corr(i)-n0):(n_corr(i)+n0)
        Y8(p)=polyval(XS0,f(p)); 
    end
    Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=fliplr(Y8(n_corr(i)-n0:n_corr(i)+n0));  
        
    end
end
% Y8(1)=0;
clear Z; clear abs_Z; clear n0; clear Thn; 
clear n_corr; clear thn; clear i; clear r;
f_dele=ncorr*Fs/Nfft;       %��Ϊȥ����խ������ģ��Ƶ��ֵ


%�����������֮���Ƶ��ֵ
abs_Y8=2*abs(Y8(nfft))/Nfft;    %��ȡ�˷�ֵ

% figure(4); 
% plot(f/1000000,20*log10(abs_Y8/2e-5),'k-'); 
% %  plot(f/1000000,abs_Y8,'k-'); 
% % plot(f,abs_Y8,'k-'); 
% grid on;
% xlabel('frequency(MHz)'); ylabel('amplitude|Y(f)|');
% xlim([0.15,4])
% title('�����������֮���Ƶ��ֵ')
% % xlim([20,20000])

y8=real(ifft(Y8,Nfft))+ave; 

lev=7;
denoise_result=wden(y8,'rigrsure','s','one',lev,'db4'); 


figure
size_font = 16;
h = gca; % ��ȡ��ǰ��ͼ�����ָ��
set(h,'FontSize',size_font);
plot(t*1000,denoise_result,'k'); 
xlabel('ʱ��(ms)'); ylabel('���ε���(mA)');
% xlim([0,1000*T])
axis([0,1000*0.04,-20,20])
% title('���ε�����������ʱ����')
% title('Time domain waveform of corona current data before denoising')
% set(gcf,'position',[left,bottem,width,height])


% figure; plot(t*1000,y8); grid on;
% xlabel('ʱ��(ms)');  ylabel('����(mA)');
% axis([0,1000*T,-20,20])
% title('ȥ����Ƶг����ʱ����')

%% ��ƽ���õ�������Ȼ����ȡƽ��ֵ
% power_f=abs_Y8.^2;  
% energy2k_20k(j)=mean(power_f(round(2*1000/Fs*Nfft):round(20*1000/Fs*Nfft)));
% energy8k_20k(j)=mean(power_f(round(8*1000/Fs*Nfft):round(20*1000/Fs*Nfft)));
% energy150k_4m(j)=mean(power_f(round(150*1000/Fs*Nfft):round(4000*1000/Fs*Nfft)));
% energy150k_500k(j)=mean(power_f(round(150*1000/Fs*Nfft):round(500*1000/Fs*Nfft)));
% 
% power_GL=10*log10(abs_Y8.^2/1e-12);  
% energy2k_20k_E(j)=mean(power_GL(round(2*1000/Fs*Nfft):round(20*1000/Fs*Nfft)));
% energy8k_20k_E(j)=mean(power_GL(round(8*1000/Fs*Nfft):round(20*1000/Fs*Nfft)));
% energy150k_4m_E(j)=mean(power_GL(round(150*1000/Fs*Nfft):round(4000*1000/Fs*Nfft)));
% energy150k_500k_E(j)=mean(power_GL(round(150*1000/Fs*Nfft):round(500*1000/Fs*Nfft)));



% energy2k_20k1(j)=10*log10(energy2k_20k(j));
% energy8k_20k1(j)=10*log10(energy8k_20k(j));
% energy150m_4m1(j)=10*log10(energy150m_4m(j));


end
% cell = {'��ѹ�ȼ�','ֱ������','2k-20k','8k-20k','0.15M-4M','0.15M-0.5M','2k-20kE','8k-20kE','0.15M-4ME','0.15M-0.5ME'};
% 
% VV = ones(1,j)*voltage_list;
% e=[VV;dc_component-dc_component0;energy2k_20k;energy8k_20k;energy150k_4m;energy150k_500k;energy2k_20k_E;energy8k_20k_E;energy150k_4m_E;energy150k_500k_E]';
% 
% filename_excel = '���������ݹ�����new.xlsx';  %�������½���һ�����
% xlRange = 'A1';
% xlswrite(filename_excel,cell,sheet,xlRange)
% xlRange = 'A2';
% xlswrite(filename_excel,e,sheet,xlRange)
% 
% 
% dc_component=dc_component';
% dc_component0=dc_component0';

