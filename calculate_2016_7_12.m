%%��������
% dc_component(j)  ��ͬ��ѹ�ȼ�ֱ��������ֵ������j��Ӧ��ͬ���ļ�
% dc_component0(j) ��������ֱ��������ֵ
%���ڼ�������֮һ��Ƶ�̵�����Ƶ�ʺʹ����20Hz-20KHz
clear;
clc; 
close all;
%��ȡexcel�ļ�
[A_third_octave_info]=xlsread('20_150k������֮һ��Ƶ��');
%  [A_third_octave_info]=xlsread('20_150k�ĵȿ�Ƶ��5K');
center_freq=A_third_octave_info(1,:)';
NN=length(center_freq);

%% ================ ȷ���ļ�·�� ====================
% �޸��ļ�·�������뱳��������������ε�����������

% path0='E:\�ظ�ѹֱ�����\����������\0\';     %��������
% path='E:\�ظ�ѹֱ�����\����������\1100\'; 
% path0='C:\Users\m1881\Desktop\���ε���������������\����\';     %��������
% path0='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120712_20m+\0\';
% path0='D:\����\shuju\20120531+\0\';
 path0='D:\����\shuju\20120712_20m+\0\';

% path='C:\Users\m1881\Desktop\���ε���������������\����\';
% path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120518+\1000\'; 
% path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20130531xiawu\1000\';
% path='E:\�ظ�ѹֱ�����\ʵ���߶�����\20120712_20m+\1000\';
path='D:\����\shuju\20130531xiawu\1000\';
% path='D:\����\shuju\20120712_20m-\-700\';

size_n=1024*1024*2;                   % ��ȡ����
% �޸��ļ���
filenames=dir([path,'*.txt']);    %�г���ǰĿ¼�·���������ʽ���ļ��к��ļ�
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
    disp(['���������е�',num2str(j),'�����ݡ�����',num2str(length(filenames)-j),'������'])      %��command ��������ʾ�����

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
x(1)=data_7(1);x(2)=data_7(3);
clear data_voltage1;
%% ��ȡ��������
% filename=strcat(path0,'2014-03-08-21_10_11.txt');      % ��������62.5MHz�ĵ�����
%  filename=strcat(path0,'20120518084929.txt');        %20120518+  ��������62.5MHz��
%  filename=strcat(path0,'20120518085138.txt');        %20120518-  ��������62.5MHz��
% filename=strcat(path0,'2012-05-31-08_56_01.txt');     %20120531 ��������62.5MHz�� 531����û�б�������
filename=strcat(path0,'2012-07-12-12_25_56.txt');   %20120712_20m+��������62.5MHz��
% filename=strcat(path0,'2012-07-12-12_28_46.txt');   %20120712_20m-��������62.5MHz��
% filename=strcat(path0,'2012-07-12-12_25_56.txt');    %20120712+
%  filename=strcat(path0,'2012-07-12-12_26_38.txt');    %20120712-
% filename=strcat(path0,'2012-07-10-17_15_24.txt');    %20120710+
% filename=strcat(path0,'2012-07-10-17_13_37.txt');    %20120710-
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

%     if Fs8~=Fs0
        Fs=1024*1024*62.5;
%     else if Fs8==Fs0
%         Fs=Fs8;
%         end
%     end

    clear Fs0; clear Fs8;

%������ȡ
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

% figure(1); plot(t*1000,X0); grid on;
% xlabel('ʱ��(ms)');  ylabel('����(mA)');
% xlim([0,1000*T])
% title('������������ʱ����')
% figure(2); plot(t*1000,X8); grid on;
% xlabel('ʱ��(ms)'); ylabel('���ε���(mA)');
% xlim([0,1000*T])
% title('���ε�����������ʱ����')

%����Ҷ�任
Nfft=2^nextpow2(N);
nfft=1:1:Nfft/2+1;
Y8=fft(X8,Nfft); 
Y0=fft(X0,Nfft);


dc_component(j)=Y8(1)/N;
dc_component0(j)=Y0(1)/N;    %��һ�������ֱ������


Y8(1)=0; Y0(1)=0;    % ȥ��ֱ������
f=Fs/2*linspace(0,1,Nfft/2+1);
abs_Y8=2*abs(Y8(nfft))/N;
abs_Y0=2*abs(Y0(nfft))/N;

clear X0;

% figure(3); subplot(211);        % �����źŵĵ��߷�ֵƵ��
% plot(f(1:1300),(abs_Y0(1:1300))); grid on;
% xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
% title('�����źŵĵ��߷�ֵƵ��')
% xlim([10,13200]);
% 
% figure(3); subplot(212);        % ��ѹ�źŵĵ��߷�ֵƵ��
%% ��ͼ
% figure(3)
% plot(f(1:1048577),(abs_Y8(1:1048577))); grid on;
% xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
% title('��ѹ�źŵĵ��߷�ֵƵ��')
% xlim([10,1048577]);



%% ============= ��Ƶг������ =================
N50=ceil(1200*Nfft/Fs);     % ��ֹƵ��1200Hz
w=zeros(20,1);
i=1;
for n=2:1:N50-1
    if abs_Y8(n)>abs_Y8(n-1)&&abs_Y8(n)>abs_Y8(n+1)   %����n��ķ�ֵ�������˵ķ�ֵ
        w(i)=n-1;                                     %����
        i=i+1;
    end
end
w(w(:)==0)=[];  % ��Ϊ0�ĳ�ʼֵ���޳���
n_w=length(w);

clear abs_Y8; clear abs_Y0; clear i;

fw=w(1)*Fs/Nfft;      
ww=40/Fs*pi;                  %�Ƕ� ��Ϊ�������������50Hz�������Լ��㲻��������50Hz������������50Hz������
r=sqrt((1-sin(ww)/cos(ww)));  %��������1
wp=pi*fw/Fs;                  %�Ƕ�ֵ��
[B,A]=Notch03(wp,r);   % �����ݲ��˲���    rӦ���Ǽ���İ뾶   �õ���B A�� �˲�����ϵ��
[H]=freqz(B,A,Nfft);   %Ȼ��������˲�����Ƶ�����ԣ�

 clear B; clear A; clear fw; clear N50;

W=ones(Nfft,1);
for n=1:1:n_w
    for i=2:1:Nfft/10
        W(i+w(n)-w(1))=H(i);
        W(Nfft-i+2)=W(i);
    end
    Y8=Y8.*W;    %�����ݲ��˲���
end

clear W; clear H; clear w;

sum_08=sum(X8);
ave=sum_08/N;  
y=real(ifft(Y8,Nfft))+ave;

%% ���ƹ�Ƶг���������ź�
% figure(4); plot(t*1000,y); grid on;
% xlabel('ʱ��(ms)');  ylabel('����(mA)');
% axis([0,1000*T,-200,200])
% title('ȥ����Ƶг����ʱ����')
%% ���ݲ��Ժ�Ľ�����в�ֵ����
for i=1:Nfft
    ff(i)=(i-1)*Fs/(Nfft+1);
end            %��������Ƶ�ʵ��Ƶ��ֵ

%����ÿһ��ķ�ֵ
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
%% Ƶ���ͼ
% figure(5)
% plot(f(1:1300),abs(abs_Y8(1:1300)));
% hold on 
% plot(f(1:1300),abs(abs_Y8_1(1:1300)),'r');
% grid on;
% xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
% title('��ѹ�źŵĵ��߷�ֵƵ��')
% xlim([10,1300]);
% clear abs_Y8;
% clear y;
%% ============= ����غ������� =================
Z=(Y8.*conj(Y0))/Nfft;
Thn=ceil(200000*Nfft/Fs);
thn=Thn:1:Nfft/2+1;  %����200k��Ƶ�ײ�������

clear Y0;

% ��̬������ֵ
r=1;
ncorr=zeros(20,1);   %��¼��ȥ����Ƶ�ʵ�
% =============�������=================
for n=1:1:7                
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
        XS0=polyfit(x_round,Y8_round,2);  %  ������ȥ���
        for p=(n_corr(i)-n0):(n_corr(i)+n0)
            Y8(p)=1*polyval(XS0,f(p)); 
        end
        Y8(Nfft-(n_corr(i)+n0)+1:Nfft-(n_corr(i)-n0)+1)=fliplr(Y8(n_corr(i)-n0:n_corr(i)+n0));  
        
    end
    
end

clear Z; clear abs_Z; clear n0; clear Thn; 
clear n_corr; clear thn; clear i; clear r;

% ncorr*Fs/Nfft       %��Ϊȥ����խ������ģ��Ƶ��ֵ
y8=real(ifft(Y8,Nfft))+ave;
abs_Y8=2*abs(Y8(nfft))/Nfft;   %���߷�ֵ
clear Y8; clear ave; 
clear abs_Y8;

%% ================ С������ ======================
lev=7;
denoise_result=wden(y8,'rigrsure','s','one',lev,'db4');


figure(8); plot(t*1000,y8); grid on;
xlabel('ʱ��(ms)'); ylabel('���ε���(mA)');
xlim([0,1000*T])
title('���ε�����������ʱ����')

%% ================ ���ݷ��� ====================== 
%��ƽ��ֵ
denoise_result_mean=mean(denoise_result);
current_T=[current_T,denoise_result_mean];
% Y8=fft(y8,Nfft);
% Y8(1)=0;    % ȥ��ֱ������
% abs_Y8=2*abs(Y8(nfft))/N;            %1024*1024����
% abs_Y8(1)=1e-6; 


% [ Yc(:,j) ] = fs_same_th( A_third_octave_info,denoise_result,Fs,1000);
 [ Yc(:,j) ] = fs_same_th( A_third_octave_info,denoise_result,Fs,2500);   %�����õ�Ҳ��û�о���С��������ź�
%  [ Yc(:,j) ] = one_3th( A_third_octave_info,denoise_result,Fs);   %�����õ�Ҳ��û�о���С��������ź�
end

current_T=current_T';

figure (6)
plot(center_freq,Yc(:,2:j),'*-')
xlabel('����Ƶ�ʵ�')
title('700KV��������')
grid on
% 
% for i=1:NN
%     AVE(i)=mean(Yc(i,:));
% end
% figure (7)
% plot(center_freq,AVE,'*-')
% xlabel('����Ƶ�ʵ�')
% title('0KV��������ƽ��')
% grid on
% AVE=AVE';
% 
dc_component=dc_component';
