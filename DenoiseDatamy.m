%% ���ε������ݽ��봦��

% ���㷨�ʺ϶��ظ�ѹֱ�����ε�������Ĵ���
% ���б��㷨����Ҫͬһ��ı����������ѹ����ͬʱ����

% ----------------------- ������� -----------------------------
% Norch03--��Ƶ�ݲ��˲���
%-------------------------------------------------------------------

clear; clc; 
close all;
%% ================ ȷ���ļ�·�� ====================
% �޸��ļ�·�������뱳��������������ε�����������
path0='D:\corona\2014\0\';     %��������
path='D:\corona\2014\800\';     %���ε���
size=1024*256*4+7;      % ��ȡ����

% �޸��ļ���
filename=strcat(path,'2014-03-08-23_00_53.txt');  %���ε���
[fid]=fopen(filename,'r+');
data_voltage=fread(fid,size,'schar'); fclose(fid);
filename=strcat(path0,'2014-03-08-21_10_11.txt'); %��������
[fid]=fopen(filename,'r+');
data_background=fread(fid,size,'schar'); fclose(fid);


%% ============== ��ȡ�����ʡ��洢��� ===================

switch data_voltage(1)
    case 0
        Fs8=1000;
	case 1
        Fs8=1000*2;
	case 2
        Fs8=1000*4;
	case 3
        Fs8=1000*8;
	case 4
        Fs8=1000*16;
	case 5
        Fs8=1000*32;
	case 6
        Fs8=1000*64;
	case 7
        Fs8=1000*128;
	case 8
        Fs8=1000*256;
	case 9
        Fs8=1000*512;
    case 10
        Fs8=1000*1000;
    case 18
        Fs8=1000*1000*62.5;
    case 21
        Fs8=1000*1000*500;
end

switch data_background(1)
    case 0
        Fs0=1000;
	case 1
        Fs0=1000*2;
	case 2
        Fs0=1000*4;
	case 3
        Fs0=1000*8;
	case 4
        Fs0=1000*16;
	case 5
        Fs0=1000*32;
	case 6
        Fs0=1000*64;
	case 7
        Fs0=1000*128;
	case 8
        Fs0=1000*256;
	case 9
        Fs0=1000*512;
    case 10
        Fs0=1000*1000;
    case 18
        Fs0=1000*1000*62.5;
    case 21
        Fs0=1000*1000*500;
end

if Fs8~=Fs0
    error('���󣡲����ʲ�ͬ��')
else if Fs8==Fs0
    Fs=Fs8;
    end
end

clear Fs0; clear Fs8;

%������ȡ
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

switch data_background(3)
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
    error('�������ݳ��Ȳ�һ�£�')
end

%% =============== �������� ==================
%����ѹ����ת��Ϊ����
N=size-7;
T=(N-1)/Fs;
t=0:1/Fs:T;

X8=10*voltage8*data_voltage(8:size)/127; 
X0=10*voltage0*data_background(8:size)/127; 

clear data_voltage; clear data_background; 
clear voltage8; clear voltage0;

figure(1); plot(t*1000,X0); grid on;
xlabel('ʱ��(ms)');  ylabel('����(mA)');
xlim([0,1000*T])
title('������������ʱ����')
figure(2); plot(t*1000,X8); grid on;
xlabel('ʱ��(ms)'); ylabel('����(mA)');
xlim([0,1000*T])
title('���ε�����������ʱ����')

%����Ҷ�任
Nfft=2^nextpow2(N);
nfft=1:1:Nfft/2+1;
Y8=fft(X8,Nfft); Y0=fft(X0,Nfft);
Y8(1)=0; Y0(1)=0;    % ȥ��ֱ������
f=Fs/2*linspace(0,1,Nfft/2+1);
abs_Y8=2*abs(Y8(nfft))/N;
abs_Y0=2*abs(Y0(nfft))/N;

% clear X0;

figure(3); subplot(211);        % �����źŵĵ��߷�ֵƵ��
plot(f,abs_Y0); grid on;
xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
title('�����źŵĵ��߷�ֵƵ��')
figure(3); subplot(212);        % ��ѹ�źŵĵ��߷�ֵƵ��
plot(f,abs_Y8); grid on;
xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
title('��ѹ�źŵĵ��߷�ֵƵ��')


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

clear W; clear H; clear w;

sum_08=sum(X8);
ave=sum_08/N;  
y=real(ifft(Y8,Nfft))+ave;

figure(4); plot(t*1000,y); grid on;
xlabel('ʱ��(ms)');  ylabel('����(mA)');
axis([0,1000*T,-150,150])
title('ȥ����Ƶг����ʱ����')

clear y;


%% ============= ����غ������� =================
Z=(Y8.*conj(Y0))/Nfft;
Thn=ceil(200000*Nfft/Fs);
thn=Thn:1:Nfft/2+1;  %����200k��Ƶ�ײ�������

clear Y0;

% ��̬������ֵ
r=1;
ncorr=zeros(20,1);   %��¼��ȥ����Ƶ�ʵ�
for n=1:1:7                % !!!!!!!!!!!!!!!�������!!!!!!!!!!!!!!!!!
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

ncorr*Fs/Nfft       %��Ϊȥ����խ������ģ��Ƶ��ֵ
y8=real(ifft(Y8,Nfft))+ave;
abs_Y8=2*abs(Y8(nfft))/Nfft;

clear Y8; clear ave; clear nfft;

figure(5); plot(f,abs_Y8); grid on;
xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
title('ȥ��խ�����ź�ĵ��߷�ֵƵ��')
xlim([0,20000]);
figure(6); plot(t*1000,y8);grid on;
xlabel('ʱ��(ms)');  ylabel('����(mA)');
axis([0,1000*T,-150,150])
title('ȥ��խ�����ź�ʱ����')

clear f; clear abs_Y8;


%% ================ С������ ======================
lev=7;
denoise_result=wden(y8,'rigrsure','s','one',lev,'db4'); 

figure(7); plot(t*1000,denoise_result);grid on;
xlabel('ʱ��(ms)');  ylabel('����(mA)');
axis([0,1000*T,-150,150])
title('С�������ʱ����')

%figure(8); plot(f,abs_Y8); grid on;
%xlabel('Ƶ��(Hz)'); ylabel('��ֵ|Y(f)|');
%title('С�������ĵ��߷�ֵƵ��')



per=sum(denoise_result.^2)/sum(X8.^2)*100 %������
err=norm(X8-denoise_result)/sqrt(N)  % �������








