function [ Yc ] = fs_same_th( fs_same_octave_info,y8,Fs,oc6)

%% %%%%%%%%%%%%%%%%%%%%
sf = Fs;      %����Ƶ��
x = y8; %������������
%% ��������֮һ��Ƶ�̵�����Ƶ��
fc=fs_same_octave_info(1,:)';
%% ȡ����Ƶ���ܵĳ���
nc=length(fc);
%%oc6=5000; 
%% ȡ�������ݵĳ���                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
n=length(x);
%% ���ڲ���ӽ�n��2���ݴη�����
nfft=2^nextpow2(n);
%% ����FFT�任    
a=fft(x,nfft);      
for j=1:nc
    fl=fc(j)-oc6;   %����Ƶ��
    fu=fc(j)+oc6;   %����Ƶ��  
   %% �������Ƶ�ʴ����۵�Ƶ����ѭ���ж�  ���ǵ�Ƶ��û�д����۵�Ƶ�ʵ� 
    if fu>sf/2
        m=j-1; 
        break;
    end
    freq(1,j)=fl;
    freq(2,j)=fu;    
    nl(j)=round(fl*nfft/sf+1);   %����Ƶ�ʶ�Ӧ�����
    nu(j)=round(fu*nfft/sf+1);   %����Ƶ�ʶ�Ӧ�����   
    %% ��ÿ������Ƶ�ʶ�Ϊͨ�����д�ͨƵ���˲�
    b=zeros(nfft,1);               %������һ��Ϊ0������ �˴�Ҳ�൱�ڽ�ֱ����������
    b(nl(j):nu(j))=a(nl(j):nu(j)); %���ض�Ҫ��������Ƶ���ڵ�ֵ��Ϊʵ�ʵģ������ľ�Ϊ0
    b(nfft-nu(j)+1:nfft-nl(j)+1)=a(nfft-nu(j)+1:nfft-nl(j)+1);   %���ԳƲ��ֵ�ֵ��Ϊʵ��ֵ
    c=ifft(b,nfft);                         %ת����ʱ��
    yc(j)=sqrt(var(real(c(1:n))));          %��׼���ƽ��  �ٿ���  �õ��Ļ��Ǳ�׼�
    Yc(j)=20*log10(yc(j)/2e-5);    
end
    nn=[nl;nu];
end

