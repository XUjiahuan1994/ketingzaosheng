function [ Yc ] = fs_same_th( fs_same_octave_info,y8,Fs,oc6)

%% %%%%%%%%%%%%%%%%%%%%
sf = Fs;      %采样频率
x = y8; %按行输入数据
%% 定义三分之一倍频程的中心频率
fc=fs_same_octave_info(1,:)';
%% 取中心频率总的长度
nc=length(fc);
%%oc6=5000; 
%% 取输入数据的长度                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
n=length(x);
%% 大于并最接近n的2的幂次方长度
nfft=2^nextpow2(n);
%% 进行FFT变换    
a=fft(x,nfft);      
for j=1:nc
    fl=fc(j)-oc6;   %下限频率
    fu=fc(j)+oc6;   %上限频率  
   %% 如果上限频率大于折叠频率则循环中断  我们的频率没有大于折叠频率的 
    if fu>sf/2
        m=j-1; 
        break;
    end
    freq(1,j)=fl;
    freq(2,j)=fu;    
    nl(j)=round(fl*nfft/sf+1);   %下限频率对应的序号
    nu(j)=round(fu*nfft/sf+1);   %上限频率对应的序号   
    %% 以每个中心频率段为通带进行带通频域滤波
    b=zeros(nfft,1);               %定义了一个为0的数组 此处也相当于将直流分量置零
    b(nl(j):nu(j))=a(nl(j):nu(j)); %将特定要处理的这个频段内的值置为实际的，其他的均为0
    b(nfft-nu(j)+1:nfft-nl(j)+1)=a(nfft-nu(j)+1:nfft-nl(j)+1);   %将对称部分的值置为实际值
    c=ifft(b,nfft);                         %转换到时域
    yc(j)=sqrt(var(real(c(1:n))));          %标准差的平方  再开根  得到的还是标准差啊
    Yc(j)=20*log10(yc(j)/2e-5);    
end
    nn=[nl;nu];
end

