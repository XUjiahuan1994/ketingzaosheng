%Design notch filter based on all-pass filter:       ���һ��ȫͨ�˲���          
%Note: the zeros and poles have the diff direction.  ���ͼ����ھ��в�ͬ�ķ��� �������ĸ�Ķ���ʽϵ���ʵ����ϵ
%H(ej0) and H(ejpi) have the same amplitude.         ����ͬ�ķ�ֵ
%w0: The notch freqency(normally between 0 and pi).  �ݲ�Ƶ�ʣ�0-pi֮��
%r:  The radian of the conjunction poles               
%[H W] = Notch1(w0, r, N) returns the frequency response at frequencies designated in vector W, in radians/sample (normally between 0 and pi).
%����ֵ��������W��ָ��Ƶ�ʵ�Ƶ����Ӧ
function [B, A] = Notch03(w0, r)
%Calculate pole freqency
tp = (1 + r*r) / (2*r);
w1 = acos(tp * cos(w0));   %arccos ֵ  �����Ǻ���ֵ
%Construct notch filter    %�����ݲ��˲���
B1 = 1; B2 = -2*cos(w0);   B3 = 1;
A1 = 1; A2 = -2*r*cos(w1); A3 = r*r;
B = [B1 B2 B3]; A = [A1 A2 A3];
% Calculate and normal coefficient    ����ϵ��
num = 2 - 2*cos(w0);
den = 1 - 2*r*cos(w1) + r*r;
coeff = num / den;
B = B ./ coeff;
