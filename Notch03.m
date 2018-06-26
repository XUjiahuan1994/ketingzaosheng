%Design notch filter based on all-pass filter:       设计一个全通滤波器          
%Note: the zeros and poles have the diff direction.  零点和极点在具有不同的方向 分子与分母的多项式系数呈倒序关系
%H(ej0) and H(ejpi) have the same amplitude.         有相同的幅值
%w0: The notch freqency(normally between 0 and pi).  陷波频率，0-pi之间
%r:  The radian of the conjunction poles               
%[H W] = Notch1(w0, r, N) returns the frequency response at frequencies designated in vector W, in radians/sample (normally between 0 and pi).
%返回值是在向量W中指定频率的频率响应
function [B, A] = Notch03(w0, r)
%Calculate pole freqency
tp = (1 + r*r) / (2*r);
w1 = acos(tp * cos(w0));   %arccos 值  反三角函数值
%Construct notch filter    %构建陷波滤波器
B1 = 1; B2 = -2*cos(w0);   B3 = 1;
A1 = 1; A2 = -2*r*cos(w1); A3 = r*r;
B = [B1 B2 B3]; A = [A1 A2 A3];
% Calculate and normal coefficient    计算系数
num = 2 - 2*cos(w0);
den = 1 - 2*r*cos(w1) + r*r;
coeff = num / den;
B = B ./ coeff;
