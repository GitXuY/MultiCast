function [lossRayleigh]=XuY_Fun_lossRayleigh(TOTAL_SUB,h)
% XuY_Fun_lossRayleigh: generate Rayleigh Channel
%   产生瑞利衰减矩阵

% Author	: XuY
% Modified	: 2014-11-29 22:44

% Create and Modify Date and History :
% -

% Error Case :
% -  

% Call 
% -

% Called
% -function [gainChannel]=XuY_Fun_ChannelGeneration(~)

% 瑞利信道说明：
% N是子载波个数
% h是平均信道功率增益，是dB形式
% 瑞利模型假设信号通过无线信道之后，其信号幅度是随机的(即“衰落”),其包络服从瑞利分布。
% 瑞利信号分为正交的两部分，而每一部分都是多个路径信号的叠加，当路径数大于一定数量的时候，他们的和就满足高斯分布。
% 而幅度就是两个正交变量和的开平方，就满足瑞利分布了。

%******************************************************
%                参量设置
%******************************************************
BAND_SUB= 0.3125; %子载波带宽，一般现在取0.3125MHz
freqSample=TOTAL_SUB*BAND_SUB;%采样频率= 子载波个数*子载波带宽
timeSample=1/freqSample; %抽样时间=1/采样频率

%每径的功率分配，Relative Path Power  db
powerRelativePath_DB=[0 -8.96 -17.37 -26.06 -34.74 -43.43];
powerRelativePath=10.^(powerRelativePath_DB/10);

%每一径的发射功率
powerPath=powerRelativePath/sum(powerRelativePath);
% Pmax*(powerRelativePath/sum(powerRelativePath))?

delayPath=[0 2380 5450 8370 13280 19270]/10;% 每径的时延ns
nSub=floor(delayPath/(timeSample*10.^3))+1;%计算每一径的时延所落的子载波

%******************************************************
%                信道生成
%******************************************************
%创建归一化的复高斯随机过程，方便后面叠加形成瑞利分布
gaussianRandom_Re=randn(6,1);%实部
gaussianRandom_Im=randn(6,1);%虚部
gaussianRandom=(gaussianRandom_Re+1i*gaussianRandom_Im)/sqrt(2);%归一化

%将衰落落到相应的子载波
tdSubAM=zeros(TOTAL_SUB,1);
for i=1:6
    tdSubAM(nSub(i),1)=(10^(h/20))*sqrt(powerPath(i))*gaussianRandom(i,1);
end

%进行FFT变换，从时域变为频率域
fdSubAM=fftshift(fft(tdSubAM));

lossRayleigh=abs(fdSubAM);

end