function [lossRayleigh]=XuY_Fun_lossRayleigh(TOTAL_SUB,h)
% XuY_Fun_lossRayleigh: generate Rayleigh Channel
%   ��������˥������

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

% �����ŵ�˵����
% N�����ز�����
% h��ƽ���ŵ��������棬��dB��ʽ
% ����ģ�ͼ����ź�ͨ�������ŵ�֮�����źŷ����������(����˥�䡱),�������������ֲ���
% �����źŷ�Ϊ�����������֣���ÿһ���ֶ��Ƕ��·���źŵĵ��ӣ���·��������һ��������ʱ�����ǵĺ;������˹�ֲ���
% �����Ⱦ����������������͵Ŀ�ƽ���������������ֲ��ˡ�

%******************************************************
%                ��������
%******************************************************
BAND_SUB= 0.3125; %���ز�����һ������ȡ0.3125MHz
freqSample=TOTAL_SUB*BAND_SUB;%����Ƶ��= ���ز�����*���ز�����
timeSample=1/freqSample; %����ʱ��=1/����Ƶ��

%ÿ���Ĺ��ʷ��䣬Relative Path Power  db
powerRelativePath_DB=[0 -8.96 -17.37 -26.06 -34.74 -43.43];
powerRelativePath=10.^(powerRelativePath_DB/10);

%ÿһ���ķ��书��
powerPath=powerRelativePath/sum(powerRelativePath);
% Pmax*(powerRelativePath/sum(powerRelativePath))?

delayPath=[0 2380 5450 8370 13280 19270]/10;% ÿ����ʱ��ns
nSub=floor(delayPath/(timeSample*10.^3))+1;%����ÿһ����ʱ����������ز�

%******************************************************
%                �ŵ�����
%******************************************************
%������һ���ĸ���˹������̣������������γ������ֲ�
gaussianRandom_Re=randn(6,1);%ʵ��
gaussianRandom_Im=randn(6,1);%�鲿
gaussianRandom=(gaussianRandom_Re+1i*gaussianRandom_Im)/sqrt(2);%��һ��

%��˥���䵽��Ӧ�����ز�
tdSubAM=zeros(TOTAL_SUB,1);
for i=1:6
    tdSubAM(nSub(i),1)=(10^(h/20))*sqrt(powerPath(i))*gaussianRandom(i,1);
end

%����FFT�任����ʱ���ΪƵ����
fdSubAM=fftshift(fft(tdSubAM));

lossRayleigh=abs(fdSubAM);

end