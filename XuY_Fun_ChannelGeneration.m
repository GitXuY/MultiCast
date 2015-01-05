function [gainChannel]=XuY_Fun_ChannelGeneration(~)
% XuY_Fun_ChannelGeneration: Generate Channel Model
%   -�����ŵ��������
%   -�û�����: �û��ڰ뾶ΪС���뾶��Բ�Ͼ��ȷֲ�
%   -�ŵ�����: �ŵ�����С�߶�ģ�ͣ������ŵ�������߶�ģ�ͣ����ɿռ�·�����ӰЧӦ��

% Author	: XuY
% Modified	: 2014-11-29 12:44

% Create and Modify Date and History :
% -2015/1/5 
%       �������ŵ�������ѭ�����棬��֤ÿ��loopʱ����С�߶�˥�䶼���������ɵ�

% Error Case :
% -

% Call 
% - function [lossRayleigh]=XuY_Fun_lossRayleigh(TOTAL_SUB,h)

% Called
% -

clc
clear all

% ����˵��
% ����С��Բ�Σ��뾶Ϊ500 m����վλ��Բ�ģ�
% �ನ�û�λ��һ��СԲ�ڣ�Բ�����վ��ľ�����[ 100m, 400m ]�ھ��ȷֲ�

%******************************************************
%                ��������
%******************************************************
N0= 10^(-14.4); %�����������ܶȣ���λW/MHz
RAD_CELL = 500;%С���뾶����λ:m
TOTAL_USER = 100;%����û���
LOOP = 10; %�ŵ�ʵ�ִ���
j = sqrt(-1);%������λ
TOTAL_SUB=64;%���ز�����
BAND_SUB= 0.3125; %���ز�����һ������ȡ0.3125MHz

lossLargeScare=zeros(1,TOTAL_USER,LOOP);%��߶��ŵ�˥������
lossSmallScare=zeros(1,TOTAL_SUB,LOOP);%С�߶��ŵ�˥������
gainChannel = zeros( TOTAL_USER, TOTAL_SUB, LOOP );%ÿ���û����վ����ŵ�����

for iLoop = 1:LOOP
    %******************************************************
    %                �û�����
    %******************************************************
    %�û��ڰ뾶ΪС���뾶��Բ�Ͼ��ȷֲ�
    locUser = ones(1,TOTAL_USER);%�û�λ�� 
    iUser=1;
    while iUser<=TOTAL_USER
        xUser = RAD_CELL*(2*rand-1);%�û����������ھ��ȷֲ�
        yUser = RAD_CELL*(2*rand-1);
        if (xUser^2+yUser^2) <= RAD_CELL^2 && (xUser^2+yUser^2)>100^2
            %�޳����������ڵ���Բ����û� && �޳����վ̫�����û�
             locUser(1,iUser) = xUser + j*yUser;
             iUser=iUser+1;
        end
    end

    %******************************************************
    %                ��߶��ŵ�����
    %******************************************************
    % �ŵ�����С�߶�ģ�ͣ������ŵ�������߶�ģ�ͣ����ɿռ�·�����ӰЧӦ��
    % lossLargeScare=zeros(1,TOTAL_USER,LOOP);%��߶��ŵ�˥������

    lossPath = 17.4 + 37.6.*log10( abs(locUser) );%·����ľ���
    lossShadow = 8.*randn(1,TOTAL_USER);%��ӰЧӦ����
    lossLargeScare(:,:,iLoop)=lossPath + lossShadow;
    lossLargeScare(:,:,iLoop) = 10.^(-lossLargeScare(:,:,iLoop)./10);%��������ԭ
    %******************************************************
    %                С�߶ȶ��ŵ�����
    %******************************************************
    % �������ŵ������ó������ɣ���Ϊ��ʱ����Ҫ�о�С�߶��ŵ��������
    %�����ŵ����,ת�ú�,����Ϊ1*TOTAL_SUB
    lossSmallScare(:,:,iLoop)=XuY_Fun_lossRayleigh(TOTAL_SUB,1)';
    %******************************************************
    %                ��С�ŵ��ϳ�
    %******************************************************
    % lossLargeScare=zeros(1,TOTAL_USER,LOOP);
    % lossSmallScare=zeros(1,TOTAL_SUB,LOOP);
    % gainChannel = zeros( TOTAL_USER, TOTAL_SUB, LOOP );
    % ǰ����ʽ����������Ƕ�������������ֱ�����
    gainChannel(:,:,iLoop)=(lossLargeScare(:,:,iGain)'*lossSmallScare(:,:,iGain))./(N0*BAND_SUB);
     
end % for iLoop = 1:LOOP

save gainChannel

end % Function

