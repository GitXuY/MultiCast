% Author	: XuY
% Modified	: 2014-12-1 13:44

% Create and Modify Date and History :
% -

% Error Case :
% -

% Call 
% -function [ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input )
% -function [ OP1_Step2_Output ] = XuY_Fun_StepOptimize1_Step2( OP1_Step2_Input )

clc
clear all
load gainChannel%�����ŵ�ģ��
% *************************************************************************
%                           ��������
% *************************************************************************
TOTAL_MULTIGROUP=2;%�ಥ������
TOTAL_USER=100;%�û�����
TOTAL_SUB=64;%���ز�����
BAND_SUB= 0.3125; %���ز�����һ������ȡ0.3125MHz
MAX_POWER_Pth=1;%����Լ�� ��λW
MIN_RATE_Rmin=1.5;% ��С��������Լ�� ��λ��

t=0;%��ʽ�滮�Ĳ�����������Ϊ0��������
% LOOP=10;%�ŵ�ѭ����

%�û���Ȥ�����û����ĸ��ಥ�����Ȥ
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);

% *************************************************************************
%                           �ֲ��Ż�����1
% *************************************************************************
% �ֲ��Ż�����1:
%       �ಥ�����û�����ѡ�� + ���ز��빦�����Ϸ���
% �÷����ĵ������ڷ�Ϊ������
%   (1) �������ز������빦�ʷ����£��ಥ�����û�����ѡ��
%   (2) �����ಥ�����û������£����ز��빦�����Ϸ����㷨

% ====================================================
% (1) �������ز������빦�ʷ����£��ಥ�����û�����ѡ��
% ====================================================
% �������
%   ��������
OP1_Step1_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;%�ಥ������
OP1_Step1_Input.TOTAL_USER=TOTAL_USER;%�û�����
OP1_Step1_Input.TOTAL_SUB=TOTAL_SUB;%���ز�����
OP1_Step1_Input.BAND_SUB=BAND_SUB;%���ز�����һ������ȡ0.3125MHz
OP1_Step1_Input.t=t;%��ʽ�滮�Ĳ�����������Ϊ0��������
OP1_Step1_Input.gainChannel=gainChannel;

%   Լ������
OP1_Step1_Input.MAX_POWER_Pth=MAX_POWER_Pth;%����Լ�� ��λW

% Step2
[ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input );
% ====================================================
% (2) �����ಥ�����û������£����ز��빦�����Ϸ����㷨
% ====================================================
% �������
%   ��������
OP1_Step2_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;
OP1_Step2_Input.TOTAL_SUB=TOTAL_SUB;
OP1_Step2_Input.omegaM=OP1_Step1_Output.omegaM;
OP1_Step2_Input.t=t;
OP1_Step2_Input.BAND_SUB=BAND_SUB;
OP1_Step2_Input.minSNR_gamma=OP1_Step1_Output.minSNR_gamma;

%   Լ������
OP1_Step2_Input.MAX_POWER_Pth=MAX_POWER_Pth;
OP1_Step2_Input.MIN_RATE_Rmin=MIN_RATE_Rmin;

% Step2
[OP1_Step2_Output]=XuY_Fun_StepOptimize1_Step2(OP1_Step2_Input);

% *************************************************************************
%                         �ֲ��Ż�����2
% *************************************************************************
% �ಥ�����û�����ѡ�������ز����� + ���ʷ���
% �÷����ĵ������ڷ�Ϊ����:
% (1) �����ಥ�����û����������ز������£����ʷ����㷨
% XuY_Fun_StepOptimize2_Step1;
% (2) �������ʷ����£��ಥ�����û�����ѡ�������ز�����
% XuY_Fun_StepOptimize2_Step2;
% *************************************************************************
%                         ������&��ͼ
% *************************************************************************