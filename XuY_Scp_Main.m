% Author	: XuY
% Modified	: 2014-12-1 13:44

% Create and Modify Date and History :
% - 2015/1/15 
%     - delete input parameter "OP1_Step1_Input.t" 
%     - add input parameter "OP1_Step1_Input.alpha_1"
%     - add input parameter "OP1_Step1_Input.alpha_2"

% Error Case :
% -

% Call 
% -function [ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input )
% -function [ OP1_Step2_Output ] = XuY_Fun_StepOptimize1_Step2( OP1_Step2_Input )

clc

% For improved performance, consider not using clear all within a script
% �������Ʊ��ʱ�ٿ�������warning
clear all

%�����ŵ�ģ��
load gainChannel
% *************************************************************************
%                           ��������
% *************************************************************************
TOTAL_MULTIGROUP=2;%�ಥ������
TOTAL_USER=100;%�û�����
TOTAL_SUB=64;%���ز�����
BAND_SUB= 0.3125; %���ز�����һ������ȡ0.3125MHz
MAX_POWER_Pth=1;%����Լ�� ��λW
MIN_RATE_Rmin=1.5;% ��С��������Լ�� ��λ��

alpha_1=1;%��վ���Ķ�Ч�õ�Ȩ��
alpha_2=1;%�û����Ķ�Ч�õ�Ȩ��
xi=0.001;%��ʾ�û�û���յ�λ����ҵ��Ķ��⹦��
betaM=zeros(1,TOTAL_MULTIGROUP);% �ɶಥ���ͼ������û���Ը����ո�ҵ����ɵ�������ʧ
% LOOP=10;%�ŵ�ѭ����

%�û���Ȥ�����û����ĸ��ಥ�����Ȥ
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);
% *************************************************************************
%                         �����������
% *************************************************************************
% �÷����ĵ������ڷ�Ϊ����:
% (1) Ŀ�꺯���ⲿmin����: XuY_Fun_UpperLimit_Min;
% �������
%   ��������
% XuY_Fun_UpperLimit_Min.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;%�ಥ������
% XuY_Fun_UpperLimit_Min.TOTAL_USER=TOTAL_USER;%�û�����
% XuY_Fun_UpperLimit_Min.TOTAL_SUB=TOTAL_SUB;%���ز�����
% XuY_Fun_UpperLimit_Min.BAND_SUB=BAND_SUB;%���ز�����һ������ȡ0.3125MHz
% % �������
Min_Output.lambda=100;
Min_Output.mu=0.1;

% (2) Ŀ�꺯���ڲ�max����: XuY_Fun_UpperLimit_Max;
% �������
%   ��������
Max_Input.TOTAL_MULTIGROUP  = TOTAL_MULTIGROUP;%�ಥ������
Max_Input.TOTAL_USER        = TOTAL_USER;%�û�����
Max_Input.TOTAL_SUB         = TOTAL_SUB;%���ز�����
Max_Input.BAND_SUB          = BAND_SUB;%���ز�����һ������ȡ0.3125MHz
Max_Input.gainChannel       = gainChannel;
Max_Input.interestUM        = interestUM;
Max_Input.alpha_1           = alpha_1;
Max_Input.alpha_2           = alpha_2;
Max_Input.xi                = xi;
Max_Input.betaM             = betaM;
%   ��������
Max_Input.lamda             = Min_Output.lambda;
Max_Input.mu                = Min_Output.mu;

% [ Max_Output ] = XuY_Fun_UpperLimit_Max( Max_Input );
% *************************************************************************
%                           �ֲ��Ż�����1
% *************************************************************************
% �ֲ��Ż�����1:
%       �ಥ�����û�����ѡ�� + ���ز��빦�����Ϸ���
% �÷�����Ϊ������
%   (1) ��ʼ�����ز������빦�ʷ��䣬���㡰��ʼ���ಥ�����û�����ѡ��
%   (2) ��������ʼ���ಥ�����û������£����㡰��ʼ�����ز��빦�ʷ���
%   (3) �ಥ�����û�����ѡ�� + ���ز��빦�����Ϸ��� ���ϵ���

% ====================================================
% (1) ���㡰��ʼ���ಥ�����û�����ѡ��
% ====================================================
% ��ʼ���ز����伯��scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
for iScheduleSub=1:TOTAL_SUB
    %����������б�֤ÿ��ֻ��һ��1����ÿ�����ز�ֻ�ܱ����䵽һ���ಥ��
    index_rho=randperm(TOTAL_MULTIGROUP);
    scheduleSub_rho(iScheduleSub,index_rho(1))=1;
end

% ��ʼ�����ز�����
powerSub_Pn=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
tempPn=rand(TOTAL_SUB,1);
%��һ��������ܹ��ʣ���֤ÿ�����ز��Ĺ���֮�͵����ܹ���
oneTempPn=MAX_POWER_Pth*tempPn./sum(tempPn(:));
% ��ÿ�����ز��Ĺ��ʶ�Ӧ��ÿ���ಥ����
for iPnN=1:TOTAL_SUB
    for iPnM=1:TOTAL_MULTIGROUP 
        if scheduleSub_rho(iPnN,iPnM)==1
            powerSub_Pn(iPnN,iPnM)=oneTempPn(iPnN,1);
        end
    end
end

% �������
%   ��������
OP1_Step1_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;%�ಥ������
OP1_Step1_Input.TOTAL_USER=TOTAL_USER;%�û�����
OP1_Step1_Input.TOTAL_SUB=TOTAL_SUB;%���ز�����
OP1_Step1_Input.BAND_SUB=BAND_SUB;%���ز�����һ������ȡ0.3125MHz
OP1_Step1_Input.gainChannel=gainChannel;
OP1_Step1_Input.alpha_1=alpha_1;
OP1_Step1_Input.alpha_2=alpha_2;
OP1_Step1_Input.xi=xi;
OP1_Step1_Input.interestUM=interestUM;
%   ��������
OP1_Step1_Input.scheduleSub_rho=scheduleSub_rho;
OP1_Step1_Input.powerSub_Pn=powerSub_Pn;
%   Լ������
OP1_Step1_Input.MAX_POWER_Pth=MAX_POWER_Pth;%����Լ�� ��λW

% Step1
[ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input );
% ====================================================
% (2) ���㡰��ʼ�����ز��빦�ʷ���
% ====================================================
% �������
%   ��������
OP1_Step2_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;
OP1_Step2_Input.TOTAL_SUB=TOTAL_SUB;
OP1_Step2_Input.alpha_1=alpha_1;
OP1_Step2_Input.alpha_2=alpha_2;
OP1_Step2_Input.BAND_SUB=BAND_SUB;
%   ��������
OP1_Step2_Input.omegaM=OP1_Step1_Output.omegaM;
OP1_Step2_Input.minSNR_gamma=OP1_Step1_Output.minSNR_gamma;
%   Լ������
OP1_Step2_Input.MAX_POWER_Pth=MAX_POWER_Pth;
OP1_Step2_Input.MIN_RATE_Rmin=MIN_RATE_Rmin;

% Step2
[OP1_Step2_Output]=XuY_Fun_StepOptimize1_Step2(OP1_Step2_Input);

% ====================================================
% (3) �ಥ�����û�����ѡ�� + ���ز��빦�����Ϸ��� ���ϵ���
% ====================================================
% while 1
% --------------�ಥ�����û�����ѡ��---------------------
%   ���µ�������
OP1_Step1_Input.scheduleSub_rho=OP1_Step2_Output.scheduleSub_rho;
OP1_Step1_Input.powerSub_Pn=OP1_Step2_Output.powerSub_Pn;
[ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input );

% --------------���ز��빦�����Ϸ���---------------------
%   ���µ�������
OP1_Step2_Input.omegaM=OP1_Step1_Output.omegaM;
OP1_Step2_Input.minSNR_gamma=OP1_Step1_Output.minSNR_gamma;
[OP1_Step2_Output]=XuY_Fun_StepOptimize1_Step2(OP1_Step2_Input);

% end % while 1

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