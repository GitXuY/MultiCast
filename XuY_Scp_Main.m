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
% 后期完善编程时再考虑上述warning
clear all

%载入信道模型
load gainChannel
% *************************************************************************
%                           参量设置
% *************************************************************************
TOTAL_MULTIGROUP=2;%多播组总数
TOTAL_USER=100;%用户总数
TOTAL_SUB=64;%子载波总数
BAND_SUB= 0.3125; %子载波带宽，一般现在取0.3125MHz
MAX_POWER_Pth=1;%功率约束 单位W
MIN_RATE_Rmin=1.5;% 最小发送速率约束 单位？

alpha_1=1;%基站功耗对效用的权重
alpha_2=1;%用户功耗对效用的权重
xi=0.001;%表示用户没接收单位数据业务的额外功耗
betaM=zeros(1,TOTAL_MULTIGROUP);% 由多播推送集合中用户不愿意接收该业务造成的收益损失
% LOOP=10;%信道循环次

%用户兴趣矩阵，用户对哪个多播组感兴趣
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);
% *************************************************************************
%                         问题上限设计
% *************************************************************************
% 该方案的迭代环节分为两步:
% (1) 目标函数外部min迭代: XuY_Fun_UpperLimit_Min;
% 输入参数
%   场景参数
% XuY_Fun_UpperLimit_Min.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;%多播组总数
% XuY_Fun_UpperLimit_Min.TOTAL_USER=TOTAL_USER;%用户总数
% XuY_Fun_UpperLimit_Min.TOTAL_SUB=TOTAL_SUB;%子载波总数
% XuY_Fun_UpperLimit_Min.BAND_SUB=BAND_SUB;%子载波带宽，一般现在取0.3125MHz
% % 输出参数
Min_Output.lambda=100;
Min_Output.mu=0.1;

% (2) 目标函数内部max迭代: XuY_Fun_UpperLimit_Max;
% 输入参数
%   场景参数
Max_Input.TOTAL_MULTIGROUP  = TOTAL_MULTIGROUP;%多播组总数
Max_Input.TOTAL_USER        = TOTAL_USER;%用户总数
Max_Input.TOTAL_SUB         = TOTAL_SUB;%子载波总数
Max_Input.BAND_SUB          = BAND_SUB;%子载波带宽，一般现在取0.3125MHz
Max_Input.gainChannel       = gainChannel;
Max_Input.interestUM        = interestUM;
Max_Input.alpha_1           = alpha_1;
Max_Input.alpha_2           = alpha_2;
Max_Input.xi                = xi;
Max_Input.betaM             = betaM;
%   迭代参数
Max_Input.lamda             = Min_Output.lambda;
Max_Input.mu                = Min_Output.mu;

% [ Max_Output ] = XuY_Fun_UpperLimit_Max( Max_Input );
% *************************************************************************
%                           分步优化方法1
% *************************************************************************
% 分步优化方法1:
%       多播推送用户集合选择 + 子载波与功率联合分配
% 该方案分为三步：
%   (1) 初始化子载波分配与功率分配，计算“初始”多播推送用户集合选择
%   (2) 给定“初始”多播推送用户集合下，计算“初始”子载波与功率分配
%   (3) 多播推送用户集合选择 + 子载波与功率联合分配 联合迭代

% ====================================================
% (1) 计算“初始”多播推送用户集合选择
% ====================================================
% 初始化载波分配集合scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
for iScheduleSub=1:TOTAL_SUB
    %利用随机序列保证每行只有一个1，即每个子载波只能被分配到一个多播组
    index_rho=randperm(TOTAL_MULTIGROUP);
    scheduleSub_rho(iScheduleSub,index_rho(1))=1;
end

% 初始化子载波功率
powerSub_Pn=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
tempPn=rand(TOTAL_SUB,1);
%归一化后乘以总功率，保证每个子载波的功率之和等于总功率
oneTempPn=MAX_POWER_Pth*tempPn./sum(tempPn(:));
% 把每个子载波的功率对应到每个多播组中
for iPnN=1:TOTAL_SUB
    for iPnM=1:TOTAL_MULTIGROUP 
        if scheduleSub_rho(iPnN,iPnM)==1
            powerSub_Pn(iPnN,iPnM)=oneTempPn(iPnN,1);
        end
    end
end

% 输入参数
%   场景参数
OP1_Step1_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;%多播组总数
OP1_Step1_Input.TOTAL_USER=TOTAL_USER;%用户总数
OP1_Step1_Input.TOTAL_SUB=TOTAL_SUB;%子载波总数
OP1_Step1_Input.BAND_SUB=BAND_SUB;%子载波带宽，一般现在取0.3125MHz
OP1_Step1_Input.gainChannel=gainChannel;
OP1_Step1_Input.alpha_1=alpha_1;
OP1_Step1_Input.alpha_2=alpha_2;
OP1_Step1_Input.xi=xi;
OP1_Step1_Input.interestUM=interestUM;
%   迭代参数
OP1_Step1_Input.scheduleSub_rho=scheduleSub_rho;
OP1_Step1_Input.powerSub_Pn=powerSub_Pn;
%   约束参数
OP1_Step1_Input.MAX_POWER_Pth=MAX_POWER_Pth;%功率约束 单位W

% Step1
[ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input );
% ====================================================
% (2) 计算“初始”子载波与功率分配
% ====================================================
% 输入参数
%   场景参数
OP1_Step2_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;
OP1_Step2_Input.TOTAL_SUB=TOTAL_SUB;
OP1_Step2_Input.alpha_1=alpha_1;
OP1_Step2_Input.alpha_2=alpha_2;
OP1_Step2_Input.BAND_SUB=BAND_SUB;
%   迭代参数
OP1_Step2_Input.omegaM=OP1_Step1_Output.omegaM;
OP1_Step2_Input.minSNR_gamma=OP1_Step1_Output.minSNR_gamma;
%   约束参数
OP1_Step2_Input.MAX_POWER_Pth=MAX_POWER_Pth;
OP1_Step2_Input.MIN_RATE_Rmin=MIN_RATE_Rmin;

% Step2
[OP1_Step2_Output]=XuY_Fun_StepOptimize1_Step2(OP1_Step2_Input);

% ====================================================
% (3) 多播推送用户集合选择 + 子载波与功率联合分配 联合迭代
% ====================================================
% while 1
% --------------多播推送用户集合选择---------------------
%   更新迭代参数
OP1_Step1_Input.scheduleSub_rho=OP1_Step2_Output.scheduleSub_rho;
OP1_Step1_Input.powerSub_Pn=OP1_Step2_Output.powerSub_Pn;
[ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input );

% --------------子载波与功率联合分配---------------------
%   更新迭代参数
OP1_Step2_Input.omegaM=OP1_Step1_Output.omegaM;
OP1_Step2_Input.minSNR_gamma=OP1_Step1_Output.minSNR_gamma;
[OP1_Step2_Output]=XuY_Fun_StepOptimize1_Step2(OP1_Step2_Input);

% end % while 1

% *************************************************************************
%                         分步优化方法2
% *************************************************************************
% 多播推送用户集合选择与子载波分配 + 功率分配
% 该方案的迭代环节分为两步:
% (1) 给定多播推送用户集合与子载波分配下，功率分配算法
% XuY_Fun_StepOptimize2_Step1;
% (2) 给定功率分配下，多播推送用户集合选择与子载波分配
% XuY_Fun_StepOptimize2_Step2;
% *************************************************************************
%                         结果输出&作图
% *************************************************************************