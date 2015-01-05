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
load gainChannel%载入信道模型
% *************************************************************************
%                           参量设置
% *************************************************************************
TOTAL_MULTIGROUP=2;%多播组总数
TOTAL_USER=100;%用户总数
TOTAL_SUB=64;%子载波总数
BAND_SUB= 0.3125; %子载波带宽，一般现在取0.3125MHz
MAX_POWER_Pth=1;%功率约束 单位W
MIN_RATE_Rmin=1.5;% 最小发送速率约束 单位？

t=0;%分式规划的参量，现在置为0，测试用
% LOOP=10;%信道循环次

%用户兴趣矩阵，用户对哪个多播组感兴趣
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);

% *************************************************************************
%                           分步优化方法1
% *************************************************************************
% 分步优化方法1:
%       多播推送用户集合选择 + 子载波与功率联合分配
% 该方案的迭代环节分为两步：
%   (1) 给定子载波分配与功率分配下，多播推送用户集合选择
%   (2) 给定多播推送用户集合下，子载波与功率联合分配算法

% ====================================================
% (1) 给定子载波分配与功率分配下，多播推送用户集合选择
% ====================================================
% 输入参数
%   场景参数
OP1_Step1_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;%多播组总数
OP1_Step1_Input.TOTAL_USER=TOTAL_USER;%用户总数
OP1_Step1_Input.TOTAL_SUB=TOTAL_SUB;%子载波总数
OP1_Step1_Input.BAND_SUB=BAND_SUB;%子载波带宽，一般现在取0.3125MHz
OP1_Step1_Input.t=t;%分式规划的参量，现在置为0，测试用
OP1_Step1_Input.gainChannel=gainChannel;

%   约束参数
OP1_Step1_Input.MAX_POWER_Pth=MAX_POWER_Pth;%功率约束 单位W

% Step2
[ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input );
% ====================================================
% (2) 给定多播推送用户集合下，子载波与功率联合分配算法
% ====================================================
% 输入参数
%   场景参数
OP1_Step2_Input.TOTAL_MULTIGROUP=TOTAL_MULTIGROUP;
OP1_Step2_Input.TOTAL_SUB=TOTAL_SUB;
OP1_Step2_Input.omegaM=OP1_Step1_Output.omegaM;
OP1_Step2_Input.t=t;
OP1_Step2_Input.BAND_SUB=BAND_SUB;
OP1_Step2_Input.minSNR_gamma=OP1_Step1_Output.minSNR_gamma;

%   约束参数
OP1_Step2_Input.MAX_POWER_Pth=MAX_POWER_Pth;
OP1_Step2_Input.MIN_RATE_Rmin=MIN_RATE_Rmin;

% Step2
[OP1_Step2_Output]=XuY_Fun_StepOptimize1_Step2(OP1_Step2_Input);

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