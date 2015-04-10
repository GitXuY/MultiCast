load gainChannel
load test_step1_new
% *************************************************************************
%                       参量设置
% *************************************************************************
%-------------------输入参数------------------------
% 多播组总数
TOTAL_MULTIGROUP=2;
% 用户总数
TOTAL_USER=100;
% 子载波总数
TOTAL_SUB=64;
% 子载波带宽，一般现在取0.3125MHz
BAND_SUB= 0.3125;
% 基站功耗对效用的权重
alpha_1=1;
% 用户功耗对效用的权重
alpha_2=1;
% 表示用户没接收单位数据业务的额外功耗
xi=0.001;
% 用户兴趣矩阵，用户对哪个多播组感兴趣
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);
% 功率约束 单位W
MAX_POWER_Pth=1;
% 最小发送速率约束 单位？
MIN_RATE_Rmin=1.5;

% omegaM;

%-----------------函数局部变量----------------------
%拉格朗日对偶乘子初始化
lambda=2;%对偶变量lambda
mu=ones(1,TOTAL_MULTIGROUP);%对偶变量mu

%迭代参数
isScheduleDone=0;% 资源分配完成标志
MAX_ITERATIONS=10000;%最大迭代次数
timesIteration=1;%初始化迭代计数器
% *************************************************************************
%                       阶段1：子载波分配
% *************************************************************************
% 每个子载波的功率相同，先把子载波分配给多播业务组

% 初始化子载波功率
powerSub_Pn=MAX_POWER_Pth/TOTAL_SUB.*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
% 初始化载波分配集合scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
% 每个多播组集合中最差接收SINR
minSNR_gamma=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
for iM=1:TOTAL_MULTIGROUP
    minSNR_gamma(:,iM)=gainChannel(criterionUserM(1,iM),:,LOOP); 
end
% 计算每个子载波在对应多播组使用时的速率
rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%每个子载波的传输速率
for iM=1:TOTAL_MULTIGROUP
    for iN=1:TOTAL_SUB
        rateSub(iN,iM)...
             =BAND_SUB*log2(1+powerSub_Pn(iN,iM)...
                  *minSNR_gamma(iN,iM));
    end
end




