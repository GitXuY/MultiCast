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
MIN_RATE_Rmin=20;
%47的时候PTS能达到1，RTS接近1，且3000次内能收敛

% omegaM;

%-----------------函数局部变量----------------------
%拉格朗日对偶乘子初始化
lambda=20;%对偶变量lambda
mu=ones(1,TOTAL_MULTIGROUP);%对偶变量mu
%-----------------初始化----------------------
% 初始化子载波功率
powerSub_Pn=MAX_POWER_Pth/TOTAL_SUB.*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
% 初始化载波分配集合scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);

%迭代参数
isScheduleDone=0;% 资源分配完成标志
MAX_ITERATIONS=3000;%最大迭代次数
timesIteration=1;%初始化迭代计数器
while isScheduleDone==0
% *************************************************************************
%                       阶段1：子载波分配
% *************************************************************************
% 每个子载波的功率相同，先把子载波分配给多播业务组
%--------------------------------------------------
%       计算每个子载波在对应多播组使用时的速率
%--------------------------------------------------
rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%每个子载波的传输速率
for iM=1:TOTAL_MULTIGROUP
    for iN=1:TOTAL_SUB
        rateSub(iN,iM)...
             =BAND_SUB*log2(1+powerSub_Pn(iN,iM)...
                  *minSNR_gamma(iN,iM));
    end
end
rateSub_greedy=rateSub;
rateM = sum(rateSub);
%--------------------------------------------------
%       贪婪算法分配子载波
%--------------------------------------------------
% 找到速率最大的子载波并定位
while sum(scheduleSub_rho(:))~=TOTAL_SUB
    [maxRate, bestSub]=max(rateSub_greedy);
    % 如果子载波在多个多播组中速率均为最大，随机选取一个多播组
    if length(bestSub) > 1
        bestMultigroup = ceil(length(bestSub)*rand(1));
    end
    scheduleSub_rho(bestSub(bestMultigroup),bestMultigroup)=1;
    rateSub_greedy(bestSub(bestMultigroup),:)=0;
end
% *************************************************************************
%                       阶段2：功率分配
% *************************************************************************
for iM=1:TOTAL_MULTIGROUP
    for iN=1:TOTAL_SUB
        if scheduleSub_rho(iN,iM)==1
            powerSub_Pn(iN,iM)=(mu(1,iM)+omegaM(1,iM))...
                *BAND_SUB./((alpha_1+lambda)*log(2))-1./minSNR_gamma(iN,iM);
            % 功率取正数
            powerSub_Pn(1,iM)=max(powerSub_Pn(1,iM),0);
        else
          powerSub_Pn(iN,iM)=0;
        end
    end
end
% *************************************************************************
%                   阶段3：对偶变量子梯度法迭代
% *************************************************************************
    % ====================================================
    %     （1）第一次迭代时使用初始化值来启动迭代过程
    % ====================================================
%     if timesIteration==1
        %Lambda的更新步长初始化
        stepLambda_c1= 50*(1e-1)/sqrt(timesIteration);
        
        %Mu的更新步长初始化
        stepMU_c2 = 1*(1e-3)/sqrt(timesIteration);
%     end
    % ====================================================
    %      （2）lambda和mu的更新
    % ====================================================
    %更新lambda
    lambda2 = lambda + stepLambda_c1*( sum(powerSub_Pn(:))-MAX_POWER_Pth);
    lambda2 = max(0,lambda2);
    
    %更新mu
    mu2=zeros(1,TOTAL_MULTIGROUP);
    for i4M=1:TOTAL_MULTIGROUP
        mu2(1,i4M) = mu(1,i4M)- stepMU_c2*(rateM(1,i4M)-MIN_RATE_Rmin);
        mu2(1,i4M) = max(0,mu2(1,i4M));
    end
    
    % ====================================================
    %        （3）接下来判断次梯度法是否达到收敛条件
    % ====================================================
    % 判断次梯度法是否收敛，由两条件判断：
    % 条件1：
        % [对偶变量变化 < 阈值] && [和功率 < 最大功率] && [速率 > 最小速率]
    % 条件2：
        % 如果迭代次数到达上限，终止次梯度法迭代
        
    % 条件1：
    %判断 [对偶变量变化 < 阈值]
    totalDualChange=0;%记录对偶变量每次变化比率之和
    DualChangeThreshold = 1*(1e-3);%对偶变量变化比率之和的阈值
    
    for iDual=1:TOTAL_MULTIGROUP
%         if mu(1,iDual)>0
            %计算对偶变量每次变化比率之和
            totalDualChange = totalDualChange+sum(abs( (lambda2-lambda)./(lambda+1*(1e-5))))...
                +sum( abs( (mu2(1,iDual)-mu(1,iDual))./( mu(1,iDual)+1*(1e-5) ) ) );
%         end
    end
    totalDualChange
    if totalDualChange < DualChangeThreshold
        isBelowDualChange=1;
    else isBelowDualChange=0;
    end
    
    %判断 [和功率 < 最大功率]
    if sum(powerSub_Pn(:))<MAX_POWER_Pth
        isBelowMAX_POWER=1;
    else isBelowMAX_POWER=0;
    end
    
    %判断 [速率 > 最小速率]
    isAboveMIN_RATE=1;
    isAboveMIN_RATE_M=zeros(1,TOTAL_MULTIGROUP);
    for irateM=1:TOTAL_MULTIGROUP
        if rateM(1,irateM)>MIN_RATE_Rmin;
            isAboveMIN_RATE_M(1,irateM)=1;
        else isAboveMIN_RATE_M(1,irateM)=0;
        end
        isAboveMIN_RATE=isAboveMIN_RATE*isAboveMIN_RATE_M(1,irateM);
    end
    
    %计算收敛标志
    isConverge=isBelowDualChange*isBelowMAX_POWER*isAboveMIN_RATE;% 收敛标志
    
    % 条件2：
    isExceedIteration=0;%达到迭代上限标志
    if timesIteration>=MAX_ITERATIONS
        %由于超出最大迭代次数而结束循环
        isExceedIteration=1;
    end
    
    % ====================================================
    %       （4）接下来根据判断结果来重新计算步长
    % ====================================================
    % 更新步长随着迭代次数的增加而减少，使得对偶变量的变化更加精细化，利于收敛
        
%     % 如果最小速率和最大功率约束都满足，缩小步长
%     if isBelowMAX_POWER==1&&isAboveMIN_RATE==1
%         stepLambda_c1 = 1*(1e-1)/sqrt(timesIteration);%Lambda的更新步长
%         stepMU_c2 = 1*(1e-4)/sqrt(timesIteration);%Mu的更新步长
%         
%     % 如果迭代步长大于3000，进一步缩小步长
%     %   注：此时一般接近于收敛，所以直接把步长缩得更小
%     elseif timesIteration>3000
%         stepLambda_c1 = 1*(1e-2)/sqrt(timesIteration);%Lambda的更新步长
%         stepMU_c2 = 1*(1e-5)/sqrt(timesIteration);%Mu的更新步长
%         
%     % 如果以上约束都不满足，正常计算步长
%     else
%         stepLambda_c1 = 5*(1e-1)/sqrt(timesIteration);%Lambda的更新步长
%         stepMU_c2 = 1*(1e-3)/sqrt(timesIteration);%Mu的更新步长
%     end
    % ====================================================
    %       （5）变量更新
    % ====================================================
    lambda = lambda2
    mu = mu2
    timesIteration=timesIteration+1;
    isScheduleDone=isExceedIteration+isConverge;
    PTS=sum(powerSub_Pn(:))/MAX_POWER_Pth
    RTS=rateM/MIN_RATE_Rmin
    
end % while 大循环结束

