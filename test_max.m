clc
clear
load gainChannel
% *************************************************************************
%                       参量设置
% *************************************************************************
TOTAL_MULTIGROUP=2;%多播组总数
TOTAL_USER=100;%用户总数
TOTAL_SUB=64;%子载波总数
BAND_SUB= 0.3125; %子载波带宽，一般现在取0.3125MHz
MAX_POWER_Pth=1;%功率约束 单位W
MIN_RATE_Rmin=30;% 最小发送速率约束 单位？

alpha_1=1;%基站功耗对效用的权重
alpha_2=1;%用户功耗对效用的权重
xi=0.001;%表示用户没接收单位数据业务的额外功耗
betaM=ones(1,TOTAL_MULTIGROUP);% 由多播推送集合中用户不愿意接收该业务造成的收益损失
% LOOP=10;%信道循环次

%用户兴趣矩阵，用户对哪个多播组感兴趣
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);
lambda=2;
mu=ones(1,TOTAL_MULTIGROUP);%对偶变量mu

%-------------------待计算参数------------------------
% 根据基准用户的参数，子载波的功分
powerSub_Pn = zeros(1,TOTAL_SUB);
% 根据基准用户的参数，能被分配到多播组的用户集合
userSet_Knm  = zeros(TOTAL_USER,TOTAL_MULTIGROUP,TOTAL_SUB);
% D(n,m) 1*多播组数
funcD=0;%D(n,m) 1*多播组数
funcD_plus=0;
% 优化变量
tmp_betaM = TOTAL_MULTIGROUP .* betaM / TOTAL_SUB;
betaMN = repmat(tmp_betaM',1,TOTAL_SUB);
% 记录多播组的收益
revenueMN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);
% 记录多播组的基准用户
criterionUSER_MN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);
% 记录多播组的基准SNR
criterionSNR_MN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);
% 记录每个用户充当基准用户时 多播组的收益
revenueUMN = zeros(TOTAL_USER,TOTAL_MULTIGROUP, TOTAL_SUB);
% 记录每个用户充当基准用户时 多播组的基准SNR
criterionSNR_UMN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);
% 记录每个子载波可获得的最大收益
maxRevenue = zeros(1,TOTAL_SUB);
% 记录每个子载波对应的最优多播组
bestMG = zeros(1,TOTAL_SUB);
% 记录子载波分配
scheduleSub_rho = zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
%每个子载波的传输速率
rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
%每个多播组的和速率
rateM=zeros(1,TOTAL_MULTIGROUP);
%每个多播组的和功率
powerM=zeros(1,TOTAL_MULTIGROUP);

%-------------------迭代参数------------------------
isScheduleDone=0;% 资源分配完成标志
MAX_ITERATIONS=3000;%最大迭代次数
timesIteration=1;%初始化迭代计数器
while isScheduleDone~=1
%迭代资源分配的大循环
% *************************************************************************
%     阶段1：假设子载波n 分配给业务m ，计算此时最佳的功分与推送用户集合
% *************************************************************************
for iSub=1:TOTAL_SUB
% 内部 max 问题可以拆分为N 个子问题，其中第 n 个子问题负责子载波n 的业务选
% 择、功分以及相应的多播推送用户选择
if sum(scheduleSub_rho(iSub,:))==0

    for iMG=1:TOTAL_MULTIGROUP
    % 假设子载波n 分配给业务m ，计算此时最佳的功分与推送用户集合,遍历 M 种可能的
    % 子载波分配策略，选出最佳的子载波分配
%      betaMN (1) = betaMN + TOTAL_MULTIGROUP * betaM(iMG) / TOTAL_SUB;

        for iUserCSNR = 1:TOTAL_USER
        % 按顺序依次选择用户作为基准用户，然后选出能使效益最大的用户作为最佳基准用户。
        
            userSet_Knm(:,iMG,iSub) = 0;    
            criterionSNR_UMN(iUserCSNR, iMG, iSub) = gainChannel(iUserCSNR,iSub,2);   %第三个参数为Loop
            % =============================================================
            %   当基准用户为 iUserCSNR 时 传输用户集合userSet_Knm
            % =============================================================
            for iUser = 1: TOTAL_USER  
                % 根据基准用户选出iMG的传输集合
                if gainChannel(iUser,iSub,1) >= criterionSNR_UMN(iUserCSNR, iMG, iSub)... 
                     && interestUM(iUser,iMG)-alpha_2*xi+mu(1,iMG) >= 0            
                    userSet_Knm(iUser,iMG,iSub) = 1;
                end  
            end 
            % =============================================================
            %     当基准用户为 iUserCSNR 时 子载波iSub功分powerSub_Pn
            % =============================================================
            temp_p1 = interestUM(:,iMG)-alpha_2*xi+mu(1,iMG);
            temp_p2 = userSet_Knm(:,iMG,iSub).* temp_p1;
            temp_p2 = temp_p2 / 300; 
            powerSub_Pn(iSub)=(BAND_SUB * sum(temp_p2(:)))...
            /(log(2)*(alpha_1+lambda))-1/criterionSNR_UMN(iUserCSNR, iMG, iSub);
            powerSub_Pn(iSub)=max(powerSub_Pn(iSub),0);
            % =============================================================
            %       当基准用户为 iUserCSNR 时 多播组的传输收益
            % =============================================================
            for iUser = 1: TOTAL_USER  
                if userSet_Knm(iUser,iMG,iSub) == 1              
                    temp_d1 = interestUM(iUser, iMG)-alpha_2*xi+mu(1,iMG);
                    temp_d2 = log2(1 + powerSub_Pn(iSub) * criterionSNR_UMN(iUserCSNR, iMG, iSub));
                    temp_d3 = betaMN(iMG,iSub) * (1 - interestUM(iUser, iMG));
                    funcD   = temp_d1 * BAND_SUB * temp_d2 - temp_d3;
                    funcD_plus = funcD_plus + funcD;
                end 
            end 

            revenueUMN(iUserCSNR,iMG, iSub) = funcD_plus - (alpha_1 + lambda) * powerSub_Pn(iSub);
            funcD_plus = 0;
            
        end % for iUserCSNR = 1:TOTAL_USER
        
        %选出使得多播组受益最大的基准用户 及其 对应的多播组收益
        [revenueMN(iMG, iSub), criterionUSER_MN(iMG, iSub)]= max(revenueUMN(:,iMG, iSub));
        criterionSNR_MN(iMG, iSub) = criterionSNR_UMN(criterionUSER_MN(iMG, iSub),iMG, iSub);
        
    end % for iMG=1:TOTAL_MULTIGROUP
    
    [ maxRevenue(iSub), bestMG(iSub) ]= max(revenueMN(:, iSub));
    
    scheduleSub_rho( iSub, bestMG(iSub)) = 1;
    
    % 计算该子载波的速率
    rateSub(iSub,bestMG(iSub))...
                    =BAND_SUB*log2(1+powerSub_Pn(iSub)...
                            *criterionSNR_MN(bestMG(iSub), iSub));
    
end
end % for iSub=1:TOTAL_SUB
        
    %计算每个多播组的和速率                    
    for irateM=1:TOTAL_MULTIGROUP
        rateM(1,irateM)=sum(rateSub(:,irateM));
    end 
        
    %计算每个多播组的和功率
    for ipowerM=1:TOTAL_MULTIGROUP
        powerM(1,ipowerM)=sum(powerSub_Pn(:,ipowerM));
    end
% *************************************************************************
%                   阶段2：对偶变量子梯度法迭代
% *************************************************************************
    % ====================================================
    %     （1）第一次迭代时使用初始化值来启动迭代过程
    % ====================================================
%     if timesIteration==1
        %Lambda的更新步长初始化
        stepLambda_c1= 5*(1e-1)/sqrt(timesIteration);
        
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
                +sum(abs( (mu2(1,iDual)-mu(1,iDual))./(mu(1,iDual)+1*(1e-5))));
%         end
    end
    
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
        
    % 如果最小速率和最大功率约束都满足，缩小步长
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
    for iM=1:TOTAL_MULTIGROUP
        if sum(scheduleSub_rho(:,iM))~=0
            betaMN(iM,:) = betaM(iM)/sum(scheduleSub_rho(:,iM));
        end
    end
    
    lambda = lambda2
    mu = mu2
    timesIteration=timesIteration+1
    isScheduleDone=isExceedIteration+isConverge;
    PTS=sum(powerSub_Pn(:))/MAX_POWER_Pth
    RTS=rateM/MIN_RATE_Rmin

end