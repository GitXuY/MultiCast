function [ Min_Output ] = XuY_Fun_UpperLimit_Min( Min_Input )
%XuY_Fun_UpperLimit_Min: minimize the maximal objective function
%   目标函数外部min迭代

% Author	: XuY
% Modified	: 2015-1-10

% Create and Modify Date and History :
% - 

% Error Case :
% -

% Call 
% -

% Called
% -XuY_Scp_Main

% *************************************************************************
%                       参量设置
% *************************************************************************
%-------------------输入参数------------------------
% 场景参数
TOTAL_MULTIGROUP=XuY_Fun_UpperLimit_Min.TOTAL_MULTIGROUP;%多播组总数
TOTAL_USER=XuY_Fun_UpperLimit_Min.TOTAL_USER;%用户总数
TOTAL_SUB=XuY_Fun_UpperLimit_Min.TOTAL_SUB;%子载波总数
BAND_SUB=XuY_Fun_UpperLimit_Min.BAND_SUB;%子载波带宽，一般现在取0.3125MHz

%-----------------函数局部变量----------------------
%拉格朗日对偶乘子初始化
lambda=2;%对偶变量lambda
mu=ones(1,TOTAL_MULTIGROUP);%对偶变量mu

%迭代参数
isScheduleDone=0;% 资源分配完成标志
MAX_ITERATIONS=10000;%最大迭代次数
timesIteration=1;%初始化迭代计数器

%迭代资源分配的大循环
while isScheduleDone==0
    
% *************************************************************************
%                       阶段1：Pn和D(n,m)的计算
% *************************************************************************
%   存储信息矩阵初始化
    funcD=zeros(1,TOTAL_MULTIGROUP);%D(n,m) 1*多播组数
    rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%每个子载波的传输速率
    powerSub_Pn=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%每个子载波的功率
    powerSub_Pn_Cal=zeros(1,TOTAL_MULTIGROUP);%计算时存储Pn值的临时变量
      
    rateM=zeros(1,TOTAL_MULTIGROUP);%每个多播组的和速率
    powerM=zeros(1,TOTAL_MULTIGROUP);%每个多播组的和功率
    
    isSubSchedule=zeros(1,TOTAL_SUB);%记录子载波是否已经被分配出去
    scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%子载波的分配函数
    
%   迭代计算    
    for iN=1:TOTAL_SUB
        for iiM=1:TOTAL_MULTIGROUP
            % 如果该子载波还未被分配
            if isSubSchedule(iN)==0
                % 子载波iN被分配到第iiM个多播组时的最优化功率
                powerSub_Pn_Cal(1,iiM)=(mu(1,iiM)+omegaM(1,iiM))...
                    *BAND_SUB./((alpha_1+lambda)*log(2))-1./minSNR_gamma(iN,iiM);
                
                % 功率取正数
                powerSub_Pn_Cal(1,iiM)=max(powerSub_Pn_Cal(1,iiM),0);
                
                % 对应该功率值，所得到的D(n,m)值
                funcD(1,iiM)=(mu(1,iiM)+omegaM(1,iiM))*BAND_SUB*...
                    log2(1+powerSub_Pn_Cal(1,iiM)*minSNR_gamma(iN,iiM))-...
                    (alpha_1+lambda)*powerSub_Pn_Cal(1,iiM);
            % 如果该子载波已经被分配
            else
                % 给该子载波“临时”功率置零
                powerSub_Pn_Cal(1,iiM)=0;
                % funcD置-1000确保该子载波不被选中
                funcD(1,iiM)=-1000;
            end % if isSubSchedule
        end % for iiM
        
        % 对第n个子载波，计算该子载波被分配到哪一个多播组的时候
        % 所得到的效益最大（D(n,m)取到最大）
        bestMultiGroupForiN = find( funcD>=max(funcD) );
        
        % 如果找到了两个以上的ms,即最大的有两个，那么随机选一个
        if length(bestMultiGroupForiN) > 1
            bestMultiGroupForiN...
                            = round(length(bestMultiGroupForiN*rand(1)));
        end
        
        % 把该子载波分配给这个多播组
        scheduleSub_rho(iN,bestMultiGroupForiN) = 1;
        
        % 标明该子载波已经被分配出去，防止重复分配
        isSubSchedule(1,iN)=1;
        
        % 把该子载波在“临时”子载波功率矩阵中的值登记到子载波功率矩阵
        powerSub_Pn(iN,bestMultiGroupForiN)...
                                = powerSub_Pn_Cal(1,bestMultiGroupForiN);
        
        % 计算该子载波的速率
        rateSub(iN,bestMultiGroupForiN)...
                    =BAND_SUB*log2(1+powerSub_Pn(iN,bestMultiGroupForiN)...
                            *minSNR_gamma(iN,bestMultiGroupForiN));
        
    end % for iN
    
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
    if timesIteration==1
        %Lambda的更新步长初始化
        stepLambda_c1= 5*(1e-1);
        
        %Mu的更新步长初始化
        stepMU_c2 = 1*(1e-3);
    end
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
    DualChangeThreshold = 5*(1e-5);%对偶变量变化比率之和的阈值
    
    for iDual=1:TOTAL_MULTIGROUP
        if mu(1,iDual)>0
            %计算对偶变量每次变化比率之和
            totalDualChange = totalDualChange+sum(abs( (lambda2-lambda)./(lambda)))...
                +sum(abs( (mu2(1,iDual)-mu(1,iDual))./(mu(1,iDual))));
        end
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
    if isBelowMAX_POWER==1&&isAboveMIN_RATE==1
        stepLambda_c1 = 1*(1e-1)/sqrt(timesIteration);%Lambda的更新步长
        stepMU_c2 = 1*(1e-4)/sqrt(timesIteration);%Mu的更新步长
        
    % 如果迭代步长大于3000，进一步缩小步长
    %   注：此时一般接近于收敛，所以直接把步长缩得更小
    elseif timesIteration>3000
        stepLambda_c1 = 1*(1e-2)/sqrt(timesIteration);%Lambda的更新步长
        stepMU_c2 = 1*(1e-5)/sqrt(timesIteration);%Mu的更新步长
        
    % 如果以上约束都不满足，正常计算步长
    else
        stepLambda_c1 = 5*(1e-1)/sqrt(timesIteration);%Lambda的更新步长
        stepMU_c2 = 1*(1e-3)/sqrt(timesIteration);%Mu的更新步长
    end
    % ====================================================
    %       （5）变量更新
    % ====================================================
    lambda = lambda2;
    mu = mu2;
    timesIteration=timesIteration+1;
    isScheduleDone=isExceedIteration+isConverge;
    PTS=sum(powerSub_Pn(:))/MAX_POWER_Pth;
    RTS=rateM/MIN_RATE_Rmin;
    
end % while 大循环结束

% *************************************************************************
%                   输出参数结构体
% *************************************************************************
% 资源分配结果
OP1_Step2_Output.powerSub_Pn=powerSub_Pn;
OP1_Step2_Output.rateSub=rateSub;
OP1_Step2_Output.rateM=rateM;
OP1_Step2_Output.powerM=powerM;
OP1_Step2_Output.scheduleSub_rho=scheduleSub_rho;
OP1_Step2_Output.PTS=PTS;
OP1_Step2_Output.RTS=RTS;

% 拉格朗日对偶乘子结果
OP1_Step2_Output.lambda=lambda;
OP1_Step2_Output.mu=mu;

% 收敛信息
OP1_Step2_Output.timesIteration=timesIteration;
OP1_Step2_Output.isExceedIteration=isExceedIteration;
OP1_Step2_Output.isAboveMIN_RATE=isAboveMIN_RATE;
OP1_Step2_Output.isBelowMAX_POWER=isBelowMAX_POWER;
OP1_Step2_Output.totalDualChange=totalDualChange;
OP1_Step2_Output.isBelowDualChange=isBelowDualChange;

end %Function
