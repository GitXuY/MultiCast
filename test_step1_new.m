clc
clear
load gainChannel
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
% 由多播推送集合中用户不愿意接收该业务造成的收益损失
betaM=10*ones(1,TOTAL_MULTIGROUP);
% 收益矩阵
revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% 暂存收益矩阵
temp_revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% 每个用户充当基准用户时的最佳可选多播推送用户集合
pushUUM=zeros(TOTAL_USER,TOTAL_USER,TOTAL_MULTIGROUP); 
% powerUE=1*(1e-3);%用户端接收功耗，单位:W
% powerBS=1*(1e-1);%基站功耗，单位:W

%-------------------step2参数------------------------
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

%-------------------输出参数------------------------
% omegaM;
% 基准用户集合
criterionUserM=zeros(1,TOTAL_MULTIGROUP);
% 最大收益
maxRev=zeros(1,TOTAL_MULTIGROUP);
% 每个多播组集合中最差接收SINR
minSNR_gamma=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
% 业务多播推送集合
pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
%******************************************************
%                迭代用户选择
%******************************************************
LOOP=3;
m=0;
for iM=1:TOTAL_MULTIGROUP
    %遍历多播业务组
    for iU=1:TOTAL_USER
        %遍历所有用户群体，每个用户充当一次基准用户
        criterionUser=iU;%设置基准用户
        isAdd=ones(1,TOTAL_USER); %标记用户是否能加入集合
        % 基准用户的平均SNR
        iU_average_SNR = 0;
        for iN=1:TOTAL_SUB
            if scheduleSub_rho(iN,iM)==1
            %遍历分配到该多播组的子载波
                iU_average_SNR = iU_average_SNR + gainChannel(criterionUser,iN,LOOP);
            end
        end
        iU_average_SNR = iU_average_SNR / TOTAL_SUB;          
        for i2U=1:TOTAL_USER
            %对每个基准用户，选择最佳可选多播推送用户集合
            %--------------------------------------------------
            %                   i2U平均SNR
            %--------------------------------------------------
            i2U_average_SNR = 0;
            for iN=1:TOTAL_SUB
                if scheduleSub_rho(iN,iM)==1
                %遍历分配到该多播组的子载波
                    i2U_average_SNR = i2U_average_SNR + gainChannel(i2U,iN,LOOP);
                end
            end
            i2U_average_SNR = i2U_average_SNR / TOTAL_SUB;  
            %--------------------------------------------------
            %     i2U平均SNR比基准用户均SNR小的不能加入集合
            %--------------------------------------------------
            if i2U_average_SNR<iU_average_SNR
                isAdd(1,i2U)=0;
            end
            %--------------------------------------------------
            %           收益小于0的不能加入集合
            %--------------------------------------------------
            tmp1=0;
            for i2N=1:TOTAL_SUB
                %多播组和速率
                tmp1=tmp1+BAND_SUB*scheduleSub_rho(i2N,iM)...
                    *log2(1+powerSub_Pn(i2N,iM)...
                    *gainChannel(criterionUser,i2N,LOOP));
            end
            %用户对接受业务不感兴趣的惩罚函数
            tmp2=interestUM(i2U,iM)-alpha_2*xi;
            %用户设备端固有功耗
            tmp3=betaM(1,iM)*(1-interestUM(i2U,iM));
            %收益=和速率-兴趣惩罚函数-设备固有功耗
            temp_revenueUM(i2U,iM)=tmp1*tmp2-tmp3;
            if temp_revenueUM(i2U,iM) < 0
                isAdd(1,i2U)=0;
            end    
        end % for i2U=1:TOTAL_USER
        %--------------------------------------------------
        %  加入新用户后更新多播推送用户集合中最差接收 SINR
        %--------------------------------------------------
        %如果加入了新用户
        if isAdd(1,i2U)==1
            % 遍历子载波
            for iN=1:TOTAL_SUB
               % 仅更新该多播组分配到的子载波
               if scheduleSub_rho(iN,iM)==1
                   if minSNR_gamma(iN,iM)==0o
                      minSNR_gamma(iN,iM)=gainChannel(i2U,iN,1); 
                   elseif gainChannel(i2U,iN,1)<minSNR_gamma(iN,iM)
                        minSNR_gamma(iN,iM)=gainChannel(i2U,iN,1);
                   end
               end
            end
        end
        %--------------------------------------------------
        %   iU用户充当基准用户时,最佳可选多播推送用户集合
        %--------------------------------------------------
        pushUUM(iU,:,iM)=isAdd(:);
        %--------------------------------------------------
        % iU用户充当基准用户时,最佳可选多播推送用户集合的收益
        %--------------------------------------------------
        for i3U=1:TOTAL_USER
            tmp1_1=0;
            if isAdd(1,i3U)==1
                for i3N=1:TOTAL_SUB
                %多播组和速率
                tmp1_1=tmp1_1+BAND_SUB*scheduleSub_rho(i3N,iM)...
                    *log2(1+powerSub_Pn(i3N,iM)...
                    *gainChannel(criterionUser,i3N,LOOP));
                end
                %用户对接受业务不感兴趣的惩罚函数
                tmp1_2=interestUM(i3U,iM)-alpha_2*xi;
                %用户设备端固有功耗
                tmp1_3=betaM(1,iM)*(1-interestUM(i3U,iM));
                %收益=和速率-兴趣惩罚函数-设备固有功耗
                revenueUM(iU,iM)=revenueUM(iU,iM)+tmp1_1*tmp1_2-tmp1_3;
            end
        end
        
    end % for iU=1:TOTAL_USER
    
    [maxRev(1,iM), criterionUserM(1,iM)] = max(revenueUM(:,iM));
    pushUM(:,iM)=pushUUM(criterionUserM(1,iM),:,iM);
    
end
% 生成多播组兴趣矩阵（omegaM）
[ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM, interestUM, alpha_2, xi );

save test_step1_new
