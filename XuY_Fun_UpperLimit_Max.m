function [ Max_Output ] = XuY_Fun_UpperLimit_Max( Max_Input )
%XuY_Fun_UpperLimit_Max: maximize objective function
%   目标函数内部max迭代

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
%   场景参数
TOTAL_MULTIGROUP    = Max_Input.TOTAL_MULTIGROUP;%多播组总数
TOTAL_USER          = Max_Input.TOTAL_USER;%用户总数
TOTAL_SUB           = Max_Input.TOTAL_SUB;%子载波总数
BAND_SUB            = Max_Input.BAND_SUB;%子载波带宽，一般现在取0.3125MHz
gainChannel         = Max_Input.gainChannel;
interestUM          = Max_Input.interestUM;% rand(TOTAL_USER,TOTAL_MULTIGROUP);
alpha_1             = Max_Input.alpha_1;
alpha_2             = Max_Input.alpha_2;
xi                  = Max_Input.xi;
betaM               = Max_Input.betaM;
%   迭代参数
lambda              = Max_Input.lamda;
mu                  = Max_Input.mu;

%-------------------待计算参数------------------------
% 根据基准用户的参数，子载波的功分
powerSub_Pn = zeros(1,TOTAL_SUB);

% 根据基准用户的参数，能被分配到多播组的用户集合
userSet_Knm  = zeros(TOTAL_USER,TOTAL_MULTIGROUP,TOTAL_SUB);

% D(n,m) 1*多播组数
funcD=0;%D(n,m) 1*多播组数
funcD_plus=0;

% 优化变量
betaMN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);

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



% *************************************************************************
%                       阶段1：
% *************************************************************************
for iSub=1:TOTAL_SUB
% 内部 max 问题可以拆分为N 个子问题，其中第 n 个子问题负责子载波n 的业务选
% 择、功分以及相应的多播推送用户选择

    for iMG=1:TOTAL_MULTIGROUP
    % 假设子载波n 分配给业务m ，计算此时最佳的功分与推送用户集合,遍历 M 种可能的
    % 子载波分配策略，选出最佳的子载波分配
    betaMN = betaMN + TOTAL_MULTIGROUP * betaM(iMG) / TOTAL_SUB;

        for iUserCSNR = 1:TOTAL_USER
        % 按顺序依次选择用户作为基准用户，然后选出能使效益最大的用户作为最佳基准用户。
        
            userSet_Knm(:,iMG,iSub) = 0;    
            criterionSNR_UMN(iUserCSNR, iMG, iSub) = gainChannel(iUserCSNR,iSub,2);   %第三个参数为Loop
            % =============================================================
            %   当基准用户为 iUserCSNR 时 传输用户集合userSet_Knm
            % =============================================================
            for iUser = 1: TOTAL_USER  
                % 根据基准用户选出iMG的传输集合
                if gainChannel(iUser,iSub,1) >= criterionSNR_UMN(iUserCSNR, iMG, iSub) && interestUM(iUser,iMG)-alpha_2*xi+mu >= 0            
                    userSet_Knm(iUser,iMG,iSub) = 1;
                end  
            end 
            % =============================================================
            %     当基准用户为 iUserCSNR 时 子载波iSub功分powerSub_Pn
            % =============================================================
            temp_p1 = interestUM-alpha_2*xi+mu;
            temp_p2 = userSet_Knm(:,:,iSub).* temp_p1;
            powerSub_Pn(iSub)=(BAND_SUB * sum(temp_p2(:)))...
            /(log(2)*(alpha_1+lambda))-1/criterionSNR_UMN(iUserCSNR, iMG, iSub);                                  
            % =============================================================
            %       当基准用户为 iUserCSNR 时 多播组的传输收益
            % =============================================================
            for iUser = 1: TOTAL_USER  
                if userSet_Knm(iUser,iMG,iSub) == 1              
                    temp_d1 = interestUM(iUser, iMG)-alpha_2*xi+mu;
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
    
end % for iSub=1:TOTAL_SUB
% *************************************************************************
%                   输出参数结构体
% *************************************************************************
Max_Output.powerSub_Pn      = powerSub_Pn;
Max_Output.userSet_Knm      = userSet_Knm;
Max_Output.funcD            = funcD;%D(n,m) 1*多播组数
Max_Output.betaMN           = betaMN;
Max_Output.revenueMN        = revenueMN;
Max_Output.criterionUSER_MN = criterionUSER_MN;
Max_Output.revenueUMN       = revenueUMN;
Max_Output.maxRevenue       = maxRevenue;
Max_Output.bestMG           = bestMG;
Max_Output.criterionSNR_MN  = criterionSNR_MN;
end % Function

