function [ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input )
%XuY_Fun_StepOptimize1_Step1: Step1 of the StepOptimize1
%   给定子载波分配与功率分配下，多播推送用户集合选择

% Author	: XuY
% Modified	: 2014-12-11 15:00

% Create and Modify Date and History :
% - 2015/1/15 
%     - change revenue
%       "tmpRev_3=t*powerUE" to "tmpRev_3=alpha_1*sum(sum(powerSub_Pn))"; 
%     - delete parameter "t" 
%     - add parameter "alpha_1"
% - 2015/4/3 
%     - remove interestUM 

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
TOTAL_MULTIGROUP=OP1_Step1_Input.TOTAL_MULTIGROUP;%多播组总数
TOTAL_USER=OP1_Step1_Input.TOTAL_USER;%用户总数
TOTAL_SUB=OP1_Step1_Input.TOTAL_SUB;%子载波总数
BAND_SUB=OP1_Step1_Input.BAND_SUB;%子载波带宽，一般现在取0.3125MHz
gainChannel=OP1_Step1_Input.gainChannel;
alpha_1=OP1_Step1_Input.alpha_1;
alpha_2=OP1_Step1_Input.alpha_2;
xi=OP1_Step1_Input.xi;
interestUM = OP1_Step1_Input.interestUM;
% alpha_2=OP1_Step1_Input.alpha_2;

% 迭代参数
scheduleSub_rho=OP1_Step1_Input.scheduleSub_rho;
powerSub_Pn=OP1_Step1_Input.powerSub_Pn;

% 约束参数
MAX_POWER_Pth=OP1_Step1_Input.MAX_POWER_Pth;%功率约束 单位W

%-----------------函数局部变量----------------------
punishBeta=10;%惩罚函数的系数
% powerUE=1*(1e-3);%用户端接收功耗，单位:W
% powerBS=1*(1e-1);%基站功耗，单位:W

% --------------------------------------------
%      初始化每个多播组集合中初始最差接收SINR
% --------------------------------------------
% 初始化的时候用一个很大的值，后面用户在迭代过程中加入集合时
% 再根据加入的用户的SINR对集合进行更新；
minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
tmp_minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);

% --------------------------------------------
%           初始化业务多播推送集合
% --------------------------------------------
pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
%                 备注：
% 如果一个用户可以被分到多个多播组，那么当最差信道的用户被分到全部多播组中时
% 可能造成几个多播组的minSNR相同，而导致把所有的子载波都分配给一个多播组
% pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% for ipushU=1:TOTAL_USER
%     %利用随机序列保证每行只有一个1，即每个子载波只能被分配到一个多播组
%     index_push=randperm(TOTAL_MULTIGROUP);
%     pushUM(ipushU,index_push(1))=1;
% end


% --------------------------------------------
%           初始化收益矩阵
% --------------------------------------------
revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
for iRevU=1:TOTAL_USER
    for iRevM=1:TOTAL_MULTIGROUP
        tmpRev_1=0;
        for iRevN=1:TOTAL_SUB
            %多播组和速率
            tmpRev_1=tmpRev_1+BAND_SUB*scheduleSub_rho(iRevN,iRevM)...
               *log2(1+powerSub_Pn(iRevN,iRevM)...
                *min(gainChannel(iRevU,iRevN,1),minSNR_gamma(iRevN,iRevM)));    
        end
        %用户对接受业务不感兴趣的惩罚函数
        tmpRev_2=punishBeta*(1-interestUM(iRevU,iRevM));
        %用户功耗
        tmpRev_3=alpha_1*sum(sum(powerSub_Pn));
        
        %收益=和速率-兴趣惩罚函数-设备固有功耗
        revenueUM(iRevU,iRevM)=tmpRev_1-tmpRev_2-tmpRev_3;
    end
end

%******************************************************
%                迭代用户选择
%******************************************************
isIterationDone=0;
while isIterationDone==0
% pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% 选择收益最大的用户—业务对(k*,m*)
[maxRev_U,maxRev_M]=find(revenueUM==max(max(revenueUM)));

% 把pushUM矩阵的每行加和，用于判断该用户是否已经接受了一个多播组业务了
pushUM_SUM=sum(pushUM,2);

% 若revenueUM(k*,m*)>0,将用户k*,加入业务m*的多播推送用户集合
if revenueUM(maxRev_U,maxRev_M)>0
    % 如果该用户已经被分配过多播组业务了
    if pushUM_SUM(maxRev_U,1)==1
        % 把收益标记为零，防止重复选择
        % 换句话说，一个用户只能接收一个多播业务
        revenueUM(maxRev_U,maxRev_M)=0;
    else      
        % 计算加入用户k*之前，多播组m*的总收益
        revPrevious=sum(revenueUM(:,maxRev_M));
        
        % 更新多播推送用户集合中最差接收 SINR
        % minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
            for iSNR_N=1:TOTAL_SUB
                if gainChannel(maxRev_U,iSNR_N,1)<minSNR_gamma(iSNR_N,maxRev_M)
                    minSNR_gamma(iSNR_N,maxRev_M)=gainChannel(maxRev_U,iSNR_N,1);
                end
            end
%此时不确定是否加入该用户，因而收益信息不能直接写入revenueUM，用临时变量来记录更新后的收益用于比较        
    temp_revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
        % 更新收益矩阵
        for iRevU=1:TOTAL_USER
            for iRevM=1:TOTAL_MULTIGROUP
                tmpRev_1=0;
                if revenueUM(iRevU,iRevM)==0
                    continue
                else
                    for iRevN=1:TOTAL_SUB
                        %多播组和速率
                        tmpRev_1=tmpRev_1+BAND_SUB*scheduleSub_rho(iRevN,iRevM)...
                            *log2(1+powerSub_Pn(iRevN,iRevM)...
                            *min(gainChannel(iRevU,iRevN,1),minSNR_gamma(iRevN,iRevM)));    
                    end
                    %用户对接受业务不感兴趣的惩罚函数
                    tmpRev_2=punishBeta*(1-interestUM(iRevU,iRevM));
                    %用户设备端固有功耗
                    tmpRev_3=alpha_1*sum(sum(powerSub_Pn));
                    %收益=和速率-兴趣惩罚函数-设备固有功耗
                    temp_revenueUM(iRevU,iRevM)=tmpRev_1-tmpRev_2-tmpRev_3;
                end
            end
        end
                       
    end % if pushUM(maxRev_U,maxRev_M)==1
    
    % 计算加入用户k*之后，多播组m*的总收益
    revAfter=sum(revenueUM(:,maxRev_M));
    
    % 如果加入新用户拉低多播组SINR从而导致整个多播组的收益下降的话
    % 那么不加入该新用户
    if revAfter >= revPrevious
        pushUM(maxRev_U,maxRev_M)=1;
        %将更新后的收益信息写入收益矩阵revenueUM
        revenueUM = temp_revenueUM;
    else
        pushUM(maxRev_U,maxRev_M)=0;
    end
    
end % if revenueUM(maxRev_U,maxRev_M)>0
revenueUM(maxRev_U,maxRev_M)=0;

% 判断迭代是否完成
 % 如果revenueUM中只剩下0和负数那么迭代完成
    if max(max(revenueUM))>0;
        isIterationDone=0;
    else
        isIterationDone=1;
    end

end % while

% 生成多播组兴趣矩阵（omegaM）
[ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM, interestUM, alpha_2, xi );

% *************************************************************************
%                   输出参数结构体
% *************************************************************************
OP1_Step1_Output.omegaM=omegaM;
OP1_Step1_Output.minSNR_gamma=minSNR_gamma;
OP1_Step1_Output.pushUM=pushUM;

end %Function