function [ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input )
%XuY_Fun_StepOptimize1_Step1: Step1 of the StepOptimize1
%   �������ز������빦�ʷ����£��ಥ�����û�����ѡ��

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
%                       ��������
% *************************************************************************
%-------------------�������------------------------
% ��������
TOTAL_MULTIGROUP=OP1_Step1_Input.TOTAL_MULTIGROUP;%�ಥ������
TOTAL_USER=OP1_Step1_Input.TOTAL_USER;%�û�����
TOTAL_SUB=OP1_Step1_Input.TOTAL_SUB;%���ز�����
BAND_SUB=OP1_Step1_Input.BAND_SUB;%���ز�����һ������ȡ0.3125MHz
gainChannel=OP1_Step1_Input.gainChannel;
alpha_1=OP1_Step1_Input.alpha_1;
alpha_2=OP1_Step1_Input.alpha_2;
xi=OP1_Step1_Input.xi;
interestUM = OP1_Step1_Input.interestUM;
% alpha_2=OP1_Step1_Input.alpha_2;

% ��������
scheduleSub_rho=OP1_Step1_Input.scheduleSub_rho;
powerSub_Pn=OP1_Step1_Input.powerSub_Pn;

% Լ������
MAX_POWER_Pth=OP1_Step1_Input.MAX_POWER_Pth;%����Լ�� ��λW

%-----------------�����ֲ�����----------------------
punishBeta=10;%�ͷ�������ϵ��
% powerUE=1*(1e-3);%�û��˽��չ��ģ���λ:W
% powerBS=1*(1e-1);%��վ���ģ���λ:W

% --------------------------------------------
%      ��ʼ��ÿ���ಥ�鼯���г�ʼ������SINR
% --------------------------------------------
% ��ʼ����ʱ����һ���ܴ��ֵ�������û��ڵ��������м��뼯��ʱ
% �ٸ��ݼ�����û���SINR�Լ��Ͻ��и��£�
minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
tmp_minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);

% --------------------------------------------
%           ��ʼ��ҵ��ಥ���ͼ���
% --------------------------------------------
pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
%                 ��ע��
% ���һ���û����Ա��ֵ�����ಥ�飬��ô������ŵ����û����ֵ�ȫ���ಥ����ʱ
% ������ɼ����ಥ���minSNR��ͬ�������°����е����ز��������һ���ಥ��
% pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% for ipushU=1:TOTAL_USER
%     %����������б�֤ÿ��ֻ��һ��1����ÿ�����ز�ֻ�ܱ����䵽һ���ಥ��
%     index_push=randperm(TOTAL_MULTIGROUP);
%     pushUM(ipushU,index_push(1))=1;
% end


% --------------------------------------------
%           ��ʼ���������
% --------------------------------------------
revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
for iRevU=1:TOTAL_USER
    for iRevM=1:TOTAL_MULTIGROUP
        tmpRev_1=0;
        for iRevN=1:TOTAL_SUB
            %�ಥ�������
            tmpRev_1=tmpRev_1+BAND_SUB*scheduleSub_rho(iRevN,iRevM)...
               *log2(1+powerSub_Pn(iRevN,iRevM)...
                *min(gainChannel(iRevU,iRevN,1),minSNR_gamma(iRevN,iRevM)));    
        end
        %�û��Խ���ҵ�񲻸���Ȥ�ĳͷ�����
        tmpRev_2=punishBeta*(1-interestUM(iRevU,iRevM));
        %�û�����
        tmpRev_3=alpha_1*sum(sum(powerSub_Pn));
        
        %����=������-��Ȥ�ͷ�����-�豸���й���
        revenueUM(iRevU,iRevM)=tmpRev_1-tmpRev_2-tmpRev_3;
    end
end

%******************************************************
%                �����û�ѡ��
%******************************************************
isIterationDone=0;
while isIterationDone==0
% pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% ѡ�����������û���ҵ���(k*,m*)
[maxRev_U,maxRev_M]=find(revenueUM==max(max(revenueUM)));

% ��pushUM�����ÿ�мӺͣ������жϸ��û��Ƿ��Ѿ�������һ���ಥ��ҵ����
pushUM_SUM=sum(pushUM,2);

% ��revenueUM(k*,m*)>0,���û�k*,����ҵ��m*�Ķಥ�����û�����
if revenueUM(maxRev_U,maxRev_M)>0
    % ������û��Ѿ���������ಥ��ҵ����
    if pushUM_SUM(maxRev_U,1)==1
        % ��������Ϊ�㣬��ֹ�ظ�ѡ��
        % ���仰˵��һ���û�ֻ�ܽ���һ���ಥҵ��
        revenueUM(maxRev_U,maxRev_M)=0;
    else      
        % ��������û�k*֮ǰ���ಥ��m*��������
        revPrevious=sum(revenueUM(:,maxRev_M));
        
        % ���¶ಥ�����û������������� SINR
        % minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
            for iSNR_N=1:TOTAL_SUB
                if gainChannel(maxRev_U,iSNR_N,1)<minSNR_gamma(iSNR_N,maxRev_M)
                    minSNR_gamma(iSNR_N,maxRev_M)=gainChannel(maxRev_U,iSNR_N,1);
                end
            end
%��ʱ��ȷ���Ƿ������û������������Ϣ����ֱ��д��revenueUM������ʱ��������¼���º���������ڱȽ�        
    temp_revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
        % �����������
        for iRevU=1:TOTAL_USER
            for iRevM=1:TOTAL_MULTIGROUP
                tmpRev_1=0;
                if revenueUM(iRevU,iRevM)==0
                    continue
                else
                    for iRevN=1:TOTAL_SUB
                        %�ಥ�������
                        tmpRev_1=tmpRev_1+BAND_SUB*scheduleSub_rho(iRevN,iRevM)...
                            *log2(1+powerSub_Pn(iRevN,iRevM)...
                            *min(gainChannel(iRevU,iRevN,1),minSNR_gamma(iRevN,iRevM)));    
                    end
                    %�û��Խ���ҵ�񲻸���Ȥ�ĳͷ�����
                    tmpRev_2=punishBeta*(1-interestUM(iRevU,iRevM));
                    %�û��豸�˹��й���
                    tmpRev_3=alpha_1*sum(sum(powerSub_Pn));
                    %����=������-��Ȥ�ͷ�����-�豸���й���
                    temp_revenueUM(iRevU,iRevM)=tmpRev_1-tmpRev_2-tmpRev_3;
                end
            end
        end
                       
    end % if pushUM(maxRev_U,maxRev_M)==1
    
    % ��������û�k*֮�󣬶ಥ��m*��������
    revAfter=sum(revenueUM(:,maxRev_M));
    
    % ����������û����Ͷಥ��SINR�Ӷ����������ಥ��������½��Ļ�
    % ��ô����������û�
    if revAfter >= revPrevious
        pushUM(maxRev_U,maxRev_M)=1;
        %�����º��������Ϣд���������revenueUM
        revenueUM = temp_revenueUM;
    else
        pushUM(maxRev_U,maxRev_M)=0;
    end
    
end % if revenueUM(maxRev_U,maxRev_M)>0
revenueUM(maxRev_U,maxRev_M)=0;

% �жϵ����Ƿ����
 % ���revenueUM��ֻʣ��0�͸�����ô�������
    if max(max(revenueUM))>0;
        isIterationDone=0;
    else
        isIterationDone=1;
    end

end % while

% ���ɶಥ����Ȥ����omegaM��
[ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM, interestUM, alpha_2, xi );

% *************************************************************************
%                   ��������ṹ��
% *************************************************************************
OP1_Step1_Output.omegaM=omegaM;
OP1_Step1_Output.minSNR_gamma=minSNR_gamma;
OP1_Step1_Output.pushUM=pushUM;

end %Function