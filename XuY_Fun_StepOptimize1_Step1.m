function [ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input )
%XuY_Fun_StepOptimize1_Step1: Step1 of the StepOptimize1
%   �������ز������빦�ʷ����£��ಥ�����û�����ѡ��

% Author	: XuY
% Modified	: 2014-12-11 15:00

% Create and Modify Date and History :
% -

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
t=OP1_Step1_Input.t;%��ʽ�滮�Ĳ�����������Ϊ0��������
gainChannel=OP1_Step1_Input.gainChannel;

%Լ������
MAX_POWER_Pth=OP1_Step1_Input.MAX_POWER_Pth;%����Լ�� ��λW

%-----------------�����ֲ�����----------------------
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);%�û��Զಥ����Ȥ����
punishBeta=10;%�ͷ�������ϵ��
powerUE=1*(1e-3);%�û��˽��չ��ģ���λ:W
powerBS=1*(1e-1);%��վ���ģ���λ:W

% --------------------------------------------
%       ��ʼ���ز����伯��scheduleSub_rho
% --------------------------------------------
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
for iScheduleSub=1:TOTAL_SUB
    %����������б�֤ÿ��ֻ��һ��1����ÿ�����ز�ֻ�ܱ����䵽һ���ಥ��
    index_rho=randperm(TOTAL_MULTIGROUP);
    scheduleSub_rho(iScheduleSub,index_rho(1))=1;
end

% --------------------------------------------
%               ��ʼ�����ز�����
% --------------------------------------------
powerSub_Pn=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
tempPn=rand(TOTAL_SUB,1);
%��һ��������ܹ��ʣ���֤ÿ�����ز��Ĺ���֮�͵����ܹ���
oneTempPn=MAX_POWER_Pth*tempPn./sum(tempPn(:));
% ��ÿ�����ز��Ĺ��ʶ�Ӧ��ÿ���ಥ����
for iPnN=1:TOTAL_SUB
    for iPnM=1:TOTAL_MULTIGROUP 
        if scheduleSub_rho(iPnN,iPnM)==1
            powerSub_Pn(iPnN,iPnM)=oneTempPn(iPnN,1);
        end
    end
end

% --------------------------------------------
%      ��ʼ��ÿ���ಥ�鼯���г�ʼ������SINR
% --------------------------------------------
% ��ʼ����ʱ����һ���ܴ��ֵ�������û��ڵ��������м��뼯��ʱ
% �ٸ��ݼ�����û���SINR�Լ��Ͻ��и��£�
minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);

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
        %�û��豸�˹��й���
        tmpRev_3=t*powerUE;
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
        %��������Ϊ�㣬��ֹ�ظ�ѡ��
        revenueUM(maxRev_U,maxRev_M)=0;
    else
      pushUM(maxRev_U,maxRev_M)=1;

        % ���¶ಥ�����û������������� SINR
        % minSNR_gamma=(1e+9).*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
            for iSNR_N=1:TOTAL_SUB
                if gainChannel(maxRev_U,iSNR_N,1)<minSNR_gamma(iSNR_N,maxRev_M)
                    minSNR_gamma(iSNR_N,maxRev_M)=gainChannel(maxRev_U,iSNR_N,1);
                end
            end

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
                    tmpRev_3=t*powerUE;
                    %����=������-��Ȥ�ͷ�����-�豸���й���
                    revenueUM(iRevU,iRevM)=tmpRev_1-tmpRev_2-tmpRev_3;
                end
            end
        end
        
    end % if pushUM(maxRev_U,maxRev_M)==1
    
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
[ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM,interestUM );

% ��������ṹ��
OP1_Step1_Output.omegaM=omegaM;
OP1_Step1_Output.minSNR_gamma=minSNR_gamma;
OP1_Step1_Output.pushUM=pushUM;

end %Function