function [ Max_Output ] = XuY_Fun_UpperLimit_Max( Max_Input )
%XuY_Fun_UpperLimit_Max: maximize objective function
%   Ŀ�꺯���ڲ�max����

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
%                       ��������
% *************************************************************************
%-------------------�������------------------------
%   ��������
TOTAL_MULTIGROUP    = Max_Input.TOTAL_MULTIGROUP;%�ಥ������
TOTAL_USER          = Max_Input.TOTAL_USER;%�û�����
TOTAL_SUB           = Max_Input.TOTAL_SUB;%���ز�����
BAND_SUB            = Max_Input.BAND_SUB;%���ز�����һ������ȡ0.3125MHz
gainChannel         = Max_Input.gainChannel;
interestUM          = Max_Input.interestUM;% rand(TOTAL_USER,TOTAL_MULTIGROUP);
alpha_1             = Max_Input.alpha_1;
alpha_2             = Max_Input.alpha_2;
xi                  = Max_Input.xi;
betaM               = Max_Input.betaM;
%   ��������
lambda              = Max_Input.lamda;
mu                  = Max_Input.mu;

%-------------------���������------------------------
% ���ݻ�׼�û��Ĳ��������ز��Ĺ���
powerSub_Pn = zeros(1,TOTAL_SUB);

% ���ݻ�׼�û��Ĳ������ܱ����䵽�ಥ����û�����
userSet_Knm  = zeros(TOTAL_USER,TOTAL_MULTIGROUP,TOTAL_SUB);

% D(n,m) 1*�ಥ����
funcD=0;%D(n,m) 1*�ಥ����
funcD_plus=0;

% �Ż�����
betaMN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);

% ��¼�ಥ�������
revenueMN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);

% ��¼�ಥ��Ļ�׼�û�
criterionUSER_MN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);

% ��¼�ಥ��Ļ�׼SNR
criterionSNR_MN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);

% ��¼ÿ���û��䵱��׼�û�ʱ �ಥ�������
revenueUMN = zeros(TOTAL_USER,TOTAL_MULTIGROUP, TOTAL_SUB);

% ��¼ÿ���û��䵱��׼�û�ʱ �ಥ��Ļ�׼SNR
criterionSNR_UMN = zeros(TOTAL_MULTIGROUP, TOTAL_SUB);

% ��¼ÿ�����ز��ɻ�õ��������
maxRevenue = zeros(1,TOTAL_SUB);

% ��¼ÿ�����ز���Ӧ�����Ŷಥ��
bestMG = zeros(1,TOTAL_SUB);



% *************************************************************************
%                       �׶�1��
% *************************************************************************
for iSub=1:TOTAL_SUB
% �ڲ� max ������Բ��ΪN �������⣬���е� n �������⸺�����ز�n ��ҵ��ѡ
% �񡢹����Լ���Ӧ�Ķಥ�����û�ѡ��

    for iMG=1:TOTAL_MULTIGROUP
    % �������ز�n �����ҵ��m �������ʱ��ѵĹ����������û�����,���� M �ֿ��ܵ�
    % ���ز�������ԣ�ѡ����ѵ����ز�����
    betaMN = betaMN + TOTAL_MULTIGROUP * betaM(iMG) / TOTAL_SUB;

        for iUserCSNR = 1:TOTAL_USER
        % ��˳������ѡ���û���Ϊ��׼�û���Ȼ��ѡ����ʹЧ�������û���Ϊ��ѻ�׼�û���
        
            userSet_Knm(:,iMG,iSub) = 0;    
            criterionSNR_UMN(iUserCSNR, iMG, iSub) = gainChannel(iUserCSNR,iSub,2);   %����������ΪLoop
            % =============================================================
            %   ����׼�û�Ϊ iUserCSNR ʱ �����û�����userSet_Knm
            % =============================================================
            for iUser = 1: TOTAL_USER  
                % ���ݻ�׼�û�ѡ��iMG�Ĵ��伯��
                if gainChannel(iUser,iSub,1) >= criterionSNR_UMN(iUserCSNR, iMG, iSub) && interestUM(iUser,iMG)-alpha_2*xi+mu >= 0            
                    userSet_Knm(iUser,iMG,iSub) = 1;
                end  
            end 
            % =============================================================
            %     ����׼�û�Ϊ iUserCSNR ʱ ���ز�iSub����powerSub_Pn
            % =============================================================
            temp_p1 = interestUM-alpha_2*xi+mu;
            temp_p2 = userSet_Knm(:,:,iSub).* temp_p1;
            powerSub_Pn(iSub)=(BAND_SUB * sum(temp_p2(:)))...
            /(log(2)*(alpha_1+lambda))-1/criterionSNR_UMN(iUserCSNR, iMG, iSub);                                  
            % =============================================================
            %       ����׼�û�Ϊ iUserCSNR ʱ �ಥ��Ĵ�������
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
        
        %ѡ��ʹ�öಥ���������Ļ�׼�û� ���� ��Ӧ�Ķಥ������
        [revenueMN(iMG, iSub), criterionUSER_MN(iMG, iSub)]= max(revenueUMN(:,iMG, iSub));
        criterionSNR_MN(iMG, iSub) = criterionSNR_UMN(criterionUSER_MN(iMG, iSub),iMG, iSub);
        
    end % for iMG=1:TOTAL_MULTIGROUP
    
    [ maxRevenue(iSub), bestMG(iSub) ]= max(revenueMN(:, iSub));
    
end % for iSub=1:TOTAL_SUB
% *************************************************************************
%                   ��������ṹ��
% *************************************************************************
Max_Output.powerSub_Pn      = powerSub_Pn;
Max_Output.userSet_Knm      = userSet_Knm;
Max_Output.funcD            = funcD;%D(n,m) 1*�ಥ����
Max_Output.betaMN           = betaMN;
Max_Output.revenueMN        = revenueMN;
Max_Output.criterionUSER_MN = criterionUSER_MN;
Max_Output.revenueUMN       = revenueUMN;
Max_Output.maxRevenue       = maxRevenue;
Max_Output.bestMG           = bestMG;
Max_Output.criterionSNR_MN  = criterionSNR_MN;
end % Function

