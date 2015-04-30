clc
clear
load gainChannel
% *************************************************************************
%                       ��������
% *************************************************************************
TOTAL_MULTIGROUP=2;%�ಥ������
TOTAL_USER=100;%�û�����
TOTAL_SUB=64;%���ز�����
BAND_SUB= 0.3125; %���ز�����һ������ȡ0.3125MHz
MAX_POWER_Pth=1;%����Լ�� ��λW
MIN_RATE_Rmin=30;% ��С��������Լ�� ��λ��

alpha_1=1;%��վ���Ķ�Ч�õ�Ȩ��
alpha_2=1;%�û����Ķ�Ч�õ�Ȩ��
xi=0.001;%��ʾ�û�û���յ�λ����ҵ��Ķ��⹦��
betaM=ones(1,TOTAL_MULTIGROUP);% �ɶಥ���ͼ������û���Ը����ո�ҵ����ɵ�������ʧ
% LOOP=10;%�ŵ�ѭ����

%�û���Ȥ�����û����ĸ��ಥ�����Ȥ
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);
lambda=2;
mu=ones(1,TOTAL_MULTIGROUP);%��ż����mu

%-------------------���������------------------------
% ���ݻ�׼�û��Ĳ��������ز��Ĺ���
powerSub_Pn = zeros(1,TOTAL_SUB);
% ���ݻ�׼�û��Ĳ������ܱ����䵽�ಥ����û�����
userSet_Knm  = zeros(TOTAL_USER,TOTAL_MULTIGROUP,TOTAL_SUB);
% D(n,m) 1*�ಥ����
funcD=0;%D(n,m) 1*�ಥ����
funcD_plus=0;
% �Ż�����
tmp_betaM = TOTAL_MULTIGROUP .* betaM / TOTAL_SUB;
betaMN = repmat(tmp_betaM',1,TOTAL_SUB);
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
% ��¼���ز�����
scheduleSub_rho = zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
%ÿ�����ز��Ĵ�������
rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
%ÿ���ಥ��ĺ�����
rateM=zeros(1,TOTAL_MULTIGROUP);
%ÿ���ಥ��ĺ͹���
powerM=zeros(1,TOTAL_MULTIGROUP);

%-------------------��������------------------------
isScheduleDone=0;% ��Դ������ɱ�־
MAX_ITERATIONS=3000;%����������
timesIteration=1;%��ʼ������������
while isScheduleDone~=1
%������Դ����Ĵ�ѭ��
% *************************************************************************
%     �׶�1���������ز�n �����ҵ��m �������ʱ��ѵĹ����������û�����
% *************************************************************************
for iSub=1:TOTAL_SUB
% �ڲ� max ������Բ��ΪN �������⣬���е� n �������⸺�����ز�n ��ҵ��ѡ
% �񡢹����Լ���Ӧ�Ķಥ�����û�ѡ��
if sum(scheduleSub_rho(iSub,:))==0

    for iMG=1:TOTAL_MULTIGROUP
    % �������ز�n �����ҵ��m �������ʱ��ѵĹ����������û�����,���� M �ֿ��ܵ�
    % ���ز�������ԣ�ѡ����ѵ����ز�����
%      betaMN (1) = betaMN + TOTAL_MULTIGROUP * betaM(iMG) / TOTAL_SUB;

        for iUserCSNR = 1:TOTAL_USER
        % ��˳������ѡ���û���Ϊ��׼�û���Ȼ��ѡ����ʹЧ�������û���Ϊ��ѻ�׼�û���
        
            userSet_Knm(:,iMG,iSub) = 0;    
            criterionSNR_UMN(iUserCSNR, iMG, iSub) = gainChannel(iUserCSNR,iSub,2);   %����������ΪLoop
            % =============================================================
            %   ����׼�û�Ϊ iUserCSNR ʱ �����û�����userSet_Knm
            % =============================================================
            for iUser = 1: TOTAL_USER  
                % ���ݻ�׼�û�ѡ��iMG�Ĵ��伯��
                if gainChannel(iUser,iSub,1) >= criterionSNR_UMN(iUserCSNR, iMG, iSub)... 
                     && interestUM(iUser,iMG)-alpha_2*xi+mu(1,iMG) >= 0            
                    userSet_Knm(iUser,iMG,iSub) = 1;
                end  
            end 
            % =============================================================
            %     ����׼�û�Ϊ iUserCSNR ʱ ���ز�iSub����powerSub_Pn
            % =============================================================
            temp_p1 = interestUM(:,iMG)-alpha_2*xi+mu(1,iMG);
            temp_p2 = userSet_Knm(:,iMG,iSub).* temp_p1;
            temp_p2 = temp_p2 / 300; 
            powerSub_Pn(iSub)=(BAND_SUB * sum(temp_p2(:)))...
            /(log(2)*(alpha_1+lambda))-1/criterionSNR_UMN(iUserCSNR, iMG, iSub);
            powerSub_Pn(iSub)=max(powerSub_Pn(iSub),0);
            % =============================================================
            %       ����׼�û�Ϊ iUserCSNR ʱ �ಥ��Ĵ�������
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
        
        %ѡ��ʹ�öಥ���������Ļ�׼�û� ���� ��Ӧ�Ķಥ������
        [revenueMN(iMG, iSub), criterionUSER_MN(iMG, iSub)]= max(revenueUMN(:,iMG, iSub));
        criterionSNR_MN(iMG, iSub) = criterionSNR_UMN(criterionUSER_MN(iMG, iSub),iMG, iSub);
        
    end % for iMG=1:TOTAL_MULTIGROUP
    
    [ maxRevenue(iSub), bestMG(iSub) ]= max(revenueMN(:, iSub));
    
    scheduleSub_rho( iSub, bestMG(iSub)) = 1;
    
    % ��������ز�������
    rateSub(iSub,bestMG(iSub))...
                    =BAND_SUB*log2(1+powerSub_Pn(iSub)...
                            *criterionSNR_MN(bestMG(iSub), iSub));
    
end
end % for iSub=1:TOTAL_SUB
        
    %����ÿ���ಥ��ĺ�����                    
    for irateM=1:TOTAL_MULTIGROUP
        rateM(1,irateM)=sum(rateSub(:,irateM));
    end 
        
    %����ÿ���ಥ��ĺ͹���
    for ipowerM=1:TOTAL_MULTIGROUP
        powerM(1,ipowerM)=sum(powerSub_Pn(:,ipowerM));
    end
% *************************************************************************
%                   �׶�2����ż�������ݶȷ�����
% *************************************************************************
    % ====================================================
    %     ��1����һ�ε���ʱʹ�ó�ʼ��ֵ��������������
    % ====================================================
%     if timesIteration==1
        %Lambda�ĸ��²�����ʼ��
        stepLambda_c1= 5*(1e-1)/sqrt(timesIteration);
        
        %Mu�ĸ��²�����ʼ��
        stepMU_c2 = 1*(1e-3)/sqrt(timesIteration);
%     end
    % ====================================================
    %      ��2��lambda��mu�ĸ���
    % ====================================================
    %����lambda
    lambda2 = lambda + stepLambda_c1*( sum(powerSub_Pn(:))-MAX_POWER_Pth);
    lambda2 = max(0,lambda2);
    
    %����mu
    mu2=zeros(1,TOTAL_MULTIGROUP);
    for i4M=1:TOTAL_MULTIGROUP
        mu2(1,i4M) = mu(1,i4M)- stepMU_c2*(rateM(1,i4M)-MIN_RATE_Rmin);
        mu2(1,i4M) = max(0,mu2(1,i4M));
    end
    
    % ====================================================
    %        ��3���������жϴ��ݶȷ��Ƿ�ﵽ��������
    % ====================================================
    % �жϴ��ݶȷ��Ƿ����������������жϣ�
    % ����1��
        % [��ż�����仯 < ��ֵ] && [�͹��� < �����] && [���� > ��С����]
    % ����2��
        % ������������������ޣ���ֹ���ݶȷ�����
        
    % ����1��
    %�ж� [��ż�����仯 < ��ֵ]
    totalDualChange=0;%��¼��ż����ÿ�α仯����֮��
    DualChangeThreshold = 1*(1e-3);%��ż�����仯����֮�͵���ֵ
    
    for iDual=1:TOTAL_MULTIGROUP
%         if mu(1,iDual)>0
            %�����ż����ÿ�α仯����֮��
            totalDualChange = totalDualChange+sum(abs( (lambda2-lambda)./(lambda+1*(1e-5))))...
                +sum(abs( (mu2(1,iDual)-mu(1,iDual))./(mu(1,iDual)+1*(1e-5))));
%         end
    end
    
    if totalDualChange < DualChangeThreshold
        isBelowDualChange=1;
    else isBelowDualChange=0;
    end
    
    %�ж� [�͹��� < �����]
    if sum(powerSub_Pn(:))<MAX_POWER_Pth
        isBelowMAX_POWER=1;
    else isBelowMAX_POWER=0;
    end
    
    %�ж� [���� > ��С����]
    isAboveMIN_RATE=1;
    isAboveMIN_RATE_M=zeros(1,TOTAL_MULTIGROUP);
    for irateM=1:TOTAL_MULTIGROUP
        if rateM(1,irateM)>MIN_RATE_Rmin;
            isAboveMIN_RATE_M(1,irateM)=1;
        else isAboveMIN_RATE_M(1,irateM)=0;
        end
        isAboveMIN_RATE=isAboveMIN_RATE*isAboveMIN_RATE_M(1,irateM);
    end
    
    %����������־
    isConverge=isBelowDualChange*isBelowMAX_POWER*isAboveMIN_RATE;% ������־
    
    % ����2��
    isExceedIteration=0;%�ﵽ�������ޱ�־
    if timesIteration>=MAX_ITERATIONS
        %���ڳ�������������������ѭ��
        isExceedIteration=1;
    end
    
    % ====================================================
    %       ��4�������������жϽ�������¼��㲽��
    % ====================================================
    % ���²������ŵ������������Ӷ����٣�ʹ�ö�ż�����ı仯���Ӿ�ϸ������������
        
    % �����С���ʺ������Լ�������㣬��С����
%     if isBelowMAX_POWER==1&&isAboveMIN_RATE==1
%         stepLambda_c1 = 1*(1e-1)/sqrt(timesIteration);%Lambda�ĸ��²���
%         stepMU_c2 = 1*(1e-4)/sqrt(timesIteration);%Mu�ĸ��²���
%         
%     % ���������������3000����һ����С����
%     %   ע����ʱһ��ӽ�������������ֱ�ӰѲ������ø�С
%     elseif timesIteration>3000
%         stepLambda_c1 = 1*(1e-2)/sqrt(timesIteration);%Lambda�ĸ��²���
%         stepMU_c2 = 1*(1e-5)/sqrt(timesIteration);%Mu�ĸ��²���
%         
%     % �������Լ���������㣬�������㲽��
%     else
%         stepLambda_c1 = 5*(1e-1)/sqrt(timesIteration);%Lambda�ĸ��²���
%         stepMU_c2 = 1*(1e-3)/sqrt(timesIteration);%Mu�ĸ��²���
%     end
    % ====================================================
    %       ��5����������
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