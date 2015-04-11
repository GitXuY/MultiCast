load gainChannel
load test_step1_new
% *************************************************************************
%                       ��������
% *************************************************************************
%-------------------�������------------------------
% �ಥ������
TOTAL_MULTIGROUP=2;
% �û�����
TOTAL_USER=100;
% ���ز�����
TOTAL_SUB=64;
% ���ز�����һ������ȡ0.3125MHz
BAND_SUB= 0.3125;
% ��վ���Ķ�Ч�õ�Ȩ��
alpha_1=1;
% �û����Ķ�Ч�õ�Ȩ��
alpha_2=1;
% ��ʾ�û�û���յ�λ����ҵ��Ķ��⹦��
xi=0.001;
% �û���Ȥ�����û����ĸ��ಥ�����Ȥ
interestUM=rand(TOTAL_USER,TOTAL_MULTIGROUP);
% ����Լ�� ��λW
MAX_POWER_Pth=1;
% ��С��������Լ�� ��λ��
MIN_RATE_Rmin=20;
%47��ʱ��PTS�ܴﵽ1��RTS�ӽ�1����3000����������

% omegaM;

%-----------------�����ֲ�����----------------------
%�������ն�ż���ӳ�ʼ��
lambda=20;%��ż����lambda
mu=ones(1,TOTAL_MULTIGROUP);%��ż����mu
%-----------------��ʼ��----------------------
% ��ʼ�����ز�����
powerSub_Pn=MAX_POWER_Pth/TOTAL_SUB.*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
% ��ʼ���ز����伯��scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);

%��������
isScheduleDone=0;% ��Դ������ɱ�־
MAX_ITERATIONS=3000;%����������
timesIteration=1;%��ʼ������������
while isScheduleDone==0
% *************************************************************************
%                       �׶�1�����ز�����
% *************************************************************************
% ÿ�����ز��Ĺ�����ͬ���Ȱ����ز�������ಥҵ����
%--------------------------------------------------
%       ����ÿ�����ز��ڶ�Ӧ�ಥ��ʹ��ʱ������
%--------------------------------------------------
rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%ÿ�����ز��Ĵ�������
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
%       ̰���㷨�������ز�
%--------------------------------------------------
% �ҵ������������ز�����λ
while sum(scheduleSub_rho(:))~=TOTAL_SUB
    [maxRate, bestSub]=max(rateSub_greedy);
    % ������ز��ڶ���ಥ�������ʾ�Ϊ������ѡȡһ���ಥ��
    if length(bestSub) > 1
        bestMultigroup = ceil(length(bestSub)*rand(1));
    end
    scheduleSub_rho(bestSub(bestMultigroup),bestMultigroup)=1;
    rateSub_greedy(bestSub(bestMultigroup),:)=0;
end
% *************************************************************************
%                       �׶�2�����ʷ���
% *************************************************************************
for iM=1:TOTAL_MULTIGROUP
    for iN=1:TOTAL_SUB
        if scheduleSub_rho(iN,iM)==1
            powerSub_Pn(iN,iM)=(mu(1,iM)+omegaM(1,iM))...
                *BAND_SUB./((alpha_1+lambda)*log(2))-1./minSNR_gamma(iN,iM);
            % ����ȡ����
            powerSub_Pn(1,iM)=max(powerSub_Pn(1,iM),0);
        else
          powerSub_Pn(iN,iM)=0;
        end
    end
end
% *************************************************************************
%                   �׶�3����ż�������ݶȷ�����
% *************************************************************************
    % ====================================================
    %     ��1����һ�ε���ʱʹ�ó�ʼ��ֵ��������������
    % ====================================================
%     if timesIteration==1
        %Lambda�ĸ��²�����ʼ��
        stepLambda_c1= 50*(1e-1)/sqrt(timesIteration);
        
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
                +sum( abs( (mu2(1,iDual)-mu(1,iDual))./( mu(1,iDual)+1*(1e-5) ) ) );
%         end
    end
    totalDualChange
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
        
%     % �����С���ʺ������Լ�������㣬��С����
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
    lambda = lambda2
    mu = mu2
    timesIteration=timesIteration+1;
    isScheduleDone=isExceedIteration+isConverge;
    PTS=sum(powerSub_Pn(:))/MAX_POWER_Pth
    RTS=rateM/MIN_RATE_Rmin
    
end % while ��ѭ������

