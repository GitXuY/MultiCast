function [ Min_Output ] = XuY_Fun_UpperLimit_Min( Min_Input )
%XuY_Fun_UpperLimit_Min: minimize the maximal objective function
%   Ŀ�꺯���ⲿmin����

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
% ��������
TOTAL_MULTIGROUP=XuY_Fun_UpperLimit_Min.TOTAL_MULTIGROUP;%�ಥ������
TOTAL_USER=XuY_Fun_UpperLimit_Min.TOTAL_USER;%�û�����
TOTAL_SUB=XuY_Fun_UpperLimit_Min.TOTAL_SUB;%���ز�����
BAND_SUB=XuY_Fun_UpperLimit_Min.BAND_SUB;%���ز�����һ������ȡ0.3125MHz

%-----------------�����ֲ�����----------------------
%�������ն�ż���ӳ�ʼ��
lambda=2;%��ż����lambda
mu=ones(1,TOTAL_MULTIGROUP);%��ż����mu

%��������
isScheduleDone=0;% ��Դ������ɱ�־
MAX_ITERATIONS=10000;%����������
timesIteration=1;%��ʼ������������

%������Դ����Ĵ�ѭ��
while isScheduleDone==0
    
% *************************************************************************
%                       �׶�1��Pn��D(n,m)�ļ���
% *************************************************************************
%   �洢��Ϣ�����ʼ��
    funcD=zeros(1,TOTAL_MULTIGROUP);%D(n,m) 1*�ಥ����
    rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%ÿ�����ز��Ĵ�������
    powerSub_Pn=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%ÿ�����ز��Ĺ���
    powerSub_Pn_Cal=zeros(1,TOTAL_MULTIGROUP);%����ʱ�洢Pnֵ����ʱ����
      
    rateM=zeros(1,TOTAL_MULTIGROUP);%ÿ���ಥ��ĺ�����
    powerM=zeros(1,TOTAL_MULTIGROUP);%ÿ���ಥ��ĺ͹���
    
    isSubSchedule=zeros(1,TOTAL_SUB);%��¼���ز��Ƿ��Ѿ��������ȥ
    scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%���ز��ķ��亯��
    
%   ��������    
    for iN=1:TOTAL_SUB
        for iiM=1:TOTAL_MULTIGROUP
            % ��������ز���δ������
            if isSubSchedule(iN)==0
                % ���ز�iN�����䵽��iiM���ಥ��ʱ�����Ż�����
                powerSub_Pn_Cal(1,iiM)=(mu(1,iiM)+omegaM(1,iiM))...
                    *BAND_SUB./((alpha_1+lambda)*log(2))-1./minSNR_gamma(iN,iiM);
                
                % ����ȡ����
                powerSub_Pn_Cal(1,iiM)=max(powerSub_Pn_Cal(1,iiM),0);
                
                % ��Ӧ�ù���ֵ�����õ���D(n,m)ֵ
                funcD(1,iiM)=(mu(1,iiM)+omegaM(1,iiM))*BAND_SUB*...
                    log2(1+powerSub_Pn_Cal(1,iiM)*minSNR_gamma(iN,iiM))-...
                    (alpha_1+lambda)*powerSub_Pn_Cal(1,iiM);
            % ��������ز��Ѿ�������
            else
                % �������ز�����ʱ����������
                powerSub_Pn_Cal(1,iiM)=0;
                % funcD��-1000ȷ�������ز�����ѡ��
                funcD(1,iiM)=-1000;
            end % if isSubSchedule
        end % for iiM
        
        % �Ե�n�����ز�����������ز������䵽��һ���ಥ���ʱ��
        % ���õ���Ч�����D(n,m)ȡ�����
        bestMultiGroupForiN = find( funcD>=max(funcD) );
        
        % ����ҵ����������ϵ�ms,����������������ô���ѡһ��
        if length(bestMultiGroupForiN) > 1
            bestMultiGroupForiN...
                            = round(length(bestMultiGroupForiN*rand(1)));
        end
        
        % �Ѹ����ز����������ಥ��
        scheduleSub_rho(iN,bestMultiGroupForiN) = 1;
        
        % ���������ز��Ѿ��������ȥ����ֹ�ظ�����
        isSubSchedule(1,iN)=1;
        
        % �Ѹ����ز��ڡ���ʱ�����ز����ʾ����е�ֵ�Ǽǵ����ز����ʾ���
        powerSub_Pn(iN,bestMultiGroupForiN)...
                                = powerSub_Pn_Cal(1,bestMultiGroupForiN);
        
        % ��������ز�������
        rateSub(iN,bestMultiGroupForiN)...
                    =BAND_SUB*log2(1+powerSub_Pn(iN,bestMultiGroupForiN)...
                            *minSNR_gamma(iN,bestMultiGroupForiN));
        
    end % for iN
    
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
    if timesIteration==1
        %Lambda�ĸ��²�����ʼ��
        stepLambda_c1= 5*(1e-1);
        
        %Mu�ĸ��²�����ʼ��
        stepMU_c2 = 1*(1e-3);
    end
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
    DualChangeThreshold = 5*(1e-5);%��ż�����仯����֮�͵���ֵ
    
    for iDual=1:TOTAL_MULTIGROUP
        if mu(1,iDual)>0
            %�����ż����ÿ�α仯����֮��
            totalDualChange = totalDualChange+sum(abs( (lambda2-lambda)./(lambda)))...
                +sum(abs( (mu2(1,iDual)-mu(1,iDual))./(mu(1,iDual))));
        end
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
    if isBelowMAX_POWER==1&&isAboveMIN_RATE==1
        stepLambda_c1 = 1*(1e-1)/sqrt(timesIteration);%Lambda�ĸ��²���
        stepMU_c2 = 1*(1e-4)/sqrt(timesIteration);%Mu�ĸ��²���
        
    % ���������������3000����һ����С����
    %   ע����ʱһ��ӽ�������������ֱ�ӰѲ������ø�С
    elseif timesIteration>3000
        stepLambda_c1 = 1*(1e-2)/sqrt(timesIteration);%Lambda�ĸ��²���
        stepMU_c2 = 1*(1e-5)/sqrt(timesIteration);%Mu�ĸ��²���
        
    % �������Լ���������㣬�������㲽��
    else
        stepLambda_c1 = 5*(1e-1)/sqrt(timesIteration);%Lambda�ĸ��²���
        stepMU_c2 = 1*(1e-3)/sqrt(timesIteration);%Mu�ĸ��²���
    end
    % ====================================================
    %       ��5����������
    % ====================================================
    lambda = lambda2;
    mu = mu2;
    timesIteration=timesIteration+1;
    isScheduleDone=isExceedIteration+isConverge;
    PTS=sum(powerSub_Pn(:))/MAX_POWER_Pth;
    RTS=rateM/MIN_RATE_Rmin;
    
end % while ��ѭ������

% *************************************************************************
%                   ��������ṹ��
% *************************************************************************
% ��Դ������
OP1_Step2_Output.powerSub_Pn=powerSub_Pn;
OP1_Step2_Output.rateSub=rateSub;
OP1_Step2_Output.rateM=rateM;
OP1_Step2_Output.powerM=powerM;
OP1_Step2_Output.scheduleSub_rho=scheduleSub_rho;
OP1_Step2_Output.PTS=PTS;
OP1_Step2_Output.RTS=RTS;

% �������ն�ż���ӽ��
OP1_Step2_Output.lambda=lambda;
OP1_Step2_Output.mu=mu;

% ������Ϣ
OP1_Step2_Output.timesIteration=timesIteration;
OP1_Step2_Output.isExceedIteration=isExceedIteration;
OP1_Step2_Output.isAboveMIN_RATE=isAboveMIN_RATE;
OP1_Step2_Output.isBelowMAX_POWER=isBelowMAX_POWER;
OP1_Step2_Output.totalDualChange=totalDualChange;
OP1_Step2_Output.isBelowDualChange=isBelowDualChange;

end %Function
