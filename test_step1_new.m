clc
clear
load gainChannel
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
MIN_RATE_Rmin=1.5;
% �ɶಥ���ͼ������û���Ը����ո�ҵ����ɵ�������ʧ
betaM=10*ones(1,TOTAL_MULTIGROUP);
% �������
revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% �ݴ��������
temp_revenueUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
% ÿ���û��䵱��׼�û�ʱ����ѿ�ѡ�ಥ�����û�����
pushUUM=zeros(TOTAL_USER,TOTAL_USER,TOTAL_MULTIGROUP); 
% powerUE=1*(1e-3);%�û��˽��չ��ģ���λ:W
% powerBS=1*(1e-1);%��վ���ģ���λ:W

%-------------------step2����------------------------
% ��ʼ���ز����伯��scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
for iScheduleSub=1:TOTAL_SUB
    %����������б�֤ÿ��ֻ��һ��1����ÿ�����ز�ֻ�ܱ����䵽һ���ಥ��
    index_rho=randperm(TOTAL_MULTIGROUP);
    scheduleSub_rho(iScheduleSub,index_rho(1))=1;
end

% ��ʼ�����ز�����
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

%-------------------�������------------------------
% omegaM;
% ��׼�û�����
criterionUserM=zeros(1,TOTAL_MULTIGROUP);
% �������
maxRev=zeros(1,TOTAL_MULTIGROUP);
% ÿ���ಥ�鼯����������SINR
minSNR_gamma=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
% ҵ��ಥ���ͼ���
pushUM=zeros(TOTAL_USER,TOTAL_MULTIGROUP);
%******************************************************
%                �����û�ѡ��
%******************************************************
LOOP=3;
m=0;
for iM=1:TOTAL_MULTIGROUP
    %�����ಥҵ����
    for iU=1:TOTAL_USER
        %���������û�Ⱥ�壬ÿ���û��䵱һ�λ�׼�û�
        criterionUser=iU;%���û�׼�û�
        isAdd=ones(1,TOTAL_USER); %����û��Ƿ��ܼ��뼯��
        % ��׼�û���ƽ��SNR
        iU_average_SNR = 0;
        for iN=1:TOTAL_SUB
            if scheduleSub_rho(iN,iM)==1
            %�������䵽�öಥ������ز�
                iU_average_SNR = iU_average_SNR + gainChannel(criterionUser,iN,LOOP);
            end
        end
        iU_average_SNR = iU_average_SNR / TOTAL_SUB;          
        for i2U=1:TOTAL_USER
            %��ÿ����׼�û���ѡ����ѿ�ѡ�ಥ�����û�����
            %--------------------------------------------------
            %                   i2Uƽ��SNR
            %--------------------------------------------------
            i2U_average_SNR = 0;
            for iN=1:TOTAL_SUB
                if scheduleSub_rho(iN,iM)==1
                %�������䵽�öಥ������ز�
                    i2U_average_SNR = i2U_average_SNR + gainChannel(i2U,iN,LOOP);
                end
            end
            i2U_average_SNR = i2U_average_SNR / TOTAL_SUB;  
            %--------------------------------------------------
            %     i2Uƽ��SNR�Ȼ�׼�û���SNRС�Ĳ��ܼ��뼯��
            %--------------------------------------------------
            if i2U_average_SNR<iU_average_SNR
                isAdd(1,i2U)=0;
            end
            %--------------------------------------------------
            %           ����С��0�Ĳ��ܼ��뼯��
            %--------------------------------------------------
            tmp1=0;
            for i2N=1:TOTAL_SUB
                %�ಥ�������
                tmp1=tmp1+BAND_SUB*scheduleSub_rho(i2N,iM)...
                    *log2(1+powerSub_Pn(i2N,iM)...
                    *gainChannel(criterionUser,i2N,LOOP));
            end
            %�û��Խ���ҵ�񲻸���Ȥ�ĳͷ�����
            tmp2=interestUM(i2U,iM)-alpha_2*xi;
            %�û��豸�˹��й���
            tmp3=betaM(1,iM)*(1-interestUM(i2U,iM));
            %����=������-��Ȥ�ͷ�����-�豸���й���
            temp_revenueUM(i2U,iM)=tmp1*tmp2-tmp3;
            if temp_revenueUM(i2U,iM) < 0
                isAdd(1,i2U)=0;
            end    
        end % for i2U=1:TOTAL_USER
        %--------------------------------------------------
        %  �������û�����¶ಥ�����û������������� SINR
        %--------------------------------------------------
        %������������û�
        if isAdd(1,i2U)==1
            % �������ز�
            for iN=1:TOTAL_SUB
               % �����¸öಥ����䵽�����ز�
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
        %   iU�û��䵱��׼�û�ʱ,��ѿ�ѡ�ಥ�����û�����
        %--------------------------------------------------
        pushUUM(iU,:,iM)=isAdd(:);
        %--------------------------------------------------
        % iU�û��䵱��׼�û�ʱ,��ѿ�ѡ�ಥ�����û����ϵ�����
        %--------------------------------------------------
        for i3U=1:TOTAL_USER
            tmp1_1=0;
            if isAdd(1,i3U)==1
                for i3N=1:TOTAL_SUB
                %�ಥ�������
                tmp1_1=tmp1_1+BAND_SUB*scheduleSub_rho(i3N,iM)...
                    *log2(1+powerSub_Pn(i3N,iM)...
                    *gainChannel(criterionUser,i3N,LOOP));
                end
                %�û��Խ���ҵ�񲻸���Ȥ�ĳͷ�����
                tmp1_2=interestUM(i3U,iM)-alpha_2*xi;
                %�û��豸�˹��й���
                tmp1_3=betaM(1,iM)*(1-interestUM(i3U,iM));
                %����=������-��Ȥ�ͷ�����-�豸���й���
                revenueUM(iU,iM)=revenueUM(iU,iM)+tmp1_1*tmp1_2-tmp1_3;
            end
        end
        
    end % for iU=1:TOTAL_USER
    
    [maxRev(1,iM), criterionUserM(1,iM)] = max(revenueUM(:,iM));
    pushUM(:,iM)=pushUUM(criterionUserM(1,iM),:,iM);
    
end
% ���ɶಥ����Ȥ����omegaM��
[ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM, interestUM, alpha_2, xi );

save test_step1_new
