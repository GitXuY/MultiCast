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
MIN_RATE_Rmin=1.5;

% omegaM;

%-----------------�����ֲ�����----------------------
%�������ն�ż���ӳ�ʼ��
lambda=2;%��ż����lambda
mu=ones(1,TOTAL_MULTIGROUP);%��ż����mu

%��������
isScheduleDone=0;% ��Դ������ɱ�־
MAX_ITERATIONS=10000;%����������
timesIteration=1;%��ʼ������������
% *************************************************************************
%                       �׶�1�����ز�����
% *************************************************************************
% ÿ�����ز��Ĺ�����ͬ���Ȱ����ز�������ಥҵ����

% ��ʼ�����ز�����
powerSub_Pn=MAX_POWER_Pth/TOTAL_SUB.*ones(TOTAL_SUB,TOTAL_MULTIGROUP);
% ��ʼ���ز����伯��scheduleSub_rho
scheduleSub_rho=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
% ÿ���ಥ�鼯����������SINR
minSNR_gamma=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);
for iM=1:TOTAL_MULTIGROUP
    minSNR_gamma(:,iM)=gainChannel(criterionUserM(1,iM),:,LOOP); 
end
% ����ÿ�����ز��ڶ�Ӧ�ಥ��ʹ��ʱ������
rateSub=zeros(TOTAL_SUB,TOTAL_MULTIGROUP);%ÿ�����ز��Ĵ�������
for iM=1:TOTAL_MULTIGROUP
    for iN=1:TOTAL_SUB
        rateSub(iN,iM)...
             =BAND_SUB*log2(1+powerSub_Pn(iN,iM)...
                  *minSNR_gamma(iN,iM));
    end
end




