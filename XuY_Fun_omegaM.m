function [ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM,interestUM )
% XuY_Fun_omegaM: calculate omegaM
%   �ಥ����Ȥ����omegaM���ļ���

% Author	: XuY
% Modified	: 2014-12-11 15:00

% Create and Modify Date and History :
% -

% Error Case :
% -

% Call 
% -

% Called
% -function [ OP1_Step1_Output ] = XuY_Fun_StepOptimize1_Step1( OP1_Step1_Input )

%omegaM
omegaM=zeros(1,TOTAL_MULTIGROUP);

for iM=1:TOTAL_MULTIGROUP
    for iU=1:TOTAL_USER
        if pushUM(iU,iM)==1&&interestUM(iU,iM)>=0.5;%����Ȥ����ֵ����Ϊ0.5
            %������û����ڵ�ǰ�����Ķಥ��iMultiGroup���Ҹ��û��Զಥ��iMultiGroup����Ȥ
            omegaM(1,iM)=omegaM(1,iM)+1;%��iMultiGroup����Ȥֵ��1
        end
    end
end
omegaM=omegaM/100;

end

