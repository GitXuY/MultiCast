function [ omegaM ] = XuY_Fun_omegaM( TOTAL_MULTIGROUP, TOTAL_USER, pushUM,interestUM, alpha_2, xi)
% XuY_Fun_omegaM: calculate omegaM
%   多播组兴趣矩阵（omegaM）的计算

% Author	: XuY
% Modified	: 2014-12-11 15:00

% Create and Modify Date and History :
% - 2015/4/3 
%     - change the formula of omegaM 


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
        if pushUM(iU,iM)==1;
            %如果该用户属于当前遍历的多播组iMultiGroup，且该用户对多播组iMultiGroup感兴趣
            omegaM(1,iM)=omegaM(1,iM)+interestUM(iU,iM)-alpha_2*xi;%则iMultiGroup的兴趣值加1
        end
    end
end
omegaM=omegaM/100;      

end

