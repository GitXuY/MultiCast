function [gainChannel]=XuY_Fun_ChannelGeneration(~)
% XuY_Fun_ChannelGeneration: Generate Channel Model
%   -产生信道增益矩阵
%   -用户撒点: 用户在半径为小区半径的圆上均匀分布
%   -信道生成: 信道考虑小尺度模型（瑞利信道），大尺度模型（自由空间路损和阴影效应）

% Author	: XuY
% Modified	: 2014-11-29 12:44

% Create and Modify Date and History :
% -2015/1/5 
%       把瑞利信道放在了循环里面，保证每次loop时，大小尺度衰落都是重新生成的

% Error Case :
% -

% Call 
% - function [lossRayleigh]=XuY_Fun_lossRayleigh(TOTAL_SUB,h)

% Called
% -

clc
clear all

% 场景说明
% 假设小区圆形，半径为500 m，基站位于圆心；
% 多波用户位于一个小圆内，圆心与基站间的距离在[ 100m, 400m ]内均匀分布

%******************************************************
%                参量设置
%******************************************************
N0= 10^(-14.4); %噪声功率谱密度，单位W/MHz
RAD_CELL = 500;%小区半径，单位:m
TOTAL_USER = 100;%最大用户数
LOOP = 10; %信道实现次数
j = sqrt(-1);%虚数单位
TOTAL_SUB=64;%子载波总数
BAND_SUB= 0.3125; %子载波带宽，一般现在取0.3125MHz

lossLargeScare=zeros(1,TOTAL_USER,LOOP);%大尺度信道衰减矩阵
lossSmallScare=zeros(1,TOTAL_SUB,LOOP);%小尺度信道衰减矩阵
gainChannel = zeros( TOTAL_USER, TOTAL_SUB, LOOP );%每个用户与基站间的信道增益

for iLoop = 1:LOOP
    %******************************************************
    %                用户撒点
    %******************************************************
    %用户在半径为小区半径的圆上均匀分布
    locUser = ones(1,TOTAL_USER);%用户位置 
    iUser=1;
    while iUser<=TOTAL_USER
        xUser = RAD_CELL*(2*rand-1);%用户在正方形内均匀分布
        yUser = RAD_CELL*(2*rand-1);
        if (xUser^2+yUser^2) <= RAD_CELL^2 && (xUser^2+yUser^2)>100^2
            %剔除在正方形内但在圆外的用户 && 剔除离基站太近的用户
             locUser(1,iUser) = xUser + j*yUser;
             iUser=iUser+1;
        end
    end

    %******************************************************
    %                大尺度信道生成
    %******************************************************
    % 信道考虑小尺度模型（瑞利信道），大尺度模型（自由空间路损和阴影效应）
    % lossLargeScare=zeros(1,TOTAL_USER,LOOP);%大尺度信道衰减矩阵

    lossPath = 17.4 + 37.6.*log10( abs(locUser) );%路径损耗矩阵
    lossShadow = 8.*randn(1,TOTAL_USER);%阴影效应矩阵
    lossLargeScare(:,:,iLoop)=lossPath + lossShadow;
    lossLargeScare(:,:,iLoop) = 10.^(-lossLargeScare(:,:,iLoop)./10);%将对数还原
    %******************************************************
    %                小尺度度信道生成
    %******************************************************
    % 把瑞利信道单独拿出来生成，因为有时候需要研究小尺度信道的相关性
    %瑞利信道损耗,转置后,矩阵为1*TOTAL_SUB
    lossSmallScare(:,:,iLoop)=XuY_Fun_lossRayleigh(TOTAL_SUB,1)';
    %******************************************************
    %                大小信道合成
    %******************************************************
    % lossLargeScare=zeros(1,TOTAL_USER,LOOP);
    % lossSmallScare=zeros(1,TOTAL_SUB,LOOP);
    % gainChannel = zeros( TOTAL_USER, TOTAL_SUB, LOOP );
    % 前两公式计算出来的是对数，所以这里直接相加
    gainChannel(:,:,iLoop)=(lossLargeScare(:,:,iGain)'*lossSmallScare(:,:,iGain))./(N0*BAND_SUB);
     
end % for iLoop = 1:LOOP

save gainChannel

end % Function

