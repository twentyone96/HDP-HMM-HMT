clc;clear;close all;

global CovA
global re_d
global re_t
global RadarLocX
global RadarLocY
global J   %采样例子数量
global gamma_a
global gamma_b
global alpha_a
global alpha_b
global PolNoMerge
global T
global Pd
global maxMes				%最大观测数
global MaxStateNum  %最大状态数，即截断状态数
global Q_coeff_P %控制U的有效性
global Q_coeff_IMM %控制U的有效性
global Q_coeff_VB %控制U的有效性
global dim
global CartX
global CartY
global PTM                   %
global IMMNUM
global TL

global a_0
global tau_0

%global nStepEnd

CovA =0.1*pi/180; % 1度换算成数字  %default 0.1*pi/180
re_d = 10; %极坐标系下的长度stard 误差  %default 10
re_t = CovA;
RadarLocX = 0; % radar location x %default 0
RadarLocY = 0;%50000; %0;0;50000; % radar location y %default 50000
%tf = 251;%260;
nStepEnd = 251;11;%21%251;
nStep = 250;20;%30%260;
J = 100;
T = 1;
Q = 10 * [1 0 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0 0;
          0 0 0 1 0 0 0 0 0;
          0 0 0 0 1 0 0 0 0;
          0 0 0 0 0 1 0 0 0;
          0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0];   % Q前面的系数可以调整
gamma_a = 1;
gamma_b = 0.1;
alpha_a = 1;
alpha_b = 0.1;
PolNoMerge = 1;
maxMes = 10;
MaxStateNum = 30; %P
TL = 20; %VB  
Q_coeff_IMM = 0.01;%7000;%6000;%5000;%0.01;
Q_coeff_P = 0.01;%7000;%6000;%5000;%0.01;
Q_coeff_VB = 100;%8000;%7000;%6000;%5000;%0.01;
dim = 2;
CartX = 50;  %default 500
CartY = 50;  %default 100\
% 
% noPlot = 0;
% PlotState = 1;
% PlotStep = 1;

% load CurveThree2.mat profile2;
% values of nTarg and T are determined by .mat file
nTarg = 2;
iCycle = 2;
IMMNUM = 4; %四个运动模型，排列：匀速、匀加速、-9、6

T=1;
%Ca CV和CT模型
F_CV=[1,T,0,0,0,0,0,0,0;
   0,1,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0,0;
   0,0,0,1,T,0,0,0,0;
   0,0,0,0,1,0,0,0,0;
   0,0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,1,0,0;
   0,0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0,0];    % 三维笛卡尔坐标系
F_CA=[1,T,T^2/2,0,0,0,0,0,0;
   0,1,T,0,0,0,0,0,0;
   0,0,1,0,0,0,0,0,0;
   0,0,0,1,T,T^2/2,0,0,0;
   0,0,0,0,1,T,0,0,0;
   0,0,0,0,0,1,0,0,0;
   0,0,0,0,0,0,1,0,0;
   0,0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0,0];    % 三维笛卡尔坐标系
H=[1,0,0,1,0,0,1,0,0];
Q_CV =[T^4/4   T^3/2   0    0       0       0       0      0      0;
    T^3/2   T^2     0       0       0       0      0      0      0;
    0       0       0      0        0       0    0     0     0;
    0       0       0     T^4/4   T^3/2     0       0    0     0;
    0       0       0     T^3/2   T^2     0       0      0       0;
    0       0       0      0        0       0    0     0     0;
    0       0       0       0       0       0    0   0     0;
    0       0       0     0      0    0   0   0     0;       
    0       0       0       0       0       0 0  0  0];    % z方向上无速度，所以Q值与x，y方向不同
Q_CA =[T^5/20   T^4/8   T^3/6       0       0       0  0  0  0; 
    T^4/8   T^3/3     T^2/2       0       0       0  0   0   0;
    T^3/6       T^2/2      T      0        0       0    0     0     0;
    0       0       0   T^5/20   T^4/8   T^3/6      0   0   0;
    0       0       0   T^4/8   T^3/3     T^2/2       0    0   0;
    0       0       0     T^3/6       T^2/2      T    0     0     0;
    0       0       0      0        0       0    0     0     0;
    0       0       0       0       0       0    0     0     0;
    0       0       0       0       0       0    0     0     0];    % z方向上无加w2速度，所以Q值与x，y方向不同
% Q=sqrt(0.01)*Q;       % 改变参数
%IMM transition Matrix
%匀速 匀加速 -9  -6
%  101  5   2   2   110
%   5   45  0   0    50
%   2   0   38  0
%   2   0   0   58
PTM = [0.918 0.045 0.0185 0.0185 ;
       0.200 0.800 0.0000 0.0000 ;
       0.050 0.000 0.9500 0.0000 ;
       0.033 0.000 0.0000 0.9670 ];
PTM = [0.927 0.036 0.0185 0.0185 ;
       0.200 0.800 0.0000 0.0000 ;
       0.050 0.000 0.9500 0.0000 ;
       0.033 0.000 0.0000 0.9670 ];

   
X = GenState(nTarg, nStep, T, F_CV, Q_CV, F_CA, Q_CA);

% T = 0.5;    % 开始这里错了，写的0.5，估计出来的速度当然就是两倍了，收敛速度还可以
% assign values of other necessary parameters
densClt = 1e-8;
%densClt = 5e-5;
Pd =0.99;
%Pd = 0.999;          % 改变参数
% Pd = .99;
q = 1;
r = sqrt(1);     % 改变参数，量测噪声方差

M = 4; % number of most probable hypos listed in each scan
N = 3; % height of the hypo tree. It seems scan depth = N-1
% M = 10; % number of most probable hypos listed in each scan
% N = 10; % height of the hypo tree. It seems scan depth = N-1
densNew = 0;
%HDP-HMM MCMC methods
RMS1 = zeros(nStep, nTarg,iCycle);
RMS2 = zeros(nStep, nTarg,iCycle);
RMS1Arv = zeros(nStep,nTarg);
RMS2Arv = zeros(nStep,nTarg);
RMSE1 = zeros(1,nTarg);
RMSE2 = zeros(1,nTarg);

% IMM methods
RMS3 = zeros(nStep, nTarg,iCycle);
RMS4 = zeros(nStep, nTarg,iCycle);
RMS3Arv = zeros(nStep,nTarg);
RMS4Arv = zeros(nStep,nTarg);
RMSE3 = zeros(1,nTarg);
RMSE4 = zeros(1,nTarg);

% HDP-HMM-VB methods
RMS5 = zeros(nStep, nTarg,iCycle);
RMS6 = zeros(nStep, nTarg,iCycle);
RMS5Arv = zeros(nStep,nTarg);
RMS6Arv = zeros(nStep,nTarg);
RMSE5 = zeros(1,nTarg);
RMSE6 = zeros(1,nTarg);

plotState1 = zeros(dim,nStep-2);
plotState2 = zeros(dim,nStep-3);
plotState3 = zeros(dim,nStep-2);
plotState4 = zeros(dim,nStep-3);
plotState5 = zeros(dim,nStep-2);
plotState6 = zeros(dim,nStep-3);
% 调用MHT函数
for i=1:iCycle
[state, meas] = FormatTrans(X, nTarg, nStep, densClt, Pd);
if i == 1
    figure(3);
    legend('真实航迹-','估计航迹o');
    if nTarg == 1
        plot(state{1}(1, 1:nStep-1), state{1}(4, 1:nStep-1), 'k-');
        hold on
        drawnow
    end
    if nTarg == 2
        plot(state{1}(1, 1:nStep), state{1}(4, 1:nStep), 'k-');
        drawnow
        hold on
        plot(state{2}(1, 1:nStep), state{2}(4, 1:nStep), 'k--');
        drawnow
        hold on
        xlabel('x direction');
        ylabel('y direction');
    end
    axis([48000,78000,-10000,20000]);
end
%tic;

[estmMHT,RMS1(:,:,i), RMS2(:,:,i), RMS3(:,:,i), RMS4(:,:,i),RMS5(:,:,i), RMS6(:,:,i),...
    RMSE1,RMSE2,RMSE3,RMSE4,RMSE5,RMSE6,plotState1,plotState2,plotState3,plotState4,plotState5,plotState6,...
    stateEstmIMMRun,stateEstmRun,stateEstmVBRun]...
    = MHT(state, meas, T, nStep, nTarg, densClt, densNew, Pd, M, N,F_CV, Q_CV, F_CA, Q_CA,...
    RMS1(:,:,i),RMS2(:,:,i),RMS3(:,:,i),RMS4(:,:,i),RMS5(:,:,i),RMS6(:,:,i),...
    RMSE1,RMSE2,RMSE3,RMSE4,RMSE5,RMSE6,i,plotState1,plotState2,plotState3,plotState4,plotState5,plotState6);
%toc;
    % 第一次循换，第一次循环的时候把目标估计画出来
    if i == 1
%         plot(plotState1(1,:),plotState1(2,:),'g*',plotState2(1,:),plotState2(2,:),'ro'); %imm
%         plot(plotState3(1,:),plotState3(2,:),'gs',plotState4(1,:),plotState4(2,:),'rp'); %hdp
%         plot(plotState5(1,:),plotState5(2,:),'g+',plotState6(1,:),plotState6(2,:),'rx'); %VB
        hold off
        figure(33);
        plot(state{1}(1, 1:nStep), state{1}(4, 1:nStep), 'k-');
        drawnow
        hold on
        plot(state{2}(1, 1:nStep), state{2}(4, 1:nStep), 'k--');
        drawnow
        hold on
        plot(stateEstmIMMRun(:,2,1),stateEstmIMMRun(:,2,2),'bo');
        plot(stateEstmRun(:,1,1),stateEstmRun(:,1,2),'bp');
        plot(stateEstmVBRun(:,1,1),stateEstmVBRun(:,1,2),'b>');
        drawnow
        hold on
        plot(stateEstmIMMRun(:,1,1),stateEstmIMMRun(:,1,2),'ro');
        plot(stateEstmRun(:,2,1),stateEstmRun(:,2,2),'rp');
        plot(stateEstmVBRun(:,2,1),stateEstmVBRun(:,2,2),'r>');
        drawnow
        hold on
        xlabel('x direction');
        ylabel('y direction');
        axis([48000,78000,-10000,20000]);
        legend('True trajectory, target 1','True trajectory, target 2',...
            'IMM-MHT target 1','HDP-HMM-MHT target 1','HDP-HMM-PMHT target 1','IMM-MHT target 2',...
            'HDP-HMM-MHT target 2','HDP-HMM-PMHT target 2');
        drawnow
        hold on
        hold off
    end
    % 这是最后一次循环，计算多轮循环的均方值
    if i == iCycle
        for j=1:nStep
            for k=1:iCycle
                RMS1Arv(j,1) = RMS1(j,1,k)/iCycle +  RMS1Arv(j,1);
                RMS1Arv(j,2) = RMS1(j,2,k)/iCycle +  RMS1Arv(j,2);
                RMS2Arv(j,1) = RMS2(j,1,k)/iCycle +  RMS2Arv(j,1);
                RMS2Arv(j,2) = RMS2(j,2,k)/iCycle +  RMS2Arv(j,2);

                RMS3Arv(j,1) = RMS3(j,1,k)/iCycle +  RMS3Arv(j,1);
                RMS3Arv(j,2) = RMS3(j,2,k)/iCycle +  RMS3Arv(j,2);
                RMS4Arv(j,1) = RMS4(j,1,k)/iCycle +  RMS4Arv(j,1);
                RMS4Arv(j,2) = RMS4(j,2,k)/iCycle +  RMS4Arv(j,2);

                RMS5Arv(j,1) = RMS5(j,1,k)/iCycle +  RMS5Arv(j,1);
                RMS5Arv(j,2) = RMS5(j,2,k)/iCycle +  RMS5Arv(j,2);
                RMS6Arv(j,1) = RMS6(j,1,k)/iCycle +  RMS6Arv(j,1);
                RMS6Arv(j,2) = RMS6(j,2,k)/iCycle +  RMS6Arv(j,2);

            end
            RMS1Arv(j,:) = sqrt(RMS1Arv(j,:));
            RMS2Arv(j,:) = sqrt(RMS2Arv(j,:));
            RMS3Arv(j,:) = sqrt(RMS3Arv(j,:));
            RMS4Arv(j,:) = sqrt(RMS4Arv(j,:));
            RMS5Arv(j,:) = sqrt(RMS5Arv(j,:));
            RMS6Arv(j,:) = sqrt(RMS6Arv(j,:));
        end
        %% iCycle==1时有效，好像没有用到
        RMSE1 = sqrt(RMSE1/((nStep-3)*iCycle))
        RMSE2 = sqrt(RMSE2/((nStep-3)*iCycle))
        RMSE3 = sqrt(RMSE3/((nStep-3)*iCycle))
        RMSE4 = sqrt(RMSE4/((nStep-3)*iCycle))
        RMSE5 = sqrt(RMSE5/((nStep-3)*iCycle))
        RMSE6 = sqrt(RMSE6/((nStep-3)*iCycle))
        
        figure(5);
        hold on
%        plot(4:nStep,RMS1Arv(4:nStep,1),'g--',4:nStep,RMS1Arv(4:nStep,2),'r-');
        plot(4:nStep,RMS1Arv(4:nStep,1),'g--',4:nStep,RMS1Arv(4:nStep,2),'r-');
        plot(4:nStep,RMS4Arv(4:nStep,1),'k--',4:nStep,RMS4Arv(4:nStep,2),'b-');
        plot(4:nStep,RMS6Arv(4:nStep,1),'m--',4:nStep,RMS6Arv(4:nStep,2),'y-');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,0,2000]),title('RMS error of position (m)')
        xlabel('time');
        ylabel('Position RMSE');
        hold off
 
        figure(6);
        hold on
%        plot(4:nStep,RMS2Arv(4:nStep,1),'g--',4:nStep,RMS2Arv(4:nStep,2),'r-');
        plot(4:nStep,RMS2Arv(4:nStep,1),'g--',4:nStep,RMS2Arv(4:nStep,2),'r-');
        plot(4:nStep,RMS3Arv(4:nStep,1),'k--',4:nStep,RMS3Arv(4:nStep,2),'b-');
        plot(4:nStep,RMS5Arv(4:nStep,1),'m--',4:nStep,RMS5Arv(4:nStep,2),'y-');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,0,2000]),title('RMS error of position (m)')
        xlabel('time');
        ylabel('Position RMSE');
        hold off
    end
end


%% 这里设计MHT_HDP的分析程序
%for i = 1 : length(estmMHT)
%    estmMHT{i} = estmMHT{i}{3}; % extract the state
%end
% calculate mean error
%first = 2; % estmMHT includes accurate initial state
%last = (nStep+1) - (N-1); 
%[errRMS, lose] = Analyse(first, last, estmMHT, state, nTarg);
%rstMHT = [errRMS lose]; % MHT tracking results
%disp('    err_x     err_vx     err_y     err_vy     lose');
%disp(rstMHT);
fprintf('MHT DONE\n\n');

%参考
