function [newTarg, newState, ZetaCD, VarepsilonAB, Xi, Varphi] = initParaVB(state)
%INITPARA 此处显示有关此函数的摘要
%   此处显示详细说明
global dim
global TL %truncation level
global TargetNum %目标个数
kappaGlobal = 0.001;%0.001;
varthetaGlobal = zeros(dim,1);%[141,141]';%
nuGlobal = 50; 
DeltaGlobal = eye(dim);
%第一级G的参数a_0
a_0 = 1;
%第二级G的参数tau_0
tau_0 = 1;
%这里假设第一级G_0和第二级G_j的截断个数都为10
% 初始化目标的运动估计和方差
%% X(1) P(2) currentState(3) cellCtr(4) 是mu...
%% cellTrans(5) 不是转移概率矩阵，是Sigma
%% kappa(6) vartheta(7)(dim=u的维度，这里选2位） nu(8) Delta(9)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(10) 每个状态保持的数量
%% cellCtrAccum(11) 每个状态U*U'
%% targPi(12) 各个目标的观测分配概率,对应论文的Pi_tm
%% cellCtr(13)是原来的5
%% lastState(14) 上一个周期所处的状态
%%%%%%%%%%%%%%%%%%%%%%%%%原始的定义
%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u的维度，这里选2位） nu(16) Delta(17)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) 每个状态保持的数量
%% cellCtrAccum(19) 每个状态U*U'
%newTarg = cell(10,TargetNum);
newTarg = cell(1,TargetNum);
for i = 1 : TargetNum
    aTarg = cell(1,17);
    X = zeros(dim*3,1);
    P = zeros(dim*3,dim*3);
%     newTarg{1,i}.X = state{i}(1:6,1);   %X;
%     newTarg{2,i}.P = [10,0,0,0,0,0;
%                       0,10,0,0,0,0;
%                       0,0,0,0,0,0;
%                       0,0,0,10,0,0;
%                       0,0,0,0,10,0;
%                       0,0,0,0,0,0];    %P;
    aTarg{1} = X;
    aTarg{1}(1) = state{i}(1);   %X;
    aTarg{1}(4) = state{i}(4);   %X;
    aTarg{2} = [1,0,0,0,0,0;
                      0,100,0,0,0,0;
                      0,0,100,0,0,0;
                      0,0,0,1,0,0;
                      0,0,0,0,100,0;
                      0,0,0,0,0,100];    %P;
    
%     newTarg{3,i}.stateP = zeros(1,1,'int8'); %当前所处的状态
%     newTarg{4,i}.mu = zeros(2,TL);    %不用
%     newTarg{5,i}.sigma = zeros(4,TL); %不用
%     newTarg{6,i}.ikappa = zeros(1,1); 
%     newTarg{7,i}.ivartheta = zeros(dim,1);
%     newTarg{8,i}.inu = zeros(1,1);
%     newTarg{9,i}.idelta = zeros(dim,dim);
%     newTarg{10,i}.cellStates = zeros(TL,1);
%     newTarg{11,i}.cellCtrAccum = zeros(TL,4);
%     newTarg{12,i}.targPi = 0;
%     newTarg{13,i}.cellCtr = zeros(2,TL);
%     newTarg{14,i}.currentState = zeros(10,1);
%     newTarg{15,i}.lastState = zeros(10,1);
%     newTarg{16,i}.probZ = zeros(TL,1);
%     newTarg{17,i}.lastProbZ = zeros(TL,1);

    aTarg{3} = zeros(1,1,'int8'); %当前所处的状态
    aTarg{4} = zeros(2,TL);    %不用
    aTarg{5} = zeros(4,TL); %不用
    aTarg{6} = zeros(1,1); 
    aTarg{7} = zeros(dim,1);
    aTarg{8} = zeros(1,1);
    aTarg{9} = zeros(dim,dim);
    aTarg{10} = zeros(TL,1);
    aTarg{11} = zeros(TL,4);
    aTarg{12} = 0;
    aTarg{13} = zeros(2,TL);
    aTarg{14} = zeros(10,1);
    aTarg{15} = zeros(10,1);
    aTarg{16} = zeros(TL,1);
    aTarg{17} = zeros(TL,1);
    newTarg{i} = aTarg;
end
for j = 1 : TargetNum
    for i = 1 : TL
%         newTarg{16,j}.probZ(i) = 1/TL;
%         newTarg{17,j}.lastProbZ(i) = 1/TL;

        newTarg{j}{16}(i) = 1/TL;
        newTarg{j}{17}(i) = 1/TL;
    end
end
% 初始化当前和上一个时刻所处状态的概率
currentState = zeros(TL,1);
lastState = zeros(TL,1);
for j = 1 : TargetNum
    for i = 1 : TL
        currentState(i) = 1/TL;
        lastState(i) = 1/TL;
    end
%     newTarg{14,j}.currentState = currentState;
%     newTarg{15,j}.lastState = lastState;

    newTarg{j}{14} = currentState;
    newTarg{j}{15} = lastState;
end

% 初始化各个目标的观测分配概率,对应论文的Pi_tm
%targPi = zeros(1,TargetNum);
for i = 1 : TargetNum
%    newTarg{12,i}.targPi = 0.5;
%     newTarg{12,i} = 0.5;
    newTarg{i}{12} = 0.5;
end

% 初始化各个目标所处的状态stateP和状态转移概率矩阵stateTranP,对应论文的Z_tm
%stateP = zeros(TL, TargetNum);

%% (1) mu ; (2) Sigma ; (3) 
% newState = cell(10,TL);
newState = cell(1,TL);
stateTranP = zeros(TL, TL);
for i = 1 : TL
    for j = 1 : TL
       stateTranP(i,j) = 1/TL;
    end
end
%初始化\Xi,\Varphi
Xi = zeros(TL,TL,TargetNum);
Varphi = zeros(TL,TL,TL);
for j = 1 : TargetNum
    for i = 1 : TL %,stick t
        for ii = 1 : TL %state_t-1
            Xi(ii,i,j) = 1/TL;
        end
    end
end
for i = 1 : TL %state_t
    for ii = 1 : TL %stick_t
        for iii = 1 : TL %state_t-1
            Varphi(iii,ii,i) = 1/TL;
        end
    end
end
%初始化\Zeta, \Zeta_c, \Zeta_d  对应文章里面的Beta
ZetaCD = zeros(TL,2);
for i = 1 : TL
    zetaVarphi = getVarphi(Varphi,i);
    ZetaCD(i,1) = 1 + zetaVarphi(1);
    ZetaCD(i,2) = a_0 + zetaVarphi(2);
end
%初始化\Varepsilon, \Varepsilon_a, \Varepsilon_b 对应文章里面的Pi
VarepsilonAB = zeros(TL,TL,2);
for i = 1 : TL
    for ii = 1 : TL
        epsilonVarphi = getXi(Xi,i,ii);
        VarepsilonAB(i,ii,1) = 1 + epsilonVarphi(1);
        VarepsilonAB(i,ii,2) = tau_0 + epsilonVarphi(2);
    end
end
for i = 1 : TL
    aState = cell(1,7);
%     newState{1,i}.mu = zeros(2,TL);    
%     newState{2,i}.sigma = zeros(2,2,TL); 
%     newState{3,i}.stateTranP = stateTranP;
    aState{1} = zeros(2,1);    
    aState{2} = zeros(2,2); 
    aState{3} = stateTranP;
    aState{4} = 0.001;%zeros(1,1); 
    aState{5} = zeros(dim,1);
    aState{6} = 10;%zeros(1,1);
    aState{7} = eye(dim,dim);
    aState{8} = zeros(2,1);
    aState{9} = zeros(1,4);
    newState{i} = aState;
end
%%%%%%%%%
% kappaGlobal = 0.001;%0.001;
% varthetaGlobal = zeros(dim,1);%[141,141]';%
% nuGlobal = 10;%50; 
% DeltaGlobal = eye(dim);

%     newState{4,i}.Xi = Xi;
%     newState{5,i}.Varphi = Varphi;
%     newState{6,i}.Zeta = Zeta{i,:};
%     newState{7,i}.Varepsilon = Varepsilon{i,:,:};
% end
% 初始化各个目标的运动控制量及方差
%Ucontrol = cell(TL, TargetNum);
for j = 1 : TL
%    mu = varthetaGlobal - 40;
    mu = varthetaGlobal - 180;
%    mu = varthetaGlobal;
    sigma = ((kappaGlobal + 1)*nuGlobal/(kappaGlobal*(nuGlobal - dim - 1)))*DeltaGlobal;
    muSample = mvnrnd(mu, sigma,1)';
%     newState{j}{1} = muSample;%mu;%muSample;
%     newState{j}{2} = sigma;
    %_1
    newState{j}{1} = mu + (j-1)*20;    
    newState{j}{2} = sigma;
end

end

