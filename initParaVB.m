function [newTarg, newState, ZetaCD, VarepsilonAB, Xi, Varphi] = initParaVB(state)
%INITPARA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
global dim
global TL %truncation level
global TargetNum %Ŀ�����
kappaGlobal = 0.001;%0.001;
varthetaGlobal = zeros(dim,1);%[141,141]';%
nuGlobal = 50; 
DeltaGlobal = eye(dim);
%��һ��G�Ĳ���a_0
a_0 = 1;
%�ڶ���G�Ĳ���tau_0
tau_0 = 1;
%��������һ��G_0�͵ڶ���G_j�Ľضϸ�����Ϊ10
% ��ʼ��Ŀ����˶����ƺͷ���
%% X(1) P(2) currentState(3) cellCtr(4) ��mu...
%% cellTrans(5) ����ת�Ƹ��ʾ�����Sigma
%% kappa(6) vartheta(7)(dim=u��ά�ȣ�����ѡ2λ�� nu(8) Delta(9)��2*2����
%% --��һ�ж���NIW�Ĳ���
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(10) ÿ��״̬���ֵ�����
%% cellCtrAccum(11) ÿ��״̬U*U'
%% targPi(12) ����Ŀ��Ĺ۲�������,��Ӧ���ĵ�Pi_tm
%% cellCtr(13)��ԭ����5
%% lastState(14) ��һ������������״̬
%%%%%%%%%%%%%%%%%%%%%%%%%ԭʼ�Ķ���
%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u��ά�ȣ�����ѡ2λ�� nu(16) Delta(17)��2*2����
%% --��һ�ж���NIW�Ĳ���
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) ÿ��״̬���ֵ�����
%% cellCtrAccum(19) ÿ��״̬U*U'
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
    
%     newTarg{3,i}.stateP = zeros(1,1,'int8'); %��ǰ������״̬
%     newTarg{4,i}.mu = zeros(2,TL);    %����
%     newTarg{5,i}.sigma = zeros(4,TL); %����
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

    aTarg{3} = zeros(1,1,'int8'); %��ǰ������״̬
    aTarg{4} = zeros(2,TL);    %����
    aTarg{5} = zeros(4,TL); %����
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
% ��ʼ����ǰ����һ��ʱ������״̬�ĸ���
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

% ��ʼ������Ŀ��Ĺ۲�������,��Ӧ���ĵ�Pi_tm
%targPi = zeros(1,TargetNum);
for i = 1 : TargetNum
%    newTarg{12,i}.targPi = 0.5;
%     newTarg{12,i} = 0.5;
    newTarg{i}{12} = 0.5;
end

% ��ʼ������Ŀ��������״̬stateP��״̬ת�Ƹ��ʾ���stateTranP,��Ӧ���ĵ�Z_tm
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
%��ʼ��\Xi,\Varphi
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
%��ʼ��\Zeta, \Zeta_c, \Zeta_d  ��Ӧ���������Beta
ZetaCD = zeros(TL,2);
for i = 1 : TL
    zetaVarphi = getVarphi(Varphi,i);
    ZetaCD(i,1) = 1 + zetaVarphi(1);
    ZetaCD(i,2) = a_0 + zetaVarphi(2);
end
%��ʼ��\Varepsilon, \Varepsilon_a, \Varepsilon_b ��Ӧ���������Pi
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
% ��ʼ������Ŀ����˶�������������
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

