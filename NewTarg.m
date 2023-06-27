function newTarg = NewTarg(Z)
global J
global MaxStateNum
global dim
%                 aTarg{3} = [aMeas(1), 0, aMeas(2), 0]'; % X
%                 aTarg{4} = diag([1 1 1 1]); % P
%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u的维度，这里选2位） nu(16) Delta(17)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) 每个状态保持的数量

	kappaGlobal = 0.001;%0.001;
	varthetaGlobal = zeros(dim,1);%[141,141]';%
	nuGlobal = 50; 
	DeltaGlobal = eye(dim);
	alpha_a = 1;
	alpha_b = 0.1;
	gamma_a = 1;
	gamma_b = 0.1;
    RCov = GenCov(Z);
    R = RCov(1:2,1:2);

	y = Z(1:2);
	newTarg = cell(1, J);
	%% 初始化每一个粒子
	for i = 1 : J
			aSample = cell(1,21);
			aSample{1} = zeros(1,6);
			aSample{2} = eye(6,6);
			aSample{3} = zeros(1,1,'int8');
			aSample{4} = zeros(1,1,'int8');
			aSample{5} = zeros(2,MaxStateNum);
			aSample{6} = zeros(MaxStateNum,MaxStateNum);
			aSample{7} = zeros(1,1);
			aSample{8} = zeros(1,1+MaxStateNum);
			aSample{9} = zeros(1,1);
			aSample{10} = zeros(MaxStateNum,MaxStateNum); 
			aSample{11} = zeros(1,1);
			aSample{12} = zeros(MaxStateNum,1);
			aSample{13} = zeros(MaxStateNum,1);
			aSample{14} = zeros(1,1);
			aSample{15} = zeros(dim,1);
			aSample{16} = zeros(1,1);
			aSample{17} = zeros(dim,dim);
			aSample{18} = zeros(MaxStateNum,1);
			aSample{19} = zeros(MaxStateNum,4);
			aSample{20} = zeros(2,1);
			aSample{21} = zeros(2,2);
	    X = [Z(1), 0, 0, Z(2), 0, 0]';
	    P = eye(6,6);
        P(1,1) = R(1,1);
        P(1,4) = R(1,2);
        P(4,1) = R(2,1);
        P(4,4) = R(2,2);
             
	    aSample{1} = X;
	    aSample{2} = P;
	    currentState = 1;
	    L = 1;
	    aSample{3} = currentState;
	    aSample{4} = L;
	    
	    mu = varthetaGlobal;
	    sigma = ((kappaGlobal + 1)*nuGlobal/(kappaGlobal*(nuGlobal - dim - 1)))*DeltaGlobal;
	    muSample = mvnrnd(mu', sigma,1)';
	    aSample{5}(:,1) = muSample;
	    
	    aSample{6}(1,1) = 1;
	   
	    gamma = random('gamma', alpha_a, alpha_b);
	    aSample{9} = gamma;
	    alpha = random('gamma', gamma_a, gamma_b);
	    aSample{7} = alpha;
        v0 = 0;
        while v0 == 0
	    v0 = random('beta', gamma, 1);
        end
	    Beta = [1-v0,v0];
	    aSample{8}(1:2) = Beta;
	    aSample{10}(1,1) = 1;
	    
	    aSample{11} = random('beta',gamma + 1, 1);
	    
	    aSample{14} = kappaGlobal;
	    aSample{15} = varthetaGlobal;
	    aSample{16} = nuGlobal;
	    aSample{17} = DeltaGlobal;
	    aSample{18}(1,1) = 1;
	    cellCtrAccumTmp = muSample*muSample';
	    cellCtrAccum = [cellCtrAccumTmp(1,:),cellCtrAccumTmp(2,:)];
	    aSample{19}(1,:) = cellCtrAccum;
	    newTarg{i} = aSample;
	end
end
