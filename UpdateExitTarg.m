function aTarg = UpdateExitTarg(Targ, Z, t)
global J
global T
global MaxStateNum
global dim
global Q_coeff_P

%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u的维度，这里选2位） nu(16) Delta(17)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) 每个状态保持的数量
%% cellCtrAccum(19)    控制量的累加\sum_1^N U*U'

		kappaGlobal = 0.001;%0.001;
		varthetaGlobal = zeros(dim,1); %[141,141]';%
		nuGlobal = 50; 
		DeltaGlobal = eye(dim);

		y = Z(1:2);
		FU = GenFU(T);
		A = FU(1:6,1:6);
		BU = GenBU(T);
		B = BU(1:6,1:2);
		C = [1 0 0 0 0 0; 0 0 0 1 0 0];
		y = Z(1:2,1);
		RCov = GenCov(Z);
		R = RCov(1:2,1:2);
		Q = Q_coeff_P*eye(6,6);

		alpha_a = 1;
		alpha_b = 0.1;
		gamma_a = 1;
		gamma_b = 0.1;
		sampleIdx = 1;
		[sampleWeight, sampleZWeight,tTarg] = CalcMargCost(Targ{3},Z);
        Targ{3} = tTarg;
		tempTarg = cell(1,J);
%         sampleSum = zeros(MaxStateNum,1);
%         sampleSum = sum(sampleZWeight,2);
%         for jjjj = 1 : MaxStateNum
%             sampleZWeight(jjjj,:) = sampleZWeight(jjjj,:)/sampleSum(jjjj);
%         end
		for i = 1 : J
% 		    for jj = 1 : sampleWeight(i)
%            thisTarg = Targ{3}{i};
            [thisTarg, sampleZW] = removeRedundant(Targ{3}{i},sampleZWeight(:,i));
            sampleZWeight(:,i) = sampleZW;
            totalL = thisTarg{4};
            X = thisTarg{1};
            P = thisTarg{2};            
            p = zeros(totalL+1, 1);
            currentS = thisTarg{3};

            %% 计算比例
            p = cumsum(sampleZWeight(1:totalL+1,i));
            u = rand()*p(totalL+1,1);
            %% j就是采样出的索引
            for jjj = 1 : totalL + 1
                if p(jjj,1) > u
                    break;
                end
            end
            L = jjj;
            thisTarg{3} = L;
            cellTrans = thisTarg{6};
            cellCtr = thisTarg{5};
            cellStates = thisTarg{18};
            cellCtrAccum = thisTarg{19};
            cellCtrLast = thisTarg{20};
            cellCtrLastS = thisTarg{21};
            %% 计算值
            gammaV = thisTarg{9};
            Beta = thisTarg{8};
            Alpha = thisTarg{7};

            K = inv(Q) + C'*inv(R)*C;
            S = inv(Q) - inv(Q)*inv(K)*inv(Q);
            LambdaTmp11 = B'*S*B;
            LambdaTmp12 = B'*S*A;
            KK = inv(Q)*A*X + C'*inv(R)*y;

            thetaTmp1 = B'*inv(Q)*inv(K)*C'*inv(R)*y;

            Lambda = LambdaTmp11;
            theta = thetaTmp1 - LambdaTmp12*X;

            Lambda = B'*S*B-B'*S*A*inv(A'*S*A+inv(P))*A'*S*B;
            theta = B'*inv(Q)*inv(K)*C'*inv(R)*y...
                -B'*S*A*inv(A'*S*A+inv(P))...
                *(inv(P)*X+A'*inv(Q)*inv(K)*C'*inv(R)*y);
            if jjj == totalL + 1 
                totalL = totalL + 1;
                thisTarg{4} = totalL;
                kappa = thisTarg{14};
                vartheta = thisTarg{15};
                nu = thisTarg{16};
                Delta = thisTarg{17};
                mu = vartheta;
                sigma = ((kappa + 1)*nu/(kappa*(nu - dim - 1)))*Delta;
                
                
                mu = cellCtrLast;  %% 影响也不是太大
                sigma = cellCtrLastS; %% 影响也不是太大

%                  mu = inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta);
%                  sigma = inv(inv(sigma)+Lambda);
                muSample = mvnrnd(mu', sigma, 1)';
%                muSample = mu; %% 影响也不是太大

                
                cellCtr(:,jjj) = cellCtr(:,jjj) + muSample;
                v0 = random('beta', gammaV, 1);
                BetaTmp = Beta(jjj);
                Beta(jjj:jjj+1) = [BetaTmp*(1-v0), BetaTmp*v0];
                cellStates(jjj) = 1;
                cellCtrAccumTmp = muSample*muSample';
                cellCtrAccum(jjj,:) =...
                    [cellCtrAccumTmp(1,:),cellCtrAccumTmp(2,:)];

            else
                kappa = thisTarg{14} + thisTarg{18}(jjj);
                vartheta = (thisTarg{14}*thisTarg{15} + thisTarg{5}(:,jjj))...
                    /kappa;
                nu = thisTarg{16} + thisTarg{18}(jjj);
                Delta = (thisTarg{16}*thisTarg{17}...
                    + [thisTarg{19}(jjj,1:2);thisTarg{19}(jjj,3:4)]...
                    + thisTarg{14}*(thisTarg{15}*thisTarg{15}')...
                    -kappa*(vartheta*vartheta'))/nu;
                mu = vartheta;
                sigma = ((kappa + 1)*nu/((kappa*(nu - dim - 1))))*Delta;

                mu = inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta);
                sigma = inv(inv(sigma)+Lambda);                
                
                muSample = mvnrnd(mu', sigma, 1)';
%                muSample = mu;
                cellCtr(:,jjj) = cellCtr(:,jjj) + muSample;
                cellStates(jjj) = cellStates(jjj) + 1;
                cellCtrAccumTmp = muSample*muSample';
                cellCtrAccum(jjj,:) = cellCtrAccum(jjj,:)...
                    + [cellCtrAccumTmp(1,:),cellCtrAccumTmp(2,:)];
            end
                        
            thisTarg{5} = cellCtr;
            cellTrans(currentS,jjj) = cellTrans(currentS,jjj) + 1;
            thisTarg{6} = cellTrans;
            thisTarg{18} = cellStates;
            thisTarg{19} = cellCtrAccum;
            
            X = A*X+B*muSample;
            P = A*P*A' + B * sigma * B' + Q;

%simple
%            X = A*X+B*mu;
%            P = A*P*A' + Q;
            innov = y - C*X;
            SPredict = C*P*C' + R;
            G = P*C'*inv(SPredict);
            X = X + G*innov;
            P = (eye(6,6) - G*C)*P;
            thisTarg{1} = X;
            thisTarg{2} = P;
            
            M = zeros(totalL);
            for j=1:totalL
                for k=1:totalL
                    if cellTrans(j,k) == 0
                        M(j,k) = 0;
                    else
                        for l=1:cellTrans(j,k)
                            M(j,k) = M(j,k)...
                                + (rand() < (Alpha * Beta(k))...
                                / (Alpha * Beta(k) + l - 1));
                        end
                    end
                end
            end
            thisTarg{10}(1:totalL,1:totalL) = M;
            
           phi = random('beta',gammaV + 1, sum(sum(M,1),2));
           thisTarg{11} = phi;
           auxVarepsilon = (gamma_a + totalL - 1)/...
               (sum(sum(M,1),2)*(gamma_b - log(phi)) + gamma_a + totalL - 1);
            gammaV = auxVarepsilon*random('gam',alpha_a + totalL,...
               alpha_b - log(phi)) + ...
               (1 - auxVarepsilon)*random('gam',alpha_a + totalL - 1,...
               alpha_b - log(phi));
           thisTarg{9} = gammaV;
           auxTau = random('beta', Alpha + 1, sum(cellTrans(1:totalL,:), 2));
           auxZeta_p = sum(cellTrans(1:totalL,:), 2)./(Alpha + sum(cellTrans(1:totalL,:), 2));
           auxZeta = rand(totalL,1);
           auxZeta = auxZeta < auxZeta_p;
           Alpha = random('gam',alpha_a + sum(sum(M,1),2) - sum(auxZeta,1),...
               alpha_b - sum(log(auxTau),1));
           thisTarg{7} = Alpha;
           %M 的列相加
           Beta = randdirichlet([sum(M(1:totalL,1:totalL), 1),gammaV]')';
           thisTarg{8} = [Beta, zeros(1,1+MaxStateNum-size(Beta,2))];
           thisTarg{14} = kappaGlobal;
           thisTarg{15} = varthetaGlobal;
           thisTarg{16} = nuGlobal;
           thisTarg{17} = DeltaGlobal;
            
           tempTarg{sampleIdx} = thisTarg;
           sampleIdx = sampleIdx + 1;
%		   end
		end
		Targ{3} = tempTarg;
		aTarg = Targ;
end
% 建议一次去掉两个状态，每次进行状态缩减都比较耗费时间
function [aTarg, sampleZW] = removeRedundant(thisTarg, sampleZW)
global MaxStateNum
    aTarg = thisTarg;
    % 最大状态等于MaxStateNum
    if aTarg{4} == MaxStateNum
        iterNum = ceil(MaxStateNum/10);
        minIndex = zeros(iterNum,1);
        tempIndex = aTarg{18};
        temp3 = aTarg{3};
        j = 1;
        [value, minIndex(j)] = min(tempIndex);
        % 如果最小的状态等于当前的状态，去掉下一个最小状态
        while  j <= iterNum
            tempIndex(minIndex(j),1) = inf;            
            if minIndex(j) ~= temp3 && j == iterNum
                break;
            end
            j = j + 1;
            [value, minIndex(j)] = min(tempIndex);
        end
        for i = 1 : iterNum
            sampleZW(minIndex(i),:) = [];
            sampleZW = [sampleZW; 0];
            aTarg{18}(minIndex(i),:) = []; 
            aTarg{18} = [aTarg{18};0]; 
            if aTarg{3} >= minIndex(i) 
                aTarg{3} = aTarg{3} - 1; %上一个状态的索引减去一
            end
            aTarg{4} = aTarg{4} - 1; %状态数减去一
            aTarg{5}(:,minIndex(i)) = []; %去掉保存的状态值 cellCtr(:,jjj)
            aTarg{5} = [aTarg{5} zeros(2,1)];
            aTarg{6}(:,minIndex(i)) = []; %去掉转移的计数值 cellTrans(currentS,jjj)，矩阵的横宽维数相等，相当于去掉一根十字线
            aTarg{6}(minIndex(i),:) = []; %去掉转移的计数值 cellTrans(currentS,jjj)，矩阵的横宽维数相等，相当于去掉一根十字线
            oneLine = zeros(MaxStateNum-1,1);
            aTarg{6} = [aTarg{6} oneLine];
            oneLine = zeros(1,MaxStateNum);
            aTarg{6} = [aTarg{6}; oneLine];
    %        aTarg{7}   %Alpha是一个变量，应该不管了
            tempBeta = aTarg{8}(minIndex(i)); %Beta是一个向量，属于dirichlet分布，故向量和为1
            aTarg{8}(:,minIndex(i)) = [];
            aTarg{8}(end) = aTarg{8}(end) + tempBeta;
            aTarg{8} = [aTarg{8} 0];
    %        aTarg{9}   %该变量没办法重新算，应该不管了
            aTarg{10}(:,minIndex(i)) = [];
            aTarg{10}(minIndex(i),:) = [];
            oneLine = zeros(MaxStateNum-1,1);
            aTarg{10} = [aTarg{10} oneLine];
            oneLine = zeros(1,MaxStateNum);
            aTarg{10} = [aTarg{10}; oneLine];
            aTarg{19}(minIndex(i),:) = [];
            aTarg{19} = [aTarg{19};zeros(1,4)];
        end
    end
end