%% 这个版本在预测步和更新步都用了复合估计
%% 参考 多机动目标跟踪算法与仿真平台开发 P64 (3) (4)
%% (3)里面的权重，用状态z的先验估计 probIdx
%% (4)为每个状态z计算预测状态的X和P

%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u的维度，这里选2位） nu(16) Delta(17)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) 每个状态保持的数量
function [sampleWeight, sampleZWeight, newTarg] = CalcMargCost(thisTarg, Z)

global T
global J
global MaxStateNum
global Q_coeff_P
global dim

    FU = GenFU(T);
    A = FU(1:6,1:6);
    BU = GenBU(T);
    B = BU(1:6,1:2);
    C = [1 0 0 0 0 0; 0 0 0 1 0 0];
    y = Z(1:2,1);
    RCov = GenCov(Z);
    R = RCov(1:2,1:2);
    Q = Q_coeff_P*eye(6,6);

    sampleWeight = zeros(J,1);
    sampleZWeight = zeros(MaxStateNum + 1, J);
    %% 先计算粒子的权重，重采样粒子
    probIdx = zeros(MaxStateNum + 1,J);
    probC = zeros(J,1);
    pMatJ1 = zeros(MaxStateNum + 1,J);
    pMax = zeros(J,1);
    testObsCtrl = zeros(dim, J);
    testMixCtrl = zeros(dim, J);
    testCtrlCtrl = zeros(dim,J);
    testX1Ctrl = zeros(dim,J);
    for i = 1 : J
        totalL = thisTarg{i}{4};    
        currentS = thisTarg{i}{3};
        p = zeros(MaxStateNum+1, 1);
        pSum = 0;
        for j = 1 : totalL + 1
            if j == totalL +1
                p(j,1) = thisTarg{i}{7} * thisTarg{i}{8}(j);
            else
                p(j,1) = thisTarg{i}{7} * thisTarg{i}{8}(j)...
                    + thisTarg{i}{6}(currentS,j);
            end
            p(j,1) = round(p(j,1),5);
            pSum = p(j,1) + pSum;
        end
        for j = 1 : totalL + 1
            probIdx(j,i) = p(j,1)/pSum;
        end
    end
    for i = 1 : J
        totalL = thisTarg{i}{4};
        X = thisTarg{i}{1};
        P = thisTarg{i}{2};

        
        cellCtrLast = thisTarg{i}{20};
        for j = 1 : totalL + 1
            if j == totalL + 1
                kappa = thisTarg{i}{14};
                vartheta = thisTarg{i}{15};
                nu = thisTarg{i}{16};
                Delta = thisTarg{i}{17};
                mu = vartheta;
                sigma = ((kappa + 1)*nu/(kappa*(nu - dim - 1)))*Delta;
            else
                kappa = thisTarg{i}{14} + thisTarg{i}{18}(j);
                vartheta = (thisTarg{i}{14}*thisTarg{i}{15} + thisTarg{i}{5}(:,j))...
                    /kappa;
                nu = thisTarg{i}{16} + thisTarg{i}{18}(j);
                Delta = (thisTarg{i}{16}*thisTarg{i}{17}...
                    + [thisTarg{i}{19}(j,1:2);thisTarg{i}{19}(j,3:4)]...
                    + thisTarg{i}{14}*(thisTarg{i}{15}*thisTarg{i}{15}')...
                    -kappa*(vartheta*vartheta'))/nu;
                mu = vartheta;
                sigma = ((kappa + 1)*nu/(kappa*(nu - dim - 1)))*Delta;
            end
            if j==1
            testCtrlCtrl(:,i) = mu;
            end

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

            mu1 = inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta);
            sigma1 = inv(inv(sigma)+Lambda);
            if j == totalL + 1
                cellCtrLast = mu1;
                cellCtrLastS = sigma1;
                thisTarg{i}{20} = cellCtrLast;
                thisTarg{i}{21} = cellCtrLastS;
            end
            mu = mu1;
            sigma = sigma1;
            
    %          x1 = (sqrt(det(inv(inv(sigma)+Lambda)))/((...
    %                 sqrt(det(sigma)))))...
    %                 * exp(-0.5*(mu'*inv(sigma)*mu + theta'*inv(Lambda)*theta...
    %                 - (inv(sigma)*mu+theta)'*inv(inv(sigma))*(inv(sigma)*mu+theta)));%...
    %          x1 = (sqrt(det(inv(inv(sigma)+Lambda)))/(((2*pi)*...
    %                 sqrt(det(sigma)))*sqrt(det((Lambda)))))...
    %                 * exp(-0.5*(mu'*inv(sigma)*mu + theta'*inv(Lambda)*theta...
    %                 - (inv(sigma)*mu+theta)'*inv(inv(sigma))*(inv(sigma)*mu+theta)));%...

    %                *(inv(sigma)*mu+theta))); %+y'*inv(R)*y); 
    %                *sqrt(inv(Q))*sqrt(inv(R))))...
    %        else
        %% 取负log
    %         x2 = 0.05*(dim*log(2*pi)+log(det(sigma))+log(det(Lambda))...
    %             -log(det(inv(inv(sigma)+Lambda)))...
    %             +dim*log(2*pi)+dim*log(2*pi)+log(det(inv(Q)))++log(det(inv(R)))...
    %             + mu'*inv(sigma)*mu + theta'*inv(Lambda)*theta...
    %             - (inv(sigma)*mu+theta)'*inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta)...
    %             +y'*inv(R)*y);

    % %         x1 = -0.5*(log(det(inv(sigma)))+log(det(Lambda))...
    % %             -log(det(inv(inv(sigma)+Lambda)))...
    % %             + mu'*inv(sigma)*mu - (A*X)'*inv(Q)*(A*X)...
    % %             + y'*R*y - KK'*inv(K)*KK + theta'*inv(Lambda)*theta...
    % %             - (inv(sigma)*mu+theta)'*inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta));

%                 x1 = 0.5*(-log(det(inv(sigma)))-log(det(Lambda))...
%                 +log(det(inv(inv(sigma)+Lambda)))...
%                 - mu'*inv(sigma)*mu - (A*X)'*inv(Q)*(A*X) - theta'*inv(Lambda)*theta...
%                 - y'*R*y + KK'*inv(K)*KK...
%                 + (inv(sigma)*mu+theta)'*inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta)...
%                 );
%                 x1 = 0.5*(-log(det(inv(sigma)))...
%                 +log(det(inv(inv(sigma)+Lambda)))...
%                 - mu'*inv(sigma)*mu - (A*X)'*inv(Q)*(A*X)...
%                 + KK'*inv(K)*KK...
%                 + (inv(sigma)*mu+theta)'*inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta)...
%                 );
%                 x1 = 0.5*(log(det(inv(sigma)))...
%                 -log(det((inv(sigma)+Lambda)))...
%                 - mu'*inv(sigma)*mu...
%                 + (inv(sigma)*mu+theta)'*inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta)...
%                 );
%simple
                x1 = 0.5*(-log(det((sigma)))...
                +log(det(inv(inv(sigma)+Lambda)))...
                - mu'*inv(sigma)*mu...
                + (inv(sigma)*mu+theta)'*inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta)...
                );
    %        end
       %%        0.5*(innov'*inv(S)*innov + dim*log(2*pi) + log(det(S)));   % 评分
            if j==1
            testX1Ctrl(:,i) = [x1;x1];
            testMixCtrl(:,i) = inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta);
            testObsCtrl(:,i) = inv(Lambda)*(theta);
            end
            if probIdx(j,i) ~= 0
                x1 = round(x1,5);
                 pMatJ1(j,i) = x1+log(probIdx(j,i));  
%                pMatJ1(j,i) = x1*(probIdx(j,i));  
            end
            sampleZWeight(j,i) = pMatJ1(j,i);       
        end
        pMax(i) = max(pMatJ1(1 : totalL + 1,i));
    end
    matJ1Max = max(pMax);
    for i = 1 : J
        for j = 1 : totalL + 1
             sampleZWeight(j,i) = sampleZWeight(j,i) - matJ1Max; 
             sampleZWeight(j,i) = exp(sampleZWeight(j,i));
             sampleZWeight(j,i) = round(sampleZWeight(j,i),5);
        end
    end
%    matJ1Max = max(pMax);   %% %%带“%%”去掉对结果影响不大，状态更多了
    pMatJ1 = sampleZWeight;%round(pMatJ1,3);
    for i = 1 : J
        totalL = thisTarg{i}{4};
%         for j = 1 : totalL + 1   %% %%带“%%”去掉对结果影响不大，状态更多了
%             pMatJ1(j,i) = (pMatJ1(j,i) - matJ1Max);  %%%%带“%%”去掉对结果影响不大，状态更多了
%             pMatJ1(j,i) = exp(pMatJ1(j,i)); %% %%带“%%”去掉对结果影响不大，状态更多了
%         end  %% %%带“%%”去掉对结果影响不大，状态更多了
        probC(i) = sum(pMatJ1(1:totalL + 1,i),1);
    end
    probC = cumsum(probC);
    for i = 1 : J
        u = rand * probC(J);
        %% j就是采样出的索引
        for jj = 1 : J
            if probC(jj) > u
                break;
            end
        end
        sampleWeight(jj) = sampleWeight(jj) + 1;
    end
    
	newTarg = cell(1, J);
    sumIndex = 1;
    tempZWeight = zeros(MaxStateNum + 1, J);
    for j = 1 : J
        for jj = 1 : sampleWeight(j)
            %这段代码完成粒子重采样
            newTarg{sumIndex} = thisTarg{j};
            %这段代码完成sampleZWeight的冲采样
            tempZWeight(:,sumIndex) = sampleZWeight(:,j);
            sumIndex = sumIndex + 1;
        end
    end
    sampleZWeight = tempZWeight;
end
function [errRMS, lose] = ShowCtrlCell(thisTarg, nStep)
global MaxStateNum
global dim
global J
%    ctrlArray = zeros(dim, MaxStateNum*J);
    ctrlArray = zeros(dim, J);
%     for i = 1 : J
%         ctrlArray(:,MaxStateNum*(i-1)+1:MaxStateNum*i) = thisTarg{i};
%     end
    for i = 1 : J
        ctrlArray(:,i) = thisTarg{i}{5}(:,nStep);
    end
    maxArray = (ctrlArray(1,:) > 0) .* sqrt(ctrlArray(1,:).^2)...
        + (ctrlArray(2,:) > 0) .* sqrt(ctrlArray(2,:).^2);
    maxCtrl = find(maxArray == max(maxArray));
    disp(maxCtrl);
end
function [errRMS, lose] = ShowCtrlArray(testMixCtrl)
%    ctrlArray = zeros(dim, MaxStateNum*J);
    ctrlArray = testMixCtrl;
    maxArray = (ctrlArray(1,:) > 0) .* sqrt(ctrlArray(1,:).^2)...
        + (ctrlArray(2,:) > 0) .* sqrt(ctrlArray(2,:).^2);
    maxCtrl = find(maxArray == max(maxArray));
    disp(maxCtrl);
end
