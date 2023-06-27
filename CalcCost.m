%% 这个版本在预测步和更新步都用了复合估计
%% 参考 多机动目标跟踪算法与仿真平台开发 P64 (3) (4)
%% (3)里面的权重，用状态z的先验估计 probIdx
%% (4)为每个状态z计算预测状态的X和P
function matJ1 = CalcCost(thisTarg, Z, noGate, rawDate)
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%这里有个巨大的问题就是最后值赋予的权重是0.05还是0.5
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
global T
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

%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u的维度，这里选2位） nu(16) Delta(17)（2*2矩阵）
%% --这一行都是NIW的参数
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) 每个状态保持的数量

% kappa = thisTarg{14};
% vartheta = thisTarg{15};
% nu = thisTarg{16};
% Delta = thisTarg{17};
%% 计算比例
    totalL = thisTarg{4};
    p = zeros(totalL+1, 1);
    pSum = 0;
    probIdx = zeros(totalL+1, 1);
    currentS = thisTarg{3};
    for j = 1 : totalL + 1
        if j == totalL + 1
            p(j,1) = thisTarg{7} * thisTarg{8}(j);
        else
            p(j,1) = thisTarg{7} * thisTarg{8}(j)...
                + thisTarg{6}(currentS,j);
        end
        pSum = p(j,1) + pSum;
    end
    u = pSum;
    for j = 1 : totalL + 1
        probIdx(j,1) = p(j,1)/u;
    end
    %% 计算值
    X = thisTarg{1};
    P = thisTarg{2};
    XPredict = zeros(size(X,1),totalL+1);
    PPredict = zeros(size(P,1),size(P,2),totalL+1);
    XComposite = zeros(size(X));
    PComposite = zeros(size(P));
    if norm(y-C*X,2) > 1000 %2000 simple
        if rawDate == 1
            matJ1 = NaN;
            return;
        else
            matJ1 = inf;
            return;
        end
    end

    for j = 1 : totalL + 1
        if j == totalL + 1
            kappa = thisTarg{14};
            vartheta = thisTarg{15};
            nu = thisTarg{16};
            Delta = thisTarg{17};
            mu = vartheta;
            sigma = ((kappa + 1)*nu/(kappa*(nu - dim - 1)))*Delta;
        else
            kappa = thisTarg{14} + thisTarg{18}(j);
            vartheta = (thisTarg{14}*thisTarg{15} + thisTarg{5}(:,j))...
                /kappa;
            nu = thisTarg{16} + thisTarg{18}(j);
            Delta = (thisTarg{16}*thisTarg{17}...
                + [thisTarg{19}(j,1:2);thisTarg{19}(j,3:4)]...
                + thisTarg{14}*(thisTarg{15}*thisTarg{15}')...
                -kappa*(vartheta*vartheta'))/nu;
            mu = vartheta;
            sigma = ((kappa + 1)*nu/(kappa*(nu - dim - 1)))*Delta;
        end

        XPredict(:,j) = A * X + B * mu;
%        PPredict(:,:,j) = A * P * A' + B * sigma * B'+ Q;
        PPredict(:,:,j) = A * P * A' + Q;

        XComposite = XComposite + probIdx(j)*XPredict(:,j);
    end
    for j=1:totalL+1
        PComposite = PComposite + probIdx(j)*(PPredict(:,:,j)+(XComposite-XPredict(:,j))...
            *(XComposite-XPredict(:,j))');
    end
%     if norm(y-C*XComposite,2) > 600 %2000 simple
%         if rawDate == 1
%             matJ1 = NaN;
%             return;
%         else
%             matJ1 = inf;
%             return;
%         end
%     end

    innov = y - C*XComposite;
    SPredict = C*PComposite*C' + R;
    %simple
%    x2 = 0.05*(innov'*inv(SPredict)*innov + dim*log(2*pi) + log(det(SPredict)));   % 评分
    x2 = 0.0005*(innov'*inv(SPredict)*innov + dim*log(2*pi) + log(det(SPredict)));   % 评分
    matJ1 = x2;
end
