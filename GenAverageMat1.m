function mat1 = GenAverageMat1(thisTarg, Z, pd)
global J
global T

dim = 2;
FU = GenFU(T);
A = FU(1:6,1:6);
BU = GenBU(T);
B = BU(1:6,1:2);
C = [1 0 0 0 0 0; 0 0 0 1 0 0];
y = Z(1:2,1);
RCov = GenCov(Z);
R = RCov(1:2,1:2);
Q = 0.01*eye(2,2);

%% X(1) P(2) currentState(3) L(4) cellCtr(5)...
%% cellTrans(6) alpha(7) beta(8) gamma(9) M(10)...
%% auxPhi(11) auxTau(12) auxZeta(13)
%% kappa(14) vartheta(15)(dim=u��ά�ȣ�����ѡ2λ�� nu(16) Delta(17)��2*2����
%% --��һ�ж���NIW�Ĳ���
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(18) ÿ��״̬���ֵ�����
% innov = H*thisTarg{3} - Z(:, i); % meas innovation ��ÿһ�����⣬ȡ��Ŀ�꣬�����裬Ȼ������ѭ��
% S = H*thisTarg{4} * H' + R; % prior cov of innov
% dim = length(thisTarg{3});
% x1 = log((1-pd)/pd); 
% x2 = 0.5*(innov'*inv(S)*innov + dim*log(2*pi) + log(det(S)));   % ����
% mat1(i, j) = x1 + x2;
mat1 = 0;
matJ1 = zeros(J,1);
 
    %% �ȼ�����Բ����
    for i = 1 : J
       noGate = 1;
       matJ1(i,1) = CalcCost(thisTarg{3}{i}, y, noGate, 0);
       mat1 = mat1 + matJ1(i,1)/J;
       if mat1 == inf
           break;
       end
    end
    x1 = log((1-pd)/pd);
%     x1 = 0;
    mat1 = x1 + mat1;
end