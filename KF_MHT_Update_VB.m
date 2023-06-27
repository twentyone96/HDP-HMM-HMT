function [newTarg,newState,ZetaCD,VarepsilonAB, Xi, Varphi]...
    = KF_MHT_Update_VB(newTarg,newState,...
    ZZ, ZetaCD, VarepsilonAB, Xi, Varphi, t)
global TL %truncation level
global TargetNum %Ŀ�����
global Q_coeff_VB
global T
global a_0
global tau_0

    FU = GenFU(T);
    A = FU(1:6,1:6);
    BU = GenBU(T);
    B = BU(1:6,1:2);
    C = [1 0 0 0 0 0; 0 0 0 1 0 0];
%     RCov = GenCov(Z);
%     R = RCov(1:2,1:2);
    Q = Q_coeff_VB*eye(6,6);
    
    a_0 = 1;
    tau_0 = 1;

%% X(1) P(2) currentState(3) cellCtr(4) ��mu...
%% cellTrans(5) ����ת�Ƹ��ʾ�����Sigma
%% kappa(6) vartheta(7)(dim=u��ά�ȣ�����ѡ2λ�� nu(8) Delta(9)��2*2����
%% --��һ�ж���NIW�Ĳ���
%% NIW(0.001, 0, 50, I_u(2*2)) 
%% cellStates(10) ÿ��״̬���ֵ�����
%% cellCtrAccum(11) ÿ��״̬U*U'
%% targPi(12) ����Ŀ��Ĺ۲�������,��Ӧ���ĵ�Pi_tm
IterNum = 2;
%TL = 10; %�����Ƿ�����룬��������Ч�ʣ�����Ϊ10���������˲��Ǳ�������25������һ��
numMeasure = size(ZZ,2); %�����ڹ۲�����
%U_E = cell(TargetNum,1); % U������
%U_P = cell(TargetNum,1); % U�ķ���
% 0.PMHTû�а취��MHT�������к�������ʼ���㣬����ȱ�㣬�Ǿ�Ĭ���Ѿ�֪��������
% ��ʼ�����ģ�����һ�����裻���û�й۲��뺽������������δ���
% �ڶ�������������PMHTû�а취����Ŀ�����������������������Ŀ��

% 1. Ϊÿһ��Ŀ�����۲��������gate����������Ϊ500m������CalcCost1.m
% 2. Ϊÿһ��Ŀ�����w_mtl,l�Ǹ�Ŀ��۲�ļ���
% 3. ����p(z_j)��ֵ
% 4. ����p(u|z_tm=i)��ֵ
% 5. ����Xi
% 6. ����Varphi

% 1. Ϊÿ��Ŀ���gate�ڵĹ۲ⵥ����ȡ����,��������һ���̶�ֵ1000����hdp-HMM�ķ���һ��
measGate = zeros(numMeasure,TargetNum); %��measGate����¼�ĸ��۲������ĸ�Ŀ��
for i = 1 : numMeasure
    y = ZZ(1:2,i);
    for j = 1 : TargetNum
%        if norm(y-C*XComposite(:,j),2) < 500
         if norm(y-C*newTarg{j}{1},2) < 1000
            measGate(i,j) = 1;
%             break;
        end
    end
end
%������Ҫ�и������ǽ����е�û�й���Ĺ۲���������һ��
actualMeasure = zeros(numMeasure,1);
for i = 1 : numMeasure
    if sum(measGate(i,:)) ~= 0
        actualMeasure(i) = 1;
    end
end
actualZ = zeros(3,sum(actualMeasure));
actualMeasGate = zeros(sum(actualMeasure), TargetNum);
j = 1;
for i = 1 : numMeasure
    if actualMeasure(i) == 1
        actualZ(:,j) = ZZ(:,i);
        actualMeasGate(j,:) = measGate(i,:);
        j = j + 1;
    end
end
Z = actualZ;
measGate = actualMeasGate;
numMeasure = sum(actualMeasure);


   XPredict = zeros(size(newTarg{1}{1},1),TL,TargetNum);
   PPredict = zeros(size(newTarg{1}{2},1),size(newTarg{1}{2},2),TL,TargetNum);
    XComposite = zeros(size(newTarg{1}{1},1),TargetNum);
    PComposite = zeros(size(newTarg{1}{2},1),size(newTarg{1}{2},2),TargetNum);
for i_iter = 1 : IterNum
   % ����Ҫ����X_t|t-1��P_t|t-1��Ҫ�ø�����⣬����ֱ����kalman�˲�
    for targetIndex = 1 : TargetNum
       X = newTarg{targetIndex}{1};
       P = newTarg{targetIndex}{2};
       for j = 1 : TL
           XPredict(:,j,targetIndex) = A * X + B * newState{j}{1};
           PPredict(:,:,j,targetIndex) = A * P * A' + B * newState{j}{2} * B' + Q;
           XComposite(:,targetIndex) = XComposite(:,targetIndex) ...
               + newTarg{targetIndex}{16}(j)*XPredict(:,j,targetIndex);
       end
       for j = 1 : TL
           PComposite(:,:,targetIndex) = PComposite(:,:,targetIndex) ...
               + newTarg{targetIndex}{16}(j)*(PPredict(:,:,j,targetIndex)...
               +(XComposite(:,targetIndex)-XPredict(:,j,targetIndex))...
               *(XComposite(:,targetIndex)-XPredict(:,j,targetIndex))');
       end
    end
   % ���numMeasure == 1,ֻ��һ��Ŀ���й۲⣬Ӧ��Ҳû�����⣬ֻ��sumMu��sumSigma�ļ�����õ�
   % wTML,���wTML==0����ʵ��������ԭ����ֵû�з����仯�������⣬�����������Ͳ�����
   % ���Ըɴ������������ĳ� numMeasure == 0���� numMeasure == 0 || numMeasure == 1
   % ��������и��õĴ�ʩ������������ÿһ�������ļ���ʱ���ж�numMeasure == 1
   if numMeasure == 0 || numMeasure == 1
        for targetIndex = 1 : TargetNum
           newTarg{targetIndex}{1} = XComposite(:,targetIndex);%XP;%XComposite(:,targetIndex);
           newTarg{targetIndex}{2} = PComposite(:,:,targetIndex);%PP;%PComposite(:,:,targetIndex);
        end
        break;
   end
   % 2.����w_tml,��Ҫ֪��ÿ��Ŀ���������ڵĹ۲�����
   % �����м����
   temp6 = zeros(numMeasure,TargetNum);
%   temp7 = zeros(numMeasure,1);
   % ��ʼ��Ŀ��-�۲�ƥ����ʾ��󣬶�Ӧw_tml
   wTML = zeros(numMeasure,TargetNum);
   %Ҫ���ǹ۲�����Ϊ0,1��3�ȼ��������������
   for targetIndex = 1 : TargetNum
       for measIndex = 1 : numMeasure
           inGate = measGate(measIndex, targetIndex);
           if inGate == 1
%               for tlIndex = 1 : TL
                   y = Z(1:2,measIndex);
                   RCov = GenCov(y);
                   R = RCov(1:2,1:2);
                   temp4 = ComputeLog(y,R,...
                       C*XComposite(:,targetIndex),...
                       [PComposite(1,1,targetIndex),PComposite(1,4,targetIndex);...
                       PComposite(4,1,targetIndex),PComposite(4,4,targetIndex)]);
%                   temp5 = temp5 + newTarg{targetIndex}{16}(tlIndex) * temp4;
%               end
               temp4 = temp4 + log(newTarg{targetIndex}{12});
               temp6(measIndex,targetIndex) = exp(temp4);
           end
%           temp7(measIndex) = temp7(measIndex) + temp6(measIndex,targetIndex);
       end   
   end
   % ����w_tml,����ÿ���۲����
   for targetIndex = 1 : TargetNum
       for measIndex = 1 : numMeasure
           if sum(temp6(measIndex,:)) == 0
               for i = 1 : TargetNum
                   if measGate(measIndex,i) == 1 && sum(temp6(:,i)) == 0
                       wTML(measIndex,i) = 1;
                       break;
                   end
               end
%               wTML(measIndex,targetIndex) = 0;
           else
               wTML(measIndex,targetIndex) = temp6(measIndex,targetIndex)...
                   /sum(temp6(measIndex,:));
           end
       end
   end
   % ȥ������ÿһ��Ŀ�궼��0����
   emptyIndex = [];
   for measIndex = 1 : numMeasure
       if wTML(measIndex,1) == 0 && wTML(measIndex,2) == 0
           emptyIndex = [emptyIndex;measIndex];
       end
   end
   for measIndex = 1 : size(emptyIndex)
       deleteIndex = emptyIndex(measIndex);
%        wTML(deleteIndex,:) = [];
%        emptyIndex = emptyIndex - 1;
%        numMeasure = numMeasure - 1;
%        Z(:,deleteIndex) = [];
   end
   % ���ر�С��ֵȥ��
   if isempty(wTML) ~= 1
       for measIndex = 1 : numMeasure
           for targetIndex = 1 : TargetNum
               if wTML(measIndex,targetIndex) < 0.01 && wTML(measIndex,targetIndex) > 0
                    wTML(measIndex,targetIndex) = 0;
                    % �����и�Сtrick��ֻ�е�target==2���������
                    wTML(measIndex,TargetNum+1-targetIndex) = 1; 
               end
           end
       end
   end
   %����Pi
   for targetIndex = 1 : TargetNum
       newTarg{targetIndex}{12} = sum(wTML(:,targetIndex))/sum(sum(wTML));
   end
   % 3. ����p(x_tm)��ֵ
   % �����sumMu��sumSigma�൱��X_t|t��P_t|t�Ĺ���ֵ
   % ����ɼ����newTarg{1/2}{1}������һ�ŵĵĹ��ƣ���X_t-1|t-1��P_t-1|t-1
   sumMu = zeros(size(newTarg{1}{1},1),TargetNum);
   sumSigma = zeros(size(newTarg{1}{2},1),size(newTarg{1}{2},2),TargetNum);
   for targetIndex = 1 : TargetNum
       % �ȼ��� y�ǲ��ֵ�x��ֵ
%        for measIndex = 1 : numMeasure
%            inGate = measGate(measIndex, targetIndex);
%            if inGate == 1
%                y = Z(1:2,measIndex);
%                RCov = GenCov(y);
%                R = RCov(1:2,1:2);
%                sumSigma(:,:,targetIndex) = sumSigma(:,:,targetIndex)...
%                    + wTML(measIndex,targetIndex)*C'/(R)*C;
%                sumMu(:,targetIndex) = sumMu(:,targetIndex)...
%                    + wTML(measIndex,targetIndex)*C'/(R)*y;
%            end
%        end
       % ������һ�εļ��㷨
       meanZ = zeros(size(Z(1:2,1)));
       meanR = zeros(2,2);
       for measIndex = 1 : numMeasure
           inGate = measGate(measIndex, targetIndex);
           if inGate == 1
               y = Z(1:2,measIndex);
               RCov = GenCov(y);
               R = RCov(1:2,1:2);
               meanZ = meanZ + (wTML(measIndex,targetIndex)*y)...
                   /sum(wTML(:,targetIndex));
               meanR = meanR + (wTML(measIndex,targetIndex)*R)...
                   /sum(wTML(:,targetIndex));
           end
       end
       % wTML�����[1 0;1 0],����һ��Ŀ��û�ж�Ӧ�۲�����������sum(wTML) == 0
        [LambdaChol,err] = cholcov(meanR,0);
        if err ~= 0
            error(message('stats:mvnpdf:BadMatrixSigma'));
        end
       sumSigma(:,:,targetIndex) = sumSigma(:,:,targetIndex)...
           + C'/(meanR)*C;
       sumMu(:,targetIndex) = sumMu(:,targetIndex)...
           + C'/(meanR)*meanZ;       
       % �ټ���(Xtm|X_t-1,m,U(z_tm))^z_tm,j
        for tlIndex = 1 : TL
           sumSigma(:,:,targetIndex) = sumSigma(:,:,targetIndex)...
               + newTarg{targetIndex}{16}(tlIndex)*inv(PPredict(:,:,tlIndex,targetIndex));
           sumMu(:,targetIndex) = sumMu(:,targetIndex)...
               + newTarg{targetIndex}{16}(tlIndex)*inv(PPredict(:,:,tlIndex,targetIndex))*XPredict(:,tlIndex,targetIndex);   
%            sumSigma(:,:,targetIndex) = sumSigma(:,:,targetIndex)...
%                + newTarg{targetIndex}{16}(tlIndex)*inv(Q);
%            sumMu(:,targetIndex) = sumMu(:,targetIndex)...
%                + newTarg{targetIndex}{16}(tlIndex)*inv(Q)...
%                * (A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1});   
        end
       % �����൱�����˼�
%        sumSigma(:,:,targetIndex) = sumSigma(:,:,targetIndex)...
%            + inv(PComposite(:,:,targetIndex));
%        sumMu(:,targetIndex) = sumMu(:,targetIndex)...
%            + (PComposite(:,:,targetIndex))\XComposite(:,targetIndex);           
%         [LambdaChol,err] = cholcov(sumSigma(:,:,targetIndex));
%         if err ~= 0
%             error(message('stats:mvnpdf:BadMatrixSigma'));
%         end
       sumSigma(:,:,targetIndex) = inv(sumSigma(:,:,targetIndex)); 
       sumMu(:,targetIndex) = sumSigma(:,:,targetIndex)...
           * sumMu(:,targetIndex);
   end
   % 4. ����p(z_j)��ֵ
   sumZ = zeros(TL,TargetNum);
   for targetIndex = 1 : TargetNum
       for tlIndex = 1 : TL
           % ����log��N(newStateU;Lambda\theta,newStateP)����ֵ
           % non-log
%           temp4 = ComputeNoLog(sumMu(:,targetIndex),PPredict(:,:,tlIndex,targetIndex),...
%               XPredict(:,tlIndex,targetIndex),...
%               A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B');
           % log
           temp4 = ComputeLog(sumMu(:,targetIndex),A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B'+Q,...
               A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1},...
               sumSigma(:,:,targetIndex));
%               sumSigma(:,:,targetIndex));
           %_1 ~ Q=5000
%            temp4 = ComputeLog(sumMu(:,targetIndex),Q,...
%                A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1},...
%                sumSigma(:,:,targetIndex)+A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B');
           sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex) + temp4;                    
           sumTemp = 0;
           for lastIndex = 1 : TL
               %�ڼ���% ElnA_1
%                if newTarg{targetIndex}{17}(lastIndex) ~= 0
%                    % log
%                    tempA = computeA(lastIndex,tlIndex,Varphi,VarepsilonAB);
%                    if tempA ~= 0
%                         sumTemp = sumTemp + log(newTarg{targetIndex}{17}(lastIndex))...
%                             + tempA;
%                    end
%                end
               % no-log_1 20210917
%                  sumTemp = sumTemp + (newTarg{targetIndex}{17}(lastIndex)...
%                         * computeANoLog(lastIndex,tlIndex,Varphi,VarepsilonAB));
               % 20210917
                 sumTemp = sumTemp + (newTarg{targetIndex}{17}(lastIndex)...
                        * computeA(lastIndex,tlIndex,Varphi,VarepsilonAB));
           end
           % log_1
%            sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex) + (sumTemp);
           % no-log_1 20210917
%           sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex) + log(sumTemp);
           % 20210917
           sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex) + (sumTemp);
           % z����E_q(u)log��N(newStateU,newStateP)��
%            sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex)...
%               -log(2*pi*sqrt(det(newState{tlIndex}{2})))...
%               - 0.5*trace(newState{tlIndex}{2});
       end
       % ���������ζ�ţ��κ�һ�ֿ��Ʒ����Ӧ�Ĺ��ƾ���۲ⶼ̫Զ����Ȼ��Ϊ0���Ǿ͸���ʼֵ
       if sum(sumZ(:,targetIndex)) == 0
           sumZ(:,targetIndex) = 1;
       end
   end
   % ��sumZ���й�һ��
   % log
   for targetIndex = 1 : TargetNum
       maxZ = max(sumZ(:,targetIndex));
       sumZ(:,targetIndex) = sumZ(:,targetIndex) - maxZ; 
   end
   sumZ = exp(sumZ);
   %%%%
   sumZ = round(sumZ,5);
   % ����ÿ��״̬���ʵ���Сֵ,û���ã��������⻹��Q
%    for targetIndex = 1 : TargetNum
%       for tlIndex = 1 : TL
%           if sumZ(tlIndex,targetIndex) < 0.0001
%             sumZ(tlIndex,targetIndex) = 0.0001;
%           end
%       end
%    end   
   for targetIndex = 1 : TargetNum
      maxZ = sum(sumZ(:,targetIndex));
      for tlIndex = 1 : TL
          sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex)/maxZ;
      end
      %����ȥ����Сֵ
%       [maxValue,maxIndex] = max(sumZ(:,targetIndex));
%       for tlIndex = 1 : TL
%           if sumZ(tlIndex,targetIndex) < 0.001 && sumZ(tlIndex,targetIndex) > 0
%             sumZ(maxIndex,targetIndex) = sumZ(maxIndex,targetIndex) + sumZ(tlIndex,targetIndex);
%             %Ȼ�����Ϊ0��Ԫ�ط���0.0001
%             sumZ(tlIndex,targetIndex) = 0.00001;
%             sumZ(maxIndex,targetIndex) = sumZ(maxIndex,targetIndex) - sumZ(tlIndex,targetIndex);
%           end
%       end
      newTarg{targetIndex}{16} = sumZ(:,targetIndex);
%       if i_iter == IterNum
%           newTarg{targetIndex}{17} = sumZ(:,targetIndex);
%       end
   end
   % 4.����p(u|z_tm=i)  
   sumLambda = zeros(2,2,TL);
   sumTheta = zeros(2,TL);
   % ���ݹ�ʽ������һ��ȷʵֻ���۲�ֵ�й�
   for tlIndex = 1 : TL
       for targetIndex = 1 : TargetNum
           sumLambda(:,:,tlIndex) = sumLambda(:,:,tlIndex)...
               + newTarg{targetIndex}{16}(tlIndex)*B'/Q*B;
           sumTheta(:,tlIndex) = sumTheta(:,tlIndex)....
               + newTarg{targetIndex}{16}(tlIndex)*B'/Q*(sumMu(:,targetIndex)-A*newTarg{targetIndex}{1});
       end
   end
   for tlIndex = 1 : TL
       tempS = newState{tlIndex}{2};
       tempM = newState{tlIndex}{1};
        [LambdaChol,err] = cholcov(sumLambda(:,:,tlIndex) + inv(tempS),0);
        if err ~= 0
            error(message('stats:mvnpdf:BadMatrixSigma'));
        end
       newState{tlIndex}{2} = inv(sumLambda(:,:,tlIndex) + inv(tempS));
       newState{tlIndex}{1} = newState{tlIndex}{2}*(sumTheta(:,tlIndex) + (tempS)\tempM);
   end
    % 5. ����\Xi,��״̬��Ǩ���ʾ��󣬶�ӦZ
%   Xi = zeros(TL,TL,TargetNum); % һTL���ϸ�ʱ��״̬����TL���ϸ�ʱ�̵�stick
   sumXi = zeros(TL,TL,TargetNum);
   % ��ʼ��Ŀ��-�۲�ƥ����ʾ��󣬶�Ӧw_tml
   for targetIndex = 1 : TargetNum
       for i = 1 : TL
           for j = 1 : TL
               tempXi = 0;
               for tlIndex = 1 : TL % tlIndx�ǵ�ǰ��״̬
                   temp4 = ComputeLog(sumMu(:,targetIndex),A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B'+Q,...
                       A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1},...
                       sumSigma(:,:,targetIndex));
%                    temp4 = ComputeLog(sumMu(:,targetIndex),Q,...
%                        A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1},...
%                        sumSigma(:,:,targetIndex)+A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B');
                   tempXi = tempXi + Varphi(i,j,tlIndex)*newTarg{targetIndex}{16}(tlIndex)*temp4;
               end
               sumXi(i,j,targetIndex) = sumXi(i,j,targetIndex) + tempXi + ...
                   getZeta2(i,j,VarepsilonAB);
               %���==0 ,��log(0)=-inf,exp(-inf)=0;
               %���Xi ��Ԫ��=0������ʲô���
               if newTarg{targetIndex}{17}(i) ~= 0
                   sumXi(i,j,targetIndex) = sumXi(i,j,targetIndex) + ...
                       log(newTarg{targetIndex}{17}(i));
               else
                   sumXi(i,j,targetIndex) = 0;
                   break;
               end
           end    
       end   
   end
   % ת����ָ��
   for targetIndex = 1 : TargetNum
       for i = 1 : TL
           tempSum = sum(sumXi(i,:,targetIndex));
           for j = 1 : TL
               if tempSum == 0
                   sumXi(i,j,targetIndex) = -inf;%1/TL;
               else
                   sumXi(i,j,targetIndex) = sumXi(i,j,targetIndex)/tempSum;
               end
           end
       end
%       sumXi(:,:,targetIndex) = sumXi(:,:,targetIndex) - max(max(sumXi(:,:,targetIndex))); 
   end   

   sumXi = exp(sumXi);
   %��һ��
   Xi = normalVarepsilon(sumXi);
% 
   % 6. ����\Varphi����ӦC
%   Varphi = zeros(TL,TL,TL); %һTL���ϸ�ʱ��״̬����TL���ϸ�ʱ�̵�stick,���ǵ�ǰʱ��״̬
   sumVarphi = zeros(TL,TL,TL);
   % ��ʼ��Ŀ��-�۲�ƥ����ʾ��󣬶�Ӧw_tml
   for tlIndex = 1 : TL % tlIndx�ǵ�ǰ��״̬
       for i = 1 : TL
           for j = 1 : TL
               tempSumW = 0;
               for targetIndex = 1 : TargetNum
                 temp4 = ComputeLog(sumMu(:,targetIndex),A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B'+Q,...
                    A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1},...
                    sumSigma(:,:,targetIndex));
%                   temp4 = ComputeLog(sumMu(:,targetIndex),Q,...
%                        A*newTarg{targetIndex}{1}+B*newState{tlIndex}{1},...
%                        sumSigma(:,:,targetIndex)+A*newTarg{targetIndex}{2}*A'+B*newState{tlIndex}{2}*B');
                   tempSumW = tempSumW + Xi(i,j,targetIndex)*newTarg{targetIndex}{16}(tlIndex)*temp4;
               end
               sumVarphi(i,j,tlIndex) = sumVarphi(i,j,tlIndex) + getZeta1(tlIndex,ZetaCD)  + tempSumW;
           end
       end
  end   
   % ת����ָ��
   for i = 1 : TL
       for j = 1 : TL
           tempSum = sum(sumVarphi(i,j,:));
           for tlIndex = 1 : TL % tlIndx�ǵ�ǰ��״̬
               if tempSum == 0
                   sumVarphi(i,j,tlIndex) = -inf;%1/TL;
               else
                   sumVarphi(i,j,tlIndex) = sumVarphi(i,j,tlIndex)/tempSum;
               end
           end
       end
   end     
   sumVarphi = exp(sumVarphi);
   %��һ��
   Varphi = normalVarphi(sumVarphi);
    
    % 7. ���²���a,b,c,d
    %��ʼ��\Zeta, \Zeta_c, \Zeta_d  ��Ӧ���������Beta
    for i = 1 : TL
        zetaVarphi = getVarphi(Varphi,i);
        ZetaCD(i,1) = 1 + zetaVarphi(1);
        ZetaCD(i,2) = a_0 + zetaVarphi(2);
    end
    %��ʼ��\Varepsilon, \Varepsilon_a, \Varepsilon_b ��Ӧ���������Pi
    for i = 1 : TL
        for ii = 1 : TL
            epsilonVarphi = getXi(Xi,i,ii);
            VarepsilonAB(i,ii,1) = 1 + epsilonVarphi(1);
            VarepsilonAB(i,ii,2) = tau_0 + epsilonVarphi(2);
        end
    end
    
    %8. ����predictive log likelihood
    %���ȼ���int(qlnp)
    qlnp = 0;
%     % W_mtl
%     for targetIndex = 1 : TargetNum
%        for measIndex = 1 : numMeasure
%            temp5 = 0;
%            inGate = measGate(measIndex, targetIndex);
%            if inGate == 1
%                for tlIndex = 1 : TL
%                    y = Z(1:2,measIndex);
%                    RCov = GenCov(y);
%                    R = RCov(1:2,1:2);
%                    [Lambda, theta] = ComputeLandT(newTarg{targetIndex}{1},newTarg{targetIndex}{2},...
%                        y,Q,C,R,B,A);
%                    % ����log��N(newStateU;Lambda\theta,newStateP)����ֵ
%                    temp4 = ComputeELikelihood(Lambda,theta,...
%                        newState{tlIndex}{1},newState{tlIndex}{2});
%                    temp5 = temp5 + newTarg{targetIndex}{16}(tlIndex) * temp4;
%                end
%                temp5 = temp5 + log(newTarg{targetIndex}{12});
%            end
%            qlnp = qlnp + temp5;
%        end   
%     end
%     % Z
%     for targetIndex = 1 : TargetNum
%        for tlIndex = 1 : TL
%            for measIndex = 1 : numMeasure
%                inGate = measGate(measIndex, targetIndex);
%                if inGate == 1
%                    y = Z(1:2,measIndex);
%                    RCov = GenCov(y);
%                    R = RCov(1:2,1:2);
%                    [Lambda, theta] = ComputeLandT(newTarg{targetIndex}{1},newTarg{targetIndex}{2},...
%                        y,Q,C,R,B,A);
%                    % ������
%                     [LambdaChol,err] = cholcov(Lambda,0);
%                     if err ~= 0
%                         error(message('stats:mvnpdf:BadMatrixSigma'));
%                     end
%                    % ����log��N(newStateU;Lambda\theta,newStateP)����ֵ
%                    temp4 = ComputeELikelihood(Lambda,theta,...
%                        newState{tlIndex}{1},newState{tlIndex}{2});
%                    sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex) + ...
%                        wTML(measIndex,targetIndex)*temp4;                    
%                end
%            end
%     %           sumZ(tlIndex,targetIndex) = exp(sumZ(tlIndex,targetIndex));
%            sumTemp = 0;
%            for lastIndex = 1 : TL
%                %�ڼ���% ElnA
%     %                sumTemp = sumTemp + log(newTarg{targetIndex}{17}(i))...
%     %                    + computeA(lastIndex,tlIndex,Varphi,VarepsilonAB);
%                sumTemp = sumTemp + (newTarg{targetIndex}{17}(i))...
%                    * computeA(lastIndex,tlIndex,Varphi,VarepsilonAB);
%            end
%            sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex) + (sumTemp);
%            % z����E_q(u)log��N(newStateU,newStateP)��
%     %            sumZ(tlIndex,targetIndex) = sumZ(tlIndex,targetIndex)...
%     %               -log(2*pi*sqrt(det(newState{tlIndex}{2})))...
%     %               - 0.5*trace(newState{tlIndex}{2});
%        end
%     end
    
    
    %��μ���int(qlnq)
    qlnq = 0;
    %������int(qlnp) - int(qlnq)

    if i_iter == IterNum
        for targetIndex = 1 : TargetNum
           newTarg{targetIndex}{17} = newTarg{targetIndex}{16};
           newTarg{targetIndex}{1} = sumMu(:,targetIndex);%XP;%XComposite(:,targetIndex);
           newTarg{targetIndex}{2} = sumSigma(:,:,targetIndex);%PP;%PComposite(:,:,targetIndex);
        end
    end
 
end
 
end

function [Lambda, theta] = ComputeLandT(X,P,y,Q,C,R,B,A)

%     X = A*X+B*mu;
%     P = A*P*A' + Q;
%     innov = y - C*X;
%     SPredict = C*P*C' + R;
%     G = P*C'*inv(SPredict);
%     X = X + G*innov;
%     P = (eye(6,6) - G*C)*P;

    K = inv(Q) + C'/(R)*C;
    S = inv(Q) - inv(Q)*inv(K)*inv(Q);
%���Ǽ򻯵��㷨
%    LambdaTmp11 = B'*S*B;
%    LambdaTmp12 = B'*S*A;
%    KK = inv(Q)*A*X + C'*inv(R)*y;
% 
%    thetaTmp1 = B'*inv(Q)*inv(K)*C'*inv(R)*y;
% 
%    Lambda = LambdaTmp11;
%    theta = thetaTmp1 - LambdaTmp12*X;

     Lambda = B'*S*B-B'*S*A/(A'*S*A+inv(P))*A'*S*B;
     theta = (B'/(Q))/(K)*C'/(R)*y...
         -(B'*S*A/(A'*S*A+inv(P)))...
         *((P)\X+(A'/(Q))/(K)*C'/(R)*y);
end
% �����˹�ֲ���logֵ
function EValue = ComputeLog(y,R,newStateU,newStateP)
    temp1 = (y - newStateU);
    temp2 = temp1*temp1';
    temp3 = temp2 + newStateP;
    EValue = -log(2*pi*sqrt(det(R))) - 0.5*trace(temp3/(R));
%    EValue = -0.5*trace(temp3/(R));
end
function EValue = ComputeLogNoConst(y,R,newStateU,newStateP)
    temp1 = (y - newStateU);
    temp2 = temp1*temp1';
    temp3 = temp2 + newStateP;
    EValue = - 0.5*trace(temp3/(R));
end
function EValue = ComputeNoLog(y,R,newStateU,newStateP)
    temp1 = (y - newStateU);
    temp2 = temp1*temp1';
    temp3 = temp2 + newStateP;
    EValue = -(2*pi*sqrt(det(R))) * exp(-0.5*trace(temp3/(R)));
end
% ��ȡ�ڶ���G��eta, eta_ij
function eta = getZeta2(i,j,Varepsilon)
% EqlnETA_ij = EqlnVarepsilon_ij + sum(1:j-1)Eqln(1-Varepsilon_ij)
    eta = 0;
    for index = 1 : j - 1
        eta = eta ...
            + psi(Varepsilon(i,index,2)) ...
            - psi(Varepsilon(i,index,1) + Varepsilon(i,index,2));
    end
    eta = psi(Varepsilon(i,j,1)) - psi(Varepsilon(i,j,1) + Varepsilon(i,j,2))...
        + eta;
end
function eta = getZeta1(j,Zeta)
    eta = 0;
    for index = 1 : j - 1
        eta = eta ...
            + psi(Zeta(index,2)) ...
            - psi(Zeta(index,1) + Zeta(index,2));
    end
    eta = psi(Zeta(j,1)) - psi(Zeta(j,1) + Zeta(j,2))...
        + eta;
end
%��һ��Varepsilon
function Xi = normalVarepsilon(sumXi)
global TargetNum
global TL
   Xi = zeros(TL,TL,TargetNum); % һTL���ϸ�ʱ��״̬����TL���ϸ�ʱ�̵�stick
   for targetIndex = 1 : TargetNum
       for i = 1 : TL
           tempSum = sum(sumXi(i,:,targetIndex));
           for j = 1 : TL
               if tempSum == 0
                   Xi(i,j,targetIndex) = 0;%1/TL;
               else
                   Xi(i,j,targetIndex) = sumXi(i,j,targetIndex)/tempSum;
               end
           end
       end
   end
end
function Varphi = normalVarphi(sumVarphi)
global TL
   Varphi = zeros(TL,TL,TL); % һTL���ϸ�ʱ��״̬����TL���ϸ�ʱ�̵�stick,���ǵ�ǰʱ��״̬
   for i = 1 : TL
       for j = 1 : TL
           tempSum = sum(sumVarphi(i,j,:));
           for tlIndex = 1 : TL % tlIndx�ǵ�ǰ��״̬
               if tempSum == 0
                   Varphi(i,j,tlIndex) = 0;%1/TL;
               else
                   Varphi(i,j,tlIndex) = sumVarphi(i,j,tlIndex)/tempSum;
               end
           end
       end
   end     
end
% ElnA
function A = computeA(lastIndex,tlIndex,Varphi,Varepsilon)
global TL
    tempA = 0;
    for sIndex = 1 : TL %stick��״̬��ǰ׺sΪstick
        tempA = tempA + Varphi(lastIndex,sIndex,tlIndex)...
            * exp(getZeta2(lastIndex,sIndex,Varepsilon));
    end
    if tempA == 0
        A = 0;
    else
        A = log(tempA);
    end
end
function A = computeANoLog(lastIndex,tlIndex,Varphi,Varepsilon)
global TL
    tempA = 0;
    for sIndex = 1 : TL %stick��״̬��ǰ׺sΪstick
        tempA = tempA + Varphi(lastIndex,sIndex,tlIndex)...
            * exp(getZeta2(lastIndex,sIndex,Varepsilon));
    end
    A = (tempA);
end
% function [Mu, Sigma] = ComputeMandS(X,P,y,Q,C,R,B,A,jjj)
%     global dim
%     kappaGlobal = 0.001;%0.001;
%     varthetaGlobal = zeros(dim,1); %[141,141]';%
%     nuGlobal = 50; 
%     DeltaGlobal = eye(dim);
% 
%     kappa = kappaGlobal + thisTarg{18}(jjj);
%     vartheta = (kappaGlobal*varthetaGlobal + thisTarg{5}(:,jjj))...
%         /kappa;
%     nu = nuGlobal + thisTarg{18}(jjj);
%     Delta = (nuGlobal*DeltaGlobal...
%         + [thisTarg{19}(jjj,1:2);thisTarg{19}(jjj,3:4)]...
%         + kappaGlobal*(varthetaGlobal*varthetaGlobal')...
%         -kappa*(vartheta*vartheta'))/nu;
%     mu = vartheta;
%     sigma = ((kappa + 1)*nu/((kappa*(nu - dim - 1))))*Delta;
% 
%     Mu = inv(inv(sigma)+Lambda)*(inv(sigma)*mu+theta);
%     Sigma = inv(inv(sigma)+Lambda);   
%     
%     muSample = mvnrnd(mu', sigma, 1)';
%     cellCtr(:,jjj) = cellCtr(:,jjj) + muSample;
%     cellStates(jjj) = cellStates(jjj) + 1;
%     cellCtrAccumTmp = muSample*muSample';
%     cellCtrAccum(jjj,:) = cellCtrAccum(jjj,:)...
%         + [cellCtrAccumTmp(1,:),cellCtrAccumTmp(2,:)];
% 
%     thisTarg{5} = cellCtr;
%     cellTrans(currentS,jjj) = cellTrans(currentS,jjj) + 1;
%     thisTarg{6} = cellTrans;
%     thisTarg{18} = cellStates;
%     thisTarg{19} = cellCtrAccum;
% end
