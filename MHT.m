function [cellEstm,RMS1,RMS2,RMS3,RMS4,RMS5,RMS6,RMSE1,RMSE2,RMSE3,RMSE4,RMSE5,RMSE6,...
    plotState1,plotState2,plotState3,plotState4,plotState5,plotState6,...
    stateEstmIMMRun,stateEstmRun,stateEstmVBRun] = MHT(state, meas, T, nStep, nTarg, ...
    densClt, densNew, Pd, M, N,F_CV, Q_CV, F_CA, Q_CA,...
    RMS1,RMS2,RMS3,RMS4,RMS5,RMS6,RMSE1,RMSE2,RMSE3,RMSE4,RMSE5,RMSE6,...
    iCycle,plotState1,plotState2,plotState3,plotState4,plotState5,plotState6)

global J
global MaxStateNum
global maxMes
global PTM
global IMMNUM
global IterNum              %iteration number for VB algorithm


global a_0
global tau_0

global TargetNum
global TL
VB_include = 1;
IMM_include = 1;
P_include = 1;
TargetNum = 2;
% VB算法状态个数

dim = 2;
%% initialization
%以下四个参数都通过参数传递进来,R由GenCov计算
% F = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
%H = [1 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0;0 0 0 0 0 0 0 0 0];
H = [1 0 0 0 0 0; 0 0 0 1 0 0];
P = 0; % zero initial covariance
% Q = q*[T^3/3 T^2/2 0 0; T^2/2 T 0 0; 0 0 T^3/3 T^2/2; 0 0 T^2/2 T];
R = diag([1^2, 1^2]);      %GenCov
%HDP-HMM-PMHT
Q = [];
% M; number of most probable hypos listed in each scan
% N; height of the hypo tree. It seems scan depth = N-1
maxNumOfHypoSeed = M^(N-1);  %16
maxNumOfHypo = M^N;	%64 % max number of hypoes kept at the same time
maxNumOfHypoAccum = 0;
for i = 1 : N - 1
    maxNumOfHypoAccum = maxNumOfHypoAccum + M^i;
end	% 1+4+16 = 21

%evolution state 初始化假设、航迹（目标）等
head = 1; %初始时刻有1个假设
rear = 1;
head_IMM = 1; %初始时刻有1个假设
rear_IMM = 1;

%cellHypo = cell(1,maxNumOfHypoAccum+maxNumOfHypo); %{{[]}}; %memory allocation
cellHypo = {{[]}};
cellHypo_IMM = {{[]}};
asso = size(maxMes,1);
prob = zeros(1,1);
hypo = {asso, prob};
%cellfun(@(v) v = hypo, cellHypo(1:maxNumOfHypoAccum+maxNumOfHypo));
%cellHypo{1,1:maxNumOfHypoAccum+maxNumOfHypo} = {hypo};

%cellTarg = cell(1,maxNumOfHypoAccum+maxNumOfHypo); %{{[]}}; %memory allocation
cellTarg = {{[]}};
cellTarg_IMM = {{[]}};
idx = zeros(1,1);
liftPoint = zeros(1,1);
targ = cell(1,19);
targ{1} = zeros(1,6);
targ{2} = eye(6,6);
targ{3} = zeros(1,1,'int8');
targ{4} = zeros(1,1,'int8');
targ{5} = zeros(2,MaxStateNum);
targ{6} = zeros(MaxStateNum,MaxStateNum);
targ{7} = zeros(1,1,'int8');
targ{8} = zeros(1,1+MaxStateNum,'int8');
targ{9} = zeros(1,1);
targ{10} = zeros(MaxStateNum,MaxStateNum,'int8'); 
targ{11} = zeros(1,1);
targ{12} = zeros(1,1);
targ{13} = zeros(MaxStateNum,1);
targ{14} = zeros(1,1);
targ{15} = zeros(dim,1);
targ{16} = zeros(1,1);
targ{17} = zeros(dim,dim);
targ{18} = zeros(MaxStateNum,1);
targ{19} = zeros(MaxStateNum,4);
%targs = cell(1,100);
%targs{1:100} = {targ};

%aTarg = {idx, liftPoint, targs};
%cellTarg{1:maxNumOfHypoAccum+maxNumOfHypo} = {aTarg};

%cellHypoSeed = cell(1,maxNumOfHypoSeed); %{{[]}}; %memory allocation{{[]}};
%cellHypoSeed{1:maxNumOfHypoSeed} = {hypo};
cellHypoSeed = {{[]}};
cellHypoSeed_IMM = {{[]}};

%cellTargSeed = cell(1,maxNumOfHypoSeed); %{{[]}}; %memory allocation{{[]}};
%cellTargSeed{1:maxNumOfHypoSeed} = {aTarg};
cellTargSeed = {{[]}};
cellTargSeed_IMM = {{[]}};

%cellEstm = cell(1,nStep*5);		% 每个周期最多M=4个估计 260*4 == 1100留点余量
cellEstm = cell(1,0);
cellEstm_IMM = cell(1,0);
%cellEstm = {aTarg};
%allHypo = cell(21,100);
%% fill in 1st cell in cellHypo: {{asso prob}  }
%% note: the value of root prob is not important
%cellHypo = {{(1:numTarget)' 0}}; 
stateEstmIMMRun = zeros(nStep-3,TargetNum,2);
stateEstmRun = zeros(nStep-3,TargetNum,2);
stateEstmVBRun = zeros(nStep,TargetNum,2);

%maxL = 20; %最多有20个状态
%% fill in 1st cell in cellTarg: {{numTarg*{idx lifePoint X P}}  }
maxLifePoint = 3;
%cellTmp = cell(1, numTarget); % number of target is known and constant 创建1*numTarget维的cell矩阵
% for i = 1 : numTarget
%     tmp = state{i};             % 第i个目标的状态
%     cellTmp{i} = {i maxLifePoint tmp(:,1) P currentState L cellCtr...
%         cellTrans alpha beta gamma M auxPhi auxTau auxZeta};
% end
%% currentState 为整数
%cellTarg = {cellTmp};       % 新的cell，cellTarg，包含一个元素，这个元素是cell，cellTmp

%% fill in 1st cell in cellEstm: {numTarg*{idx startTime matX}}
%cellEstm = cell(1, numTarget); 
% for i = 1 : numTarget
%     tmp = state{i};
%     cellEstm{i} = {i 0 tmp(:, 1)};  % 第一个时刻的状态
% end

%% evolution stage
% head = 1;
% rear = 1;

% figure(3);
% legend('真实航迹-','估计航迹o');
% if nTarg == 1
%     plot(state{1}(1, 1:nStep-1), state{1}(4, 1:nStep-1), 'k-');
%     hold on
%     drawnow
% end
% if nTarg == 2
%     plot(state{1}(1, 1:nStep-1), state{1}(4, 1:nStep-1), 'k-');
%     hold on
%     drawnow
%     plot(state{2}(1, 1:nStep-1), state{2}(4, 1:nStep-1), 'k-');
%     hold on
%     drawnow
%     xlabel('x direction');
%     ylabel('y direction');
% end
% axis([48000,78000,-10000,20000]);
% RMS1 = zeros(nStep, nTarg);
% RMS2 = zeros(nStep, nTarg);
% IMM methods
% RMS3 = zeros(nStep, nTarg);
% RMS4 = zeros(nStep, nTarg);
for t = 1 : nStep
    % hypoes and estims formed last step as seeds to generate new ones
    tStart = tic;
    disp(t);
if P_include    
    cellHypoSeed = {cellHypo{head:rear}};   % 时刻4时，为6~21的16个seed
end
if IMM_include    
    cellHypoSeed_IMM = {cellHypo_IMM{head_IMM:rear_IMM}};   % 时刻4时，为6~21的16个seed
end
       %% HDP-HMM的方法,采用100个粒子,在cellTarg里面的每个Targ由100个粒子组成
if P_include    
    cellTargSeed = {cellTarg{head:rear}};
end
if IMM_include    
    cellTargSeed_IMM = {cellTarg_IMM{head_IMM:rear_IMM}};
end
if IMM_include    
    if ~isempty(cellTargSeed_IMM{(rear_IMM - head_IMM + 1)}{1}) %% 是否有需要预测
    %% predict stage of Kalman filter
       cellTargSeed_IMM = KF_MHT_Predict_IMM(cellTargSeed_IMM, F_CV,F_CA, Q_CV,Q_CA,T);  % 卡尔曼滤波预测
    %% 针对HDP粒子方法,不需要预测步,可以直接计算predicit likelihood,PL
    end
end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% generate new hypoes
      Z = meas{t}; % get measurement of this time step 新的量测
    if isempty(Z)
        error('There in no measurement at step %d.', t);
    else
if P_include        
        cellHypoNew = GenHypo(cellHypoSeed, cellTargSeed, M,...
            Z, H, R, Pd, densNew, densClt);  %% 时刻4时，为16个seed，64个new hypo
end
if IMM_include
        cellHypoNew_IMM = GenHypo_IMM(cellHypoSeed_IMM, cellTargSeed_IMM, M,...
            Z, H, R, Pd, densNew, densClt);  %% 时刻4时，为16个seed，64个new hypo
end            
	        %% 要将Target_update与prue掉个个，被删除的假设对应的target不更新
if P_include            
	        cellHypo = [cellHypo, cellHypoNew]; %% 时刻4，前面的1+4+16个假设加上新生的64个假设
end
if IMM_include
	        cellHypo_IMM = [cellHypo_IMM, cellHypoNew_IMM]; %% 时刻4，前面的1+4+16个假设加上新生的64个假设
end
    end
if P_include    
    if t < N
        if t == 1
            chooseIdx = [];
            chooseIdx_IMM = [];
            idx = 0;
            idx_IMM = 0;
        end
        chooseIdx = [1:idx+M^t];
    else
        chooseIdx = GetMinProbIdx(cellHypo,M,N,t);
    end
end
    %% update stage of Kalman filter
    %% 新增加一个参数chooseIdx,在chooseIdx范围中的targ才更新
if P_include    
    cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,...
        M, Z, H, R, maxLifePoint,chooseIdx,t); %% expand cellTarg and cellHypo. NOTE THE DIFFERENCE BETWEEN {} & []
end
if IMM_include
    cellTargNew_IMM = KF_MHT_Update_IMM(cellTargSeed_IMM, cellHypoNew_IMM,...
        M, Z, H, R, maxLifePoint,T); %% expand cellTarg and cellHypo. NOTE THE DIFFERENCE BETWEEN {} & []
end
if P_include    
    cellTarg = [cellTarg, cellTargNew];
end
if IMM_include
    cellTarg_IMM = [cellTarg_IMM, cellTargNew_IMM];
end
% 
    %% prune and update
    if t < N
if P_include        
        head = rear + 1;
        rear = rear + M^t;
end
if IMM_include
        head_IMM = rear_IMM + 1;
        rear_IMM = rear_IMM + M^t;
end
    else %%  positions of head and rear are now fixed
if P_include        
        [cellEstm, cellHypo, cellTarg] = Prune(cellEstm,...
            cellHypo, cellTarg, M, N, t, maxLifePoint);
end
if IMM_include
        [cellEstm_IMM, cellHypo_IMM, cellTarg_IMM] = Prune_IMM(cellEstm_IMM,...
            cellHypo_IMM, cellTarg_IMM, M, N, t, maxLifePoint);
end
    end
    % VB Algorithm
if VB_include    
    if(t == 1)  %第一次先做初始化
       % 初始化所有用到的参数
       [newTarg, newState, ZetaCD, VarepsilonAB, Xi, Varphi] = initParaVB(state);
    else
       [newTarg,newState,ZetaCD,VarepsilonAB,Xi, Varphi] ...
           = KF_MHT_Update_VB(newTarg,newState,...
            Z,ZetaCD, VarepsilonAB,Xi, Varphi, t);   
    end
end    
    tElapsed = toc(tStart);
    disp(tElapsed);
    %画图HDP ---- 在图三画观测，没必要，如果跟踪连续，在调试状态下可以
%     for i = 1 : size(Z,2)
%         plot(Z(1, i), Z(2, i), 'ro');
%         hold on
%         drawnow
%     end
    %画图IMM
% 	for i = 1 : size(cellEstm_IMM, 2)
%         stateEstmRaw = cellEstm_IMM{i}{3}{1};
% 	    modeP = cellEstm_IMM{i}{3}{4};
%         XX = zeros(3*dim,1);
%         for k=1:IMMNUM
%              XX = XX + stateEstmRaw(3*dim*(k-1)+1:3*dim*k)*modeP(k);
%         end
% %     RMS = zeros(dim, nStep, nTarg);
%             if i == 2
% %                RMS(t,i) = norm(stateEstm-state{i}(1,t));
%                 RMS3(t,i) = ((stateEstm(1)-state{i}(1,t-2))^2+...
%                     (stateEstm(4)-state{i}(4,t-2))^2);
%                 RMS4(t,i) = ((stateEstm(1)-state{1}(1,t-2))^2+...
%                     (stateEstm(4)-state{1}(4,t-2))^2);
%                 if t >=4
%                     RMSE3(i) = (stateEstm(1)-state{i}(1,t-2))^2+...
%                     (stateEstm(4)-state{i}(4,t-2))^2 + RMSE3(i);
%                     RMSE4(i) = (stateEstm(1)-state{1}(1,t-2))^2+...
%                     (stateEstm(4)-state{1}(4,t-2))^2 + RMSE4(i);
%                 end
% if iCycle == 1
%                plot(stateEstm(1, :), stateEstm(4, :), 'ro');
%                 hold on
%                 drawnow
% end
%             end
%             if i == 1
%                 RMS3(t,i) = ((stateEstm(1)-state{i}(1,t-2))^2+...
%                     (stateEstm(4)-state{i}(4,t-2))^2);
%                 RMS4(t,i) = ((stateEstm(1)-state{2}(1,t-2))^2+...
%                     (stateEstm(4)-state{2}(4,t-2))^2);
%                 if t >=4
%                     RMSE3(i) = (stateEstm(1)-state{i}(1,t-2))^2+...
%                     (stateEstm(4)-state{i}(4,t-2))^2 + RMSE3(i);
%                     RMSE4(i) = (stateEstm(1)-state{2}(1,t-2))^2+...
%                     (stateEstm(4)-state{2}(4,t-2))^2 + RMSE4(i);
%                 end
% if iCycle == 1
%                 plot(stateEstm(1, :), stateEstm(4, :), 'g*');
%                 hold on
%                 drawnow
%  end
%             end
%             if i >= 3
% if iCycle == 1
%                 plot(stateEstm(1, :), stateEstm(4, :), 'bx');
%                 hold on
%                 drawnow
% end
%             end
%     end
    %画图IMM
if IMM_include    
	for i = 1 : size(cellEstm_IMM, 2)
		cellIMM = cellEstm_IMM{i}{3}; %% 第三个参数是100个粒子的数组
	    X = cellIMM{1};
	    modeP = cellIMM{4};
        stateEstmIMM = zeros(6,1);
  		XX = zeros(3*dim,1);
			for k=1:IMMNUM
			  XX = XX + X(3*dim*(k-1)+1:3*dim*k,1)*modeP(k);
			end
			stateEstmIMM = XX;
%     RMS = zeros(dim, nStep, nTarg);
            if i == 2
%                RMS(t,i) = norm(stateEstm-state{i}(1,t));
                RMS3(t,i) = ((stateEstmIMM(1)-state{i}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{i}(4,t-2))^2);
                RMS4(t,i) = ((stateEstmIMM(1)-state{1}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{1}(4,t-2))^2);
                if t >=4
                    RMSE3(i) = (stateEstmIMM(1)-state{i}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{i}(4,t-2))^2 + RMSE3(i);
                    RMSE4(i) = (stateEstmIMM(1)-state{1}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{1}(4,t-2))^2 + RMSE4(i);
                end
                plotState3(:,t-2) = [stateEstmIMM(1);stateEstmIMM(4)];
                if iCycle == 1
                    plot(stateEstmIMM(1, :), stateEstmIMM(4, :), 'ro');
                    % 定义变量stateEstmIMM是为了在RunIt里面一次性完成所有数据的绘图
                    stateEstmIMMRun(t-2,i,:) = [stateEstmIMM(1, :), stateEstmIMM(4, :)];
                    hold on
                    drawnow
                end
            end
            if i == 1
                RMS3(t,i) = ((stateEstmIMM(1)-state{i}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{i}(4,t-2))^2);
                RMS4(t,i) = ((stateEstmIMM(1)-state{2}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{2}(4,t-2))^2);
                if t >=4
                    RMSE3(i) = (stateEstmIMM(1)-state{i}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{i}(4,t-2))^2 + RMSE3(i);
                    RMSE4(i) = (stateEstmIMM(1)-state{2}(1,t-2))^2+...
                    (stateEstmIMM(4)-state{2}(4,t-2))^2 + RMSE4(i);
                end
                plotState4(:,t-2) = [stateEstmIMM(1);stateEstmIMM(4)];
                if iCycle == 1
                    plot(stateEstmIMM(1, :), stateEstmIMM(4, :), 'bo');
                    % 定义变量stateEstmIMM是为了在RunIt里面一次性完成所有数据的绘图
                    stateEstmIMMRun(t-2,i,:) = [stateEstmIMM(1, :), stateEstmIMM(4, :)];
                    hold on
                    drawnow
                end
            end
            if i >= 3
%                 plot(stateEstmIMM(1, :), stateEstmIMM(4, :), 'b+');
%                 hold on
%                 drawnow
            end
    end
end
if P_include
    % HDP-HMM粒子滤波算法
	for i = 1 : size(cellEstm, 2)
		cellSample = cellEstm{i}{3}; %% 第三个参数是100个粒子的数组
        stateEstm = zeros(6,1);
        if size(cellSample,2) == J
            for j = 1 : J
                stateEstm = stateEstm + cellSample{j}{1}/J;
            end
        else
                stateEstm = cellSample{1}{1};
        end
%     RMS = zeros(dim, nStep, nTarg);
        % 注释，由于MHT算法的滞后效果，第一个目标估计在第三个周期，及第三个
        % 周期得到第一个周期的一个目标估计，第四个周期得到第二个周期的两个
        % 目标估计，这样就有了i==1实际上在判断是否是第三个周期
        if i == 2
%                RMS(t,i) = norm(stateEstm-state{i}(1,t));
            RMS1(t,i) = ((stateEstm(1)-state{i}(1,t-2))^2+...
                (stateEstm(4)-state{i}(4,t-2))^2);
            RMS2(t,i) = ((stateEstm(1)-state{1}(1,t-2))^2+...
                (stateEstm(4)-state{1}(4,t-2))^2);
            if t >=4
                RMSE1(i) = (stateEstm(1)-state{i}(1,t-2))^2+...
                (stateEstm(4)-state{i}(4,t-2))^2 + RMSE1(i);
                RMSE2(i) = (stateEstm(1)-state{1}(1,t-2))^2+...
                (stateEstm(4)-state{1}(4,t-2))^2 + RMSE2(i);
            end
            plotState2(:,t-2) = [stateEstm(1);stateEstm(4)];
            if iCycle == 1
                plot(stateEstm(1, :), stateEstm(4, :), 'rp');
                stateEstmRun(t-2,i,:) = [stateEstm(1, :), stateEstm(4, :)]; 
                hold on
                drawnow
            end
        end
        if i == 1
            RMS1(t,i) = ((stateEstm(1)-state{i}(1,t-2))^2+...
                (stateEstm(4)-state{i}(4,t-2))^2);
            RMS2(t,i) = ((stateEstm(1)-state{2}(1,t-2))^2+...
                (stateEstm(4)-state{2}(4,t-2))^2);
            if t >=4
                RMSE1(i) = (stateEstm(1)-state{i}(1,t-2))^2+...
                (stateEstm(4)-state{i}(4,t-2))^2 + RMSE1(i);
                RMSE2(i) = (stateEstm(1)-state{2}(1,t-2))^2+...
                (stateEstm(4)-state{2}(4,t-2))^2 + RMSE2(i);
            end
            plotState1(:,t-2) = [stateEstm(1);stateEstm(4)];
            if iCycle == 1
                plot(stateEstm(1, :), stateEstm(4, :), 'bp');
                stateEstmRun(t-2,i,:) = [stateEstm(1, :), stateEstm(4, :)]; 
                hold on
                drawnow
             end
        end
            if i >= 3
                if iCycle == 1
%                 plot(stateEstm(1, :), stateEstm(4, :), 'bx');
%                 hold on
%                 drawnow
                end
            end
    end
end    
if VB_include
    %画图HDP-PMHT VB
	for i = 1 : size(newTarg, 2)
        stateEstmVB = zeros(6,1);
        stateEstmVB = newTarg{i}{1};
%     RMS = zeros(dim, nStep, nTarg);
        if i == 2
%                RMS(t,i) = norm(stateEstm-state{i}(1,t));
            RMS5(t,i) = ((stateEstmVB(1)-state{i}(1,t))^2+...
                (stateEstmVB(4)-state{i}(4,t))^2);
            RMS6(t,i) = ((stateEstmVB(1)-state{1}(1,t))^2+...
                (stateEstmVB(4)-state{1}(4,t))^2);
            if t >=4
                RMSE5(i) = (stateEstmVB(1)-state{i}(1,t))^2+...
                (stateEstmVB(4)-state{i}(4,t))^2 + RMSE5(i);
                RMSE6(i) = (stateEstmVB(1)-state{1}(1,t))^2+...
                (stateEstmVB(4)-state{1}(4,t))^2 + RMSE6(i);
            end
            plotState5(:,t) = [stateEstmVB(1);stateEstmVB(4)];
            if iCycle == 1
                plot(stateEstmVB(1, :), stateEstmVB(4, :), 'r<');
                plot(state{i}(1,t),state{i}(4,t),'r*');
                stateEstmVBRun(t,i,:) = [stateEstmVB(1, :), stateEstmVB(4, :)];
                hold on
                drawnow
            end
        end
        if i == 1
            RMS5(t,i) = ((stateEstmVB(1)-state{i}(1,t))^2+...
                (stateEstmVB(4)-state{i}(4,t))^2);
            RMS6(t,i) = ((stateEstmVB(1)-state{2}(1,t))^2+...
                (stateEstmVB(4)-state{2}(4,t))^2);
            if t >=4
                RMSE5(i) = (stateEstmVB(1)-state{i}(1,t))^2+...
                (stateEstmVB(4)-state{i}(4,t))^2 + RMSE5(i);
                RMSE6(i) = (stateEstmVB(1)-state{2}(1,t))^2+...
                (stateEstmVB(4)-state{2}(4,t))^2 + RMSE6(i);
            end
            plotState6(:,t) = [stateEstmVB(1);stateEstmVB(4)];
            if iCycle == 1
                plot(stateEstmVB(1, :), stateEstmVB(4, :), 'b>');
                plot(state{i}(1,t),state{i}(4,t),'b*');
                stateEstmVBRun(t,i,:) = [stateEstmVB(1, :), stateEstmVB(4, :)];
                hold on
                drawnow
            end
        end
            if i >= 3
                if iCycle == 1
%                 plot(stateEstm(1, :), stateEstm(4, :), 'bx');
%                 hold on
%                 drawnow
                end
            end
    end
end
if 0
    %画柱状图
    if t == nStep
        hold off
        totalNum = 15;
        cellStatesNum = zeros(totalNum,nTarg);
        for i = 1 : nTarg
            cellSample = cellEstm{i}{3}; %% 第三个参数是100个粒子的数组
            for j = 1 : J
                cStates = cellSample{j}{18};
                cStates = cStates > 4;
                tempCell = cumsum(cStates);
                cellStatesNum(tempCell(end),i) = ...
                    cellStatesNum(tempCell(end),i) + 1;
            end            
        end
        figure(4);
        y = cellStatesNum;
%       y=[300 311;390 425; 312 321; 250 185; 550 535; 420 432; 410 520;];
        b=bar(y);
        ch = get(b,'children');
        set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'})
%        set(ch,'FaceVertexCData',[1 0 1;0 0 0;])
        legend('第一个目标','第二个目标');
        xlabel('Number of clusters ');
        ylabel('Counts over 250 cycles');
        hold off
        figure(7);
        y = cellStatesNum;
%       y=[300 311;390 425; 312 321; 250 185; 550 535; 420 432; 410 520;];
        b=bar(y);
        ch = get(b,'children');
        set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'})
%        set(ch,'FaceVertexCData',[1 0 1;0 0 0;])
        legend('第一个目标','第二个目标');
        xlabel('Number of clusters ');
        ylabel('Counts over 250 cycles');
        hold off
    end
end
    %统计误差
%     RMS = zeros(dim, nStep, nTarg);
if 0
    if t == nStep
        figure(5);
        hold on
        plot(4:nStep,RMS1(4:nStep,1),'g--',4:nStep,RMS1(4:nStep,2),'r-');
%        plot(4:nStep,RMS3(4:nStep,1),'k--',4:nStep,RMS3(4:nStep,2),'b-');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,0,2000]),title('RMS error of position (m)')
        xlabel('time');
        ylabel('Position RMSE');
        hold off
    end
    if t == nStep
        figure(6);
        hold on
        plot(4:nStep,RMS2(4:nStep,1),'g--',4:nStep,RMS2(4:nStep,2),'r-');
%        plot(4:nStep,RMS4(4:nStep,1),'k--',4:nStep,RMS4(4:nStep,2),'b-');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,0,2000]),title('RMS error of position (m)')
        xlabel('time');
        ylabel('Position RMSE');
        hold off
    end
    if t == nStep
        figure(7);
        hold on
        plot(4:nStep,RMS5(4:nStep,1),'m--',4:nStep,RMS5(4:nStep,2),'y-');
%        plot(4:nStep,RMS4(4:nStep,1),'k--',4:nStep,RMS4(4:nStep,2),'b-');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,0,2000]),title('RMS error of position (m)')
        xlabel('time');
        ylabel('Position RMSE');
        hold off
    end
end
end       




