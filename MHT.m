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
% VB�㷨״̬����

dim = 2;
%% initialization
%�����ĸ�������ͨ���������ݽ���,R��GenCov����
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

%evolution state ��ʼ�����衢������Ŀ�꣩��
head = 1; %��ʼʱ����1������
rear = 1;
head_IMM = 1; %��ʼʱ����1������
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

%cellEstm = cell(1,nStep*5);		% ÿ���������M=4������ 260*4 == 1100��������
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

%maxL = 20; %�����20��״̬
%% fill in 1st cell in cellTarg: {{numTarg*{idx lifePoint X P}}  }
maxLifePoint = 3;
%cellTmp = cell(1, numTarget); % number of target is known and constant ����1*numTargetά��cell����
% for i = 1 : numTarget
%     tmp = state{i};             % ��i��Ŀ���״̬
%     cellTmp{i} = {i maxLifePoint tmp(:,1) P currentState L cellCtr...
%         cellTrans alpha beta gamma M auxPhi auxTau auxZeta};
% end
%% currentState Ϊ����
%cellTarg = {cellTmp};       % �µ�cell��cellTarg������һ��Ԫ�أ����Ԫ����cell��cellTmp

%% fill in 1st cell in cellEstm: {numTarg*{idx startTime matX}}
%cellEstm = cell(1, numTarget); 
% for i = 1 : numTarget
%     tmp = state{i};
%     cellEstm{i} = {i 0 tmp(:, 1)};  % ��һ��ʱ�̵�״̬
% end

%% evolution stage
% head = 1;
% rear = 1;

% figure(3);
% legend('��ʵ����-','���ƺ���o');
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
    cellHypoSeed = {cellHypo{head:rear}};   % ʱ��4ʱ��Ϊ6~21��16��seed
end
if IMM_include    
    cellHypoSeed_IMM = {cellHypo_IMM{head_IMM:rear_IMM}};   % ʱ��4ʱ��Ϊ6~21��16��seed
end
       %% HDP-HMM�ķ���,����100������,��cellTarg�����ÿ��Targ��100���������
if P_include    
    cellTargSeed = {cellTarg{head:rear}};
end
if IMM_include    
    cellTargSeed_IMM = {cellTarg_IMM{head_IMM:rear_IMM}};
end
if IMM_include    
    if ~isempty(cellTargSeed_IMM{(rear_IMM - head_IMM + 1)}{1}) %% �Ƿ�����ҪԤ��
    %% predict stage of Kalman filter
       cellTargSeed_IMM = KF_MHT_Predict_IMM(cellTargSeed_IMM, F_CV,F_CA, Q_CV,Q_CA,T);  % �������˲�Ԥ��
    %% ���HDP���ӷ���,����ҪԤ�ⲽ,����ֱ�Ӽ���predicit likelihood,PL
    end
end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% generate new hypoes
      Z = meas{t}; % get measurement of this time step �µ�����
    if isempty(Z)
        error('There in no measurement at step %d.', t);
    else
if P_include        
        cellHypoNew = GenHypo(cellHypoSeed, cellTargSeed, M,...
            Z, H, R, Pd, densNew, densClt);  %% ʱ��4ʱ��Ϊ16��seed��64��new hypo
end
if IMM_include
        cellHypoNew_IMM = GenHypo_IMM(cellHypoSeed_IMM, cellTargSeed_IMM, M,...
            Z, H, R, Pd, densNew, densClt);  %% ʱ��4ʱ��Ϊ16��seed��64��new hypo
end            
	        %% Ҫ��Target_update��prue����������ɾ���ļ����Ӧ��target������
if P_include            
	        cellHypo = [cellHypo, cellHypoNew]; %% ʱ��4��ǰ���1+4+16���������������64������
end
if IMM_include
	        cellHypo_IMM = [cellHypo_IMM, cellHypoNew_IMM]; %% ʱ��4��ǰ���1+4+16���������������64������
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
    %% ������һ������chooseIdx,��chooseIdx��Χ�е�targ�Ÿ���
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
    if(t == 1)  %��һ��������ʼ��
       % ��ʼ�������õ��Ĳ���
       [newTarg, newState, ZetaCD, VarepsilonAB, Xi, Varphi] = initParaVB(state);
    else
       [newTarg,newState,ZetaCD,VarepsilonAB,Xi, Varphi] ...
           = KF_MHT_Update_VB(newTarg,newState,...
            Z,ZetaCD, VarepsilonAB,Xi, Varphi, t);   
    end
end    
    tElapsed = toc(tStart);
    disp(tElapsed);
    %��ͼHDP ---- ��ͼ�����۲⣬û��Ҫ����������������ڵ���״̬�¿���
%     for i = 1 : size(Z,2)
%         plot(Z(1, i), Z(2, i), 'ro');
%         hold on
%         drawnow
%     end
    %��ͼIMM
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
    %��ͼIMM
if IMM_include    
	for i = 1 : size(cellEstm_IMM, 2)
		cellIMM = cellEstm_IMM{i}{3}; %% ������������100�����ӵ�����
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
                    % �������stateEstmIMM��Ϊ����RunIt����һ��������������ݵĻ�ͼ
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
                    % �������stateEstmIMM��Ϊ����RunIt����һ��������������ݵĻ�ͼ
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
    % HDP-HMM�����˲��㷨
	for i = 1 : size(cellEstm, 2)
		cellSample = cellEstm{i}{3}; %% ������������100�����ӵ�����
        stateEstm = zeros(6,1);
        if size(cellSample,2) == J
            for j = 1 : J
                stateEstm = stateEstm + cellSample{j}{1}/J;
            end
        else
                stateEstm = cellSample{1}{1};
        end
%     RMS = zeros(dim, nStep, nTarg);
        % ע�ͣ�����MHT�㷨���ͺ�Ч������һ��Ŀ������ڵ��������ڣ���������
        % ���ڵõ���һ�����ڵ�һ��Ŀ����ƣ����ĸ����ڵõ��ڶ������ڵ�����
        % Ŀ����ƣ�����������i==1ʵ�������ж��Ƿ��ǵ���������
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
    %��ͼHDP-PMHT VB
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
    %����״ͼ
    if t == nStep
        hold off
        totalNum = 15;
        cellStatesNum = zeros(totalNum,nTarg);
        for i = 1 : nTarg
            cellSample = cellEstm{i}{3}; %% ������������100�����ӵ�����
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
        legend('��һ��Ŀ��','�ڶ���Ŀ��');
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
        legend('��һ��Ŀ��','�ڶ���Ŀ��');
        xlabel('Number of clusters ');
        ylabel('Counts over 250 cycles');
        hold off
    end
end
    %ͳ�����
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




