%为每个目标生成轨迹
function [Xset] = GenState(nTarg, nStep, T, F_CV, Q_CV, F_CA, Q_CA)

Xset = cell(1,nTarg);

targ = GenTarg(nTarg,nStep);
%q=10;
%GAMA=[T^2/2;           %过程噪声增益
%      T;    
%      T^2/2;           %过程噪声增益
%      T;];
%PNoise=sqrt(q)*randn(1,100);
%每个target可以单独通过一个训练来构建真实航迹
for i=1:nTarg 
    k = 0;
    if i == 1 || i == 2
        X=zeros(9,nStep);
        X(:,targ(i).start)=targ(i).inistate;
        for j= (targ(i).start+1) : (targ(i).end)
            k = floor(j/10) + 1;
            if k == 1 || k == 3 ||...
               k == 6 || k == 8 ||...
               k == 12 || k == 14 ||...
               k == 17 || k == 19 ||...
               k == 23 || k == 25 || k == 26
                noise = mvnrnd([0; 0; 0; 0; 0; 0;0;0;0], Q_CV, 1)';
                X(:,j)=F_CV*X(:,j-1);%+noise;
            elseif  k == 2 || k == 7  || k == 13 || k == 18 || k == 24
                noise = mvnrnd([0; 0; 0; 0; 0; 0;0;0;0], Q_CA, 1)';
                X(:,j)=F_CA*X(:,j-1);%+noise;
                X(3,j) = 25;
                X(6,j) = 25;
            elseif k == 4 || k == 5 || k == 15 || k == 16
                QCT = GenQCT(-9*pi/180,T);
                noise = mvnrnd([0; 0; 0; 0; 0; 0;0;0;0], QCT, 1)';
                X(:,j)=GenFCT(-9*pi/180,T)*X(:,j-1);%+noise;
            elseif k == 9 || k == 10 || k == 11 || k == 20 || k == 21 || k == 22
                QCT = GenQCT(6*pi/180,T);
                noise = mvnrnd([0; 0; 0; 0; 0; 0;0;0;0], QCT, 1)';
                X(:,j)=GenFCT(6*pi/180,T)*X(:,j-1);%+noise;
            end
        end
    end
    Xset{i}=X; %X为x_position_1,x_vet_2,y_pos_3,y_vet_4,z_pos_5,z_vet_6
end
%目标沿着45度,先10s匀速,在10s加速(25m^2,25m^2),在10s匀速:共30s
%在拐弯180度,?s  20s 9度
%在匀速10s,在减速(-25m^2,-25m^2),在10s匀速:共30s
%在拐弯180度,?s 30s 6度
%在匀速10s,在减速(25m^2,25m^2),在10s匀速:共30s
%在拐弯180度,?s 20s 9度
%在匀速10s,在减速(-25m^2,-25m^2),在10s匀速:共30s
%在拐弯180度,?s 30s 6度
%在匀速10s,在减速(25m^2,25m^2),在10s匀速:共30s
%一共30s*5+?s*4=?s

    figure(1);
    for i=1:nStep
        plot(Xset{1}(1,i),Xset{1}(4,i),'.b');
        hold on
        if nTarg == 2
            plot(Xset{2}(1,i),Xset{2}(4,i),'.r');
    %        plot(Xset{3}(1,i),Xset{3}(3,i));       %第三个目标
            hold on
        end
    end
    axis([48000,78000,-10000,20000]);
    legend('target 1','target 2');
    drawnow
    title('true trajectory');
    xlabel('x direction');
    ylabel('y direction');
    hold off

% 量测噪声只要加在x，y上就行了

end

%
function [targ] = GenTarg(nTarg,nStep)
for i=1:nTarg
    if i==1
        targ(1).inistate = [50000;141.4;0;10000;141.4;0;0;0;0];   % 飞行高度50km，z方向速度为0
        targ(1).start = 1;
        targ(1).end = nStep;
        %可以在增加一个cell，代表机动的时刻或者什么的
    end
    
    if i==2
        targ(2).inistate = [51000;141.4;0;9000;141.4;0;0;0;0];
        targ(2).start = 1;
        targ(2).end = nStep;
    end
    
    if i==3
        targ(3).inistate = [250;5.0;-200.0;8.5;50;0];
        targ(3).start = 20;
        targ(3).end = 80;
    end
end

end