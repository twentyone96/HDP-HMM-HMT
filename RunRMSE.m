clc;clear;close all;
global iCycle
% 15与11重复了，15变成115-,22变成116-,26变成117-
% 7与5重复了，7变成118-,8变成119-,12变成120-,16变成121-,19变成122-,23变成123-,29变成124
global nStep
global nTarg
nStepEnd = 251;11;%21%251;
nStep = 250;20;%30%260;
dim = 2;
nTarg = 2;
iCycle = 124; %RMS3 target2的第十五个有问题%%标志
T=1;
%HDP-HMM MCMC methods
RMS1 = zeros(nStep, nTarg,iCycle);
RMS1Arv = zeros(3,nStep,nTarg); % 1 代表前半段，2代表后半段，3代表总的
RMS1T1 = 0;
RMS1T2 = 0;
RMS1T1F = 0; %first segementation
RMS1T2F = 0;
RMS1T1S = 0; %second segementation
RMS1T2S = 0;

% IMM methods
RMS3 = zeros(nStep, nTarg,iCycle);
RMS3Arv = zeros(3,nStep,nTarg);
RMS3T1 = 0;
RMS3T2 = 0;
RMS3T1F = 0;
RMS3T2F = 0;
RMS3T1S = 0;
RMS3T2S = 0;
% HDP-HMM-VB methods
RMS5 = zeros(nStep, nTarg,iCycle);
RMS5Arv = zeros(3,nStep,nTarg);
RMS5T1 = 0;
RMS5T2 = 0;
RMS5T1F = 0;
RMS5T2F = 0;
RMS5T1S = 0;
RMS5T2S = 0;

RMS1 = RMS1F();
RMS3 = RMS3F();
RMS5 = RMS5F();
%%%RMS1Start

%%%RMS3Start


%%%RMS5Start

% i=iCycle;
    % 这是最后一次循环，计算多轮循环的均方值
       % 总的
       for j=1:nStep
            for k=1:iCycle
                RMS1Arv(3,j,1) = RMS1(j,1,k)/(iCycle-0) +  RMS1Arv(3,j,1);
                RMS1Arv(3,j,2) = RMS1(j,2,k)/(iCycle-0) +  RMS1Arv(3,j,2);

                RMS3Arv(3,j,1) = RMS3(j,1,k)/(iCycle-0) +  RMS3Arv(3,j,1);
                RMS3Arv(3,j,2) = RMS3(j,2,k)/(iCycle-0) +  RMS3Arv(3,j,2);

                RMS5Arv(3,j,1) = RMS5(j,1,k)/(iCycle-0) +  RMS5Arv(3,j,1);
                RMS5Arv(3,j,2) = RMS5(j,2,k)/(iCycle-0) +  RMS5Arv(3,j,2);
            end
            RMS1Arv(3,j,:) = sqrt(RMS1Arv(3,j,:));
            RMS3Arv(3,j,:) = sqrt(RMS3Arv(3,j,:));
            RMS5Arv(3,j,:) = sqrt(RMS5Arv(3,j,:));
       end
       %前半段
       for j=1:110
            for k=1:iCycle
                RMS1Arv(1,j,1) = RMS1(j,1,k)/(iCycle-0) +  RMS1Arv(1,j,1);
                RMS1Arv(1,j,2) = RMS1(j,2,k)/(iCycle-0) +  RMS1Arv(1,j,2);

                RMS3Arv(1,j,1) = RMS3(j,1,k)/(iCycle-0) +  RMS3Arv(1,j,1);
                RMS3Arv(1,j,2) = RMS3(j,2,k)/(iCycle-0) +  RMS3Arv(1,j,2);

                RMS5Arv(1,j,1) = RMS5(j,1,k)/(iCycle-0) +  RMS5Arv(1,j,1);
                RMS5Arv(1,j,2) = RMS5(j,2,k)/(iCycle-0) +  RMS5Arv(1,j,2);
            end
            RMS1Arv(1,j,:) = sqrt(RMS1Arv(1,j,:));
            RMS3Arv(1,j,:) = sqrt(RMS3Arv(1,j,:));
            RMS5Arv(1,j,:) = sqrt(RMS5Arv(1,j,:));
        end
       for j=111:nStep
            for k=1:iCycle
                RMS1Arv(2,j,1) = RMS1(j,1,k)/(iCycle-0) +  RMS1Arv(2,j,1);
                RMS1Arv(2,j,2) = RMS1(j,2,k)/(iCycle-0) +  RMS1Arv(2,j,2);

                RMS3Arv(2,j,1) = RMS3(j,1,k)/(iCycle-0) +  RMS3Arv(2,j,1);
                RMS3Arv(2,j,2) = RMS3(j,2,k)/(iCycle-0) +  RMS3Arv(2,j,2);

                RMS5Arv(2,j,1) = RMS5(j,1,k)/(iCycle-0) +  RMS5Arv(2,j,1);
                RMS5Arv(2,j,2) = RMS5(j,2,k)/(iCycle-0) +  RMS5Arv(2,j,2);
            end
            RMS1Arv(2,j,:) = sqrt(RMS1Arv(2,j,:));
            RMS3Arv(2,j,:) = sqrt(RMS3Arv(2,j,:));
            RMS5Arv(2,j,:) = sqrt(RMS5Arv(2,j,:));
        end
       
       for j = 4 : nStep - 2
            RMS1T1 = RMS1Arv(3,j,1)/(nStep-5) + RMS1T1;
            RMS3T1 = RMS3Arv(3,j,1)/(nStep-5) + RMS3T1;
            RMS5T1 = RMS5Arv(3,j,1)/(nStep-5) + RMS5T1;
            RMS1T2 = RMS1Arv(3,j,2)/(nStep-5) + RMS1T2;
            RMS3T2 = RMS3Arv(3,j,2)/(nStep-5) + RMS3T2;
            RMS5T2 = RMS5Arv(3,j,2)/(nStep-5) + RMS5T2;
        end
       for j = 4 : 110
            RMS1T1F = RMS1Arv(1,j,1)/(110-3) + RMS1T1F;
            RMS3T1F = RMS3Arv(1,j,1)/(110-3) + RMS3T1F;
            RMS5T1F = RMS5Arv(1,j,1)/(110-3) + RMS5T1F;
            RMS1T2F = RMS1Arv(1,j,2)/(110-3) + RMS1T2F;
            RMS3T2F = RMS3Arv(1,j,2)/(110-3) + RMS3T2F;
            RMS5T2F = RMS5Arv(1,j,2)/(110-3) + RMS5T2F;
        end
       for j = 111 : nStep - 2
            RMS1T1S = RMS1Arv(2,j,1)/(140-2) + RMS1T1S;
            RMS3T1S = RMS3Arv(2,j,1)/(140-2) + RMS3T1S;
            RMS5T1S = RMS5Arv(2,j,1)/(140-2) + RMS5T1S;
            RMS1T2S = RMS1Arv(2,j,2)/(140-2) + RMS1T2S;
            RMS3T2S = RMS3Arv(2,j,2)/(140-2) + RMS3T2S;
            RMS5T2S = RMS5Arv(2,j,2)/(140-2) + RMS5T2S;
        end
        figure(5);
        hold on
        p1 = plot(4:nStep,RMS1Arv(3,4:nStep,1),'g-');
        p2 = plot(4:nStep,RMS3Arv(3,4:nStep,1),'k-');
        p3 = plot(4:nStep,RMS5Arv(3,4:nStep,1),'m-');
        p4 = plot(4:110,RMS1Arv(1,4:110,1),'g--');
        p5 = plot(4:110,RMS3Arv(1,4:110,1),'k--');
        p6 = plot(4:110,RMS5Arv(1,4:110,1),'m--');
        p7 = plot(111:248,RMS1Arv(2,111:248,1),'g--');
        p8 = plot(111:248,RMS3Arv(2,111:248,1),'k--');
        p9 = plot(111:248,RMS5Arv(2,111:248,1),'m--');
        % 画出机动状态的分割线
        line([10,10], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([20,20], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([30,30], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([50,50], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([60,60], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([70,70], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([80,80], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([110,110], [0,300], 'color', [.5,.5,.5], 'linestyle', '--');
        line([120,120], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([130,130], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([140,140], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([160,160], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([170,170], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([180,180], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([190,190], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([220,220], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([230,230], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([240,240], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,50,350]),title('Target 1')
        xlabel('time');
        ylabel('Position RMSE (m)');
        ax1 = gca;
        ax2 = axes( 'Position',get(ax1,'Position'),'Visible','off');
        ax3 = axes( 'Position',get(ax1,'Position'),'Visible','off');
        str1 = strcat('IMM-MHT target 1',32,32,32,32,32,num2str(roundn(RMS1T1,-2)),'m');
        str2 = strcat('HDP-HMM-MHT target 1',32,32,num2str(roundn(RMS3T1,-2)),'m');
        str3 = strcat('HDP-HMM-PMHT target 1',32,num2str(roundn(RMS5T1,-2)),'m');
        legend(ax1,[p1,p2,p3],str1, str2, str3);         
        drawnow
        hold on
        str4 = strcat('IMM-MHT, first half',32,32,32,32,32,num2str(roundn(RMS1T1F,-2)),'m');
        str5 = strcat('HDP-HMM-MHT, first half',32,32,num2str(roundn(RMS3T1F,-2)),'m');
        str6 = strcat('HDP-HMM-PMHT, first half',32,num2str(roundn(RMS5T1F,-2)),'m');
%        legend(ax2,[p4,p5,p6],str4, str5, str6);%,'Position',[200 600 0.1 0.2]);         
        strings={str4;str5;str6};
        annotation('textbox',[0.3,0.7,0.42,0.15],'LineStyle','-',...
   'LineWidth',.5,'String',strings);
        drawnow
        hold on
        str7 = strcat('IMM-MHT, second half',32,32,32,32,32,num2str(roundn(RMS1T1S,-2)),'m');
        str8 = strcat('HDP-HMM-MHT, second half',32,32,num2str(roundn(RMS3T1S,-2)),'m');
        str9 = strcat('HDP-HMM-PMHT, second half',32,num2str(roundn(RMS5T1S,-2)),'m');
        strings={str7;str8;str9};
        annotation('textbox',[0.3,0.7,0.44,0.15],'LineStyle','-',...
   'LineWidth',.5,'String',strings);
 %       legend(ax3,[p7,p8,p9],str7, str8, str9,'Position');%,[6 16 0.1 0.2]);         

%         legend('IMM-MHT target 1',...
%             'HDP-HMM-MHT target 1',...
%             'HDP-HMM-PMHT target 1');
        drawnow
        hold on
        hold off
        figure(55);
        hold on
        p1 = plot(4:nStep,RMS1Arv(3,4:nStep,2),'g-');
        p2 = plot(4:nStep,RMS3Arv(3,4:nStep,2),'k-');
        p3 = plot(4:nStep,RMS5Arv(3,4:nStep,2),'m-');
        % 画出机动状态的分割线
        line([10,10], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([20,20], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([30,30], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([50,50], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([60,60], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([70,70], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([80,80], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([110,110], [0,300], 'color', [.5,.5,.5], 'linestyle', '--');
        line([120,120], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([130,130], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([140,140], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([160,160], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([170,170], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([180,180], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([190,190], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([220,220], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([230,230], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([240,240], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        axis([0,nStep,50,350]),title('Target 2')
        xlabel('time');
        ylabel('Position RMSE (m)');
        hold off
        ax1 = gca;
        ax2 = axes( 'Position',get(ax1,'Position'),'Visible','off');
        ax3 = axes( 'Position',get(ax1,'Position'),'Visible','off');
        str1 = strcat('IMM-MHT target 2',32,32,32,32,32,num2str(roundn(RMS1T2,-2)),'m');
        str2 = strcat('HDP-HMM-MHT target 2',32,32,num2str(roundn(RMS3T2,-2)),'m');
        str3 = strcat('HDP-HMM-PMHT target 2',32,num2str(roundn(RMS5T2,-2)),'m');
        legend(ax1,[p1,p2,p3],str1, str2, str3);     
        drawnow
        hold on

        str4 = strcat('IMM-MHT, first half',32,32,32,32,32,num2str(roundn(RMS1T2F,-2)),'m');
        str5 = strcat('HDP-HMM-MHT, first half',32,32,num2str(roundn(RMS3T2F,-2)),'m');
        str6 = strcat('HDP-HMM-PMHT, first half',32,num2str(roundn(RMS5T2F,-2)),'m');
%        legend(ax2,[p1,p2,p3],str4, str5, str6);%,'Position',[0.2 0.6 0.1 0.2]);         
        strings={str4;str5;str6};
        annotation('textbox',[0.3,0.7,0.42,0.15],'LineStyle','-',...
   'LineWidth',.5,'String',strings);
        drawnow
        hold on
        
        str7 = strcat('IMM-MHT, second half',32,32,32,32,32,num2str(roundn(RMS1T2S,-2)),'m');
        str8 = strcat('HDP-HMM-MHT, second half',32,32,num2str(roundn(RMS3T2S,-2)),'m');
        str9 = strcat('HDP-HMM-PMHT, second half',32,num2str(roundn(RMS5T2S,-2)),'m');
%        legend(ax3,[p1,p2,p3],str7, str8, str9);%,'Position',[0.2 0.6 0.1 0.2]);         
        strings={str7;str8;str9};
        annotation('textbox',[0.3,0.7,0.44,0.15],'LineStyle','-',...
   'LineWidth',.5,'String',strings);
        drawnow
        hold on
        hold off
        figure(6);
        hold on
%        plot(4:nStep,RMS2Arv(4:nStep,1),'g--',4:nStep,RMS2Arv(4:nStep,2),'r-');
% % %         plot(4:nStep,RMS1Arv(3,4:nStep,1),'g--',4:nStep,RMS1Arv(3,4:nStep,2),'r-');
% % %         plot(4:nStep,RMS3Arv(3,4:nStep,1),'k--',4:nStep,RMS3Arv(3,4:nStep,2),'b-');
% % %         plot(4:nStep,RMS5Arv(3,4:nStep,1),'m--',4:nStep,RMS5Arv(3,4:nStep,2),'y-');
        plot(4:nStep,(RMS1Arv(3,4:nStep,1)+RMS1Arv(3,4:nStep,2))/2,'g-');
        plot(4:nStep,(RMS3Arv(3,4:nStep,1)+RMS3Arv(3,4:nStep,2))/2,'k-');
        plot(4:nStep,(RMS5Arv(3,4:nStep,1)+RMS5Arv(3,4:nStep,2))/2,'m-');
%         for k=4:nStep
%           plot(k,RMS(k,1),'g*',k,RMS(k,2),'m*');
%         end
        % 画出机动状态的分割线
        line([10,10], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([20,20], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([30,30], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([50,50], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([60,60], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([70,70], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([80,80], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([110,110], [0,300], 'color', [.5,.5,.5], 'linestyle', '--');
        line([120,120], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([130,130], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([140,140], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([160,160], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([170,170], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([180,180], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([190,190], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([220,220], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([230,230], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        line([240,240], [0,200], 'color', [.5,.5,.5], 'linestyle', '--');
        axis([0,nStep,80,250]),title('RMS error of position (m)')
        xlabel('time');
        ylabel('Position RMSE');
        hold off
        legend('IMM-MHT', 'HDP-HMM-MHT', 'HDP-HMM-PMHT');
        drawnow
        hold on
        hold off

fprintf('MHT DONE\n\n');

%参考
