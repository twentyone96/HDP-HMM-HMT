% Analyse performs basic tracking result analysis.
% 
% input: 
%   first & last - only compare real target states and estimations 
%       during time step [first : last].
%   estm - 1*nTarg (or nTarg*1) cell array. Each cell contains a 
%       4*? matrix. Each column of the matrix is estimation for 
%       [x vx y vy]'. 
%   state - real target state recorded in a similar way as estm. 
%   nTarg - number of targets
%   
% output: 
%   errRMS - nTarg*4 vector recording root mean square error of
%       [x vx y vy] for each target. 
%   lose - nTarg*1 binary vector indicating whether losing track 
%       event occurs for each target. 

function [errRMS, lose] = Analyse(first, last, estm, state, nTarg)

errRMS = zeros(nTarg, 4);
lose = zeros(nTarg, 1);

for i = 1 : nTarg
    a = estm{i};
    b = state{i}; 

	% check if losing track on the middle way
	if length(a) < last
        lose(i) = 1; 
        errRMS(i, :) = NaN;
        continue;
	end

	% check if the estimation deviates from the real track too
	% much (losing track). The method is to compare mean position
	% error and mean step displacement of the target: if the 
	% former is larger, then we think the tracking fails. 
	errPos = sqrt((a(1, first:last) - b(1, first:last)).^2 + ...
		(a(3, first:last) - b(3, first:last)).^2);
	difPos = diff(b([1 3], first:last), 1, 2);
	stepPos = sqrt(difPos(1, :).^2 + difPos(2, :).^2);
	if mean(errPos) > mean(stepPos)
		lose(i) = 1;
		errRMS(i, :) = NaN;
		continue;
	end
	
	err = a(:, first:last) - b(:, first:last); % err[x vx y vy]'
    
    errrms=zeros(4,97);
    if i==1
        for n=1:last
            errn = a(:, first:n) - b(:, first:n);
            errrms(:,n)=(sqrt(mean(errn.^2, 2)))';
        end
%        k=1:last-1;
%        figure(5);
%        plot(k,errrms(1,k));
%    legend('估计误差');
%    xlabel('Time(s)');ylabel('VelocityX Var');
%    title('目标1x位置估计误差结果');
%        figure(6);
%        plot(k,errrms(2,k));
%    legend('估计误差');
%    xlabel('Time(s)');ylabel('VelocityX Var');
%    title('目标1x速度估计误差结果');
        
    end
    
	errRMS(i, :) = (sqrt(mean(err.^2, 2)))';
end


function [errRMS, lose] = ShowCtrl(thisTarg, nStep)
global MaxStateNum
global dim
global J
%    ctrlArray = zeros(dim, MaxStateNum*J);
    ctrlArray = zeros(dim, J);
%     for i = 1 : J
%         ctrlArray(:,MaxStateNum*(i-1)+1:MaxStateNum*i) = thisTarg{i};
%     end
    for i = 1 : J
        ctrlArray(:,MaxStateNum*(i-1)+nStep) = thisTarg{i}(nStemp);
    end
    maxArray = (ctrlArray(1,:) > 0) .* sqrt(ctrlArray(1,:))...
        + (ctrlArray(2,:) > 0) .* sqrt(ctrlArray(2,:))
    maxCtrl = find(maxArray == max(maxArray));
    disp(maxCtrl);
end






