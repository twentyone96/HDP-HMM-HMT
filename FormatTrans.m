% FormatTrans transforms the 3D profile to 2D states and meas.
% Suppose all targets have same lifetime, or the number of steps 
% 'nStep' is constant for all targets

function [state, meas, nStep] = ...
	FormatTrans(X, nTarg, nStep, densClt, Pd)

global maxMes

state = cell(1, nTarg);
xMax = -Inf; xMin = Inf;    % ??????
yMax = -Inf; yMin = Inf;
zMax = -Inf; zMin = Inf;

for i = 1 : nTarg
	%targ = profile{i};
	state{i} = X{i};
    j=find(state{i}(1,:)>0);
	xMax = max(xMax, max(state{i}(1, j)));  
	xMin = min(xMin, min(state{i}(1, j)));
	yMax = max(yMax, max(state{i}(4, j)));
	yMin = min(yMin, min(state{i}(4, j)));
    zMax = max(zMax, max(state{i}(7, j)));
	zMin = min(zMin, min(state{i}(7, j)));
end

meas = cell(1, nStep);         
poissClt = densClt*(xMax-xMin)*(yMax-yMin); 
xIntervel = xMax-xMin;
yIntervel = yMax-yMin;
%x轴和y轴都放大1.1倍
xIntervel = 0.1*xIntervel/2;
yIntervel = 0.1*yIntervel/2;
for i = 1 : nStep               
%    thisMeas = [];
    iStep = i;
    thisMeas = GenMeasurementsPol(state, nTarg, iStep);
%    nClt = poissrnd(poissClt);           
   nClt = 5;
    % make sure nClt > 0 to avoid the case that meas{i} is empty
    while nClt == 0
        nClt = poissrnd(poissClt);            
    end
    if (nClt + size(thisMeas,1)) >= maxMes
    	nClt = maxMes - size(thisMeas,1);
    end
%     nClt = 0;
    cltMeas = [unifrnd(xMin - xIntervel, xMax + xIntervel, 1, nClt); ...
        unifrnd(yMin - yIntervel, yMax + yIntervel, 1, nClt); ...
        unifrnd(zMin, zMax, 1, nClt)]; 
%simaple
     thisMeas = [thisMeas, cltMeas];
	
	meas{i} = thisMeas;
end

% the philosophy of this step is that N points seperate N-1 
% parts. The first meas is for initialization, not estimation.
nStep = nStep - 1;

% plot if allowed
	figure(2);
	for t = 1 : nTarg
        j = find (state{t}(1,:)>0);
		plot(state{t}(1, j), state{t}(4, j),'k*');
		hold on
	end
	for i = 1 : nStep
		measMat = meas{i};
		if size(measMat, 2) > 0
			plot(measMat(1, :), measMat(2, :), 'r*');
			hold on
		end
	end
    drawnow
    axis([49000,73000,-8000,20000])
    legend('Track 1','Track 2','Measurements')
    xlabel('x direction');
    ylabel('y direction');
    hold off
