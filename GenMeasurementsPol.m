% Bayesian approach extended object and cluster tracking using Random
% Matrices, apex E provides the related procedures
function [z_measure] = GenMeasurementsPol(x_state,nTarg, iStep)   %x_state{:}([1 4 7],i))
%GENMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here
  global Pd
  global RadarLocX
  global RadarLocY
  %noise factor 
  global re_d
  global re_t
  global CartX
  global CartY
  
  x_stateP = zeros(3,nTarg); %1-x;2-y;3-z

  for i=1:nTarg
    x_stateP(:,i) = x_state{i}([1 4 7],iStep);
  end
  % cov of x_state
  x_stateSize = nTarg;  %size(x_state);
  x_stateCov = ones(1,x_stateSize);
  % all pairs in list has been removed
  NoMerg = 0;
    
  while ~NoMerg && (x_stateSize >= 2)
    PuList = CalcPuList(x_stateP,x_stateSize);
    [x_stateP, NoMerg, x_stateCov] = MergeMeasure(x_stateP, PuList, x_stateCov);
    x_stateSize = size(x_stateP,2);
  end
  % using cov to rebuild the measurements
  x_statePSize = size(x_stateP);
  x_stateLeft = x_statePSize(2);
  z_measureP = zeros(2, x_stateLeft);
  for k = 1 : x_stateLeft
    [z_measureP(2,k),z_measureP(1,k)] = cart2pol(x_stateP(1,k) - RadarLocX,x_stateP(2,k) - RadarLocY);
    %% 这里控制观测的误差
    z_measureP(1,k) = z_measureP(1,k) + x_stateCov(k) * normrnd(0,re_d^2);
    z_measureP(2,k) = z_measureP(2,k) + x_stateCov(k) * normrnd(0,re_t^2);
    [x_stateP(1,k),x_stateP(2,k)] = pol2cart(z_measureP(2,k),z_measureP(1,k));
    x_stateP(1,k) = z_measureP(1,k)*cos(z_measureP(2,k))+RadarLocX;
    x_stateP(2,k) = z_measureP(1,k)*sin(z_measureP(2,k))+RadarLocY;

    x_stateP(1,k) = x_stateP(1,k)+normrnd(0,CartX);
    x_stateP(2,k) = x_stateP(2,k)+normrnd(0,CartY);
    
  end
  % using Pd to rebuild the measurements
  k = 1;
  while k <= x_stateLeft
    if(x_stateCov(k) == 1)
      tempPd = rand;
      if 0%(tempPd > Pd)   % && x_stateLeft > 2)   %可以考虑所有目标都没有被检测到的情况
        x_stateLeft = x_stateLeft - 1;
        x_stateP(:,x_stateLeft) = [];
      else
        k = k + 1;
      end
    else
      k = k + 1;
    end
  end
  z_measure = x_stateP;
end

function [PuList] = CalcPuList(x_stateP,nTarg)
% calculate Pu of any two measures
  global RadarLocX
  global RadarLocY
  N_target = nTarg;
  
  RangeRes = 100;
  AngleRes = 0.1*pi/180;
  
  measureNum = nTarg;
  z_measureP = zeros(2, N_target); %极坐标观测值，第一个值是距离，第二个值是角度
  idx = 0;
  
  for k = 1 : measureNum
   z_measureP(1,k) = ((x_stateP(1,k) - RadarLocX)^2 + (x_stateP(2,k) - RadarLocY)^2)^(1/2);
   z_measureP(2,k) = atan((x_stateP(2,k) - RadarLocY)/(x_stateP(1,k)-RadarLocX));
  end
  %para1:Pu value; para2, para3:two measure index
  PuList = zeros(3,measureNum*(measureNum-1)*(1/2));
  for k = 1 : measureNum
    for i = k+1 : measureNum
      idx = idx + 1;
      minus1 = abs(z_measureP(1,k)-z_measureP(1,i));
      Bi1 = minus1/RangeRes;
      Squre1 = -0.5*Bi1^2;
      exp1 = exp(Squre1);
      minus2 = abs(z_measureP(2,k)-z_measureP(2,i));
      Bi2 = minus2/AngleRes;
      Squre2 = -0.5*Bi2^2;
      exp2 = exp(Squre2);
      PuList(1,idx) = exp1 * exp2;
      PuList(2,idx) = k;
      PuList(3,idx) = i;
    end
  end
  %PuList = sort(PuList, 4, 'ascend'); %可以不需要对PuList排序
end

function [x_stateP, NoMerg, x_stateCov] = MergeMeasure(x_stateP, PuList, x_stateCov)
% Merge any two measures according Pu   
  PuListNumV = size(PuList);
  PuListNum = PuListNumV(2);
  FReturn = 0;
  global PolNoMerge
  
  while PuListNum && ~FReturn
    PuTemp = rand;
    PuIdx = randi(PuListNum);
    PuIdxValue = PuList(:,PuIdx);
    FReturn = 0; % return the function
    PuIdxValueTemp = PuIdxValue(1);
    if (PuTemp > PuIdxValueTemp || PolNoMerge)%resolvable
%    if 1
      % remove the pair
      PuList(:,PuIdx) = [];
    else
      %if(PuTemp < Pu) %unresovable
      z_measureIdxA = PuIdxValue(2);
      z_measureIdxB = PuIdxValue(3);
      x_stateP(1,z_measureIdxA(1)) = (x_stateP(1,z_measureIdxA(1)) + x_stateP(1,z_measureIdxB(1)))/2;
      x_stateP(2,z_measureIdxA(1)) = (x_stateP(2,z_measureIdxA(1)) + x_stateP(2,z_measureIdxB(1)))/2;
      x_stateCov(z_measureIdxA(1)) = x_stateCov(z_measureIdxA(1)) + x_stateCov(z_measureIdxB(1));
      % merge the pair
      x_stateP(:,z_measureIdxB)=[]; %remove the z_measureIdxB row
      x_stateCov(z_measureIdxB)=[];
      FReturn = 1;
    end
    PuListNumV = size(PuList);
    PuListNum = PuListNumV(2);
  end
  
  if FReturn
    NoMerg = 0;
  else
    NoMerg = 1;
  end
end