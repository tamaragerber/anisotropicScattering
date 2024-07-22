%linenum = 20220625
%          20220627
%          20220628
%          20220630
%          20220701
%          20220704
%          20220705
%          20220706

function prepareProfiles(linenum)

  path2file = 'radarData/';                            % adjust path to radar data

  files = dir([path2file,num2str(linenum),'_*.h5']);
  numSegments = length(files);

  % extract coordinates and Indices every 5 km
  %------------------------------------------------------------------------
  data.lat = [];
  data.lon = [];
  data.Lat = {};
  data.Lon = {};
  
  for i= 1:numSegments
    fname = files(i,1).name;
    lat = h5read([path2file,fname],'/processed/combined/Integrator_2/lat');
    data.lat = [data.lat;lat];
    lon = h5read([path2file,fname],'/processed/combined/Integrator_2/lon');
    data.lon = [data.lon;lon];
    data.Lat{i} = lat;
    data.Lon{i} = lon;
  end
  
  wgs84 = wgs84Ellipsoid;
  data.dist = [0;cumsum(distance(data.lat(1:end-1),data.lon(1:end-1),data.lat(2:end),data.lon(2:end),wgs84))]; % distance along profile
  interv = 0:5000:data.dist(end);
  [~, data.Idx] = min(abs(data.dist-interv)); 
  
  
  % Initialize variables
  %------------------------------------------------------------------------
  cumulativeLengths = cumsum(cellfun(@length, data.Lat));
  Indices = cell(1, numSegments);

  % Loop through each segment and calculate indices
  % -----------------------------------------------------------------------
  for i = 1:numSegments

      endIdx = cumulativeLengths(i);
      p = data.Idx(data.Idx <= endIdx);
      if i > 1
          p = p - cumulativeLengths(i-1);
      end
      p(p <= 0) = [];
      Indices{i} = p;
  end
  
  % plot radar line from all segments and analysis points
  % -----------------------------------------------------------------------
  figure()
  plot(data.lon,data.lat,'.');hold on
  plot(data.lon(data.Idx),data.lat(data.Idx),'ko');
  title(['profile: ',num2str(linenum)])
  xlabel('Longitude')
  ylabel('Latitude')
  set(gca,'fontsize',14');
  legend('radar line','analysis points')
  box on
  grid on
  axis square
  
  for i=1:numSegments
    Lat = data.Lat{1,i};
    Lon = data.Lon{1,i};
    plot(Lon(Indices{i}),Lat(Indices{i}),'r*')
  end
 

  %% initialize parameters
  epsilonPerpendicular = 3.15;                                              % perpendicular component of real part of dielectric permittivity
  epsilonParallel = 3.1840;                                                 % parallel component of real part of dielectric permittivity
  lightSpeed = 299792458;                                                   % Speed of light [m a^-1]
  centerFrequency = 330e6;                                                  % Central frequency [Hz]
  epsilonAverage = mean([epsilonPerpendicular epsilonParallel]);            % average permittivity
  deltaEpsilon = 0.034;                                                     % eiPar-eiPer; Dielectric anisotropy
  dzPowerNorm = 50;                                                         % vertical averaging distance
  cropIdx = 14190;                                                          % cropping index
  dZ = 0.336;                                                               % vertical resolution

  %% point Analysis
  % looping through all segment files
  % -----------------------------------------------------------------------
  for pcount = 1:numSegments
    clear fullDrivingDir drivingDirection
    idx = Indices{pcount};
    if isempty(idx)
      continue                                                              % continue if no analysis point in current segment
    end

    fname = files(pcount,1).name;
    
    % loading
    %-----------------------------------------------------------------------
    [VVlat, VVlon, time, HHchirp,VVchirp,HVchirp,VHchirp,~] = readRadar([path2file,fname]);
   
    % calculate driving direction
    % ---------------------------------------------------------------------
    fullDrivingDir = azimuth(VVlat(1:end-1),VVlon(1:end-1),VVlat(2:end),VVlon(2:end),wgs84);  % heading, clockwise from North
   
    % convert time to depth
    %----------------------------------------------------------------------
    z = time.*1e6.*168./2; z= z(1:cropIdx);

    % initialize S (containing all points in Segment 'pcount')
    %----------------------------------------------------------------------
    clear S
    S = struct('idx',idx,'lat',VVlat(idx),'lon',VVlon(idx));
    S.z= z(1:cropIdx);
    S.time = time(1:cropIdx);
    
    % horizontal averaging, # of traces (averaging over 2*trav traces)
    %----------------------------------------------------------------------
    trav = 10;                                                              % averaging over 2*trav traces
    
    % count through all analysis points of Segment 'pcount'
    %----------------------------------------------------------------------
    for I = 1:length(idx)
      i = idx(I);                                                           % trace index of analysis point
      
      % calculate scattering matrix & driving direction, averaged over 2xtrav traces
      % -------------------------------------------------------------------
      [shhs,svvs,shvs,svhs,drivingDirection,drivingSigma] = computeS(HHchirp,VVchirp,HVchirp,VHchirp,i,trav,fullDrivingDir,I);

      % synthesize azimuthal response
      %----------------------------------------------------------------------
      synt= synthesizeResponse(shhs,svvs,shvs,svhs,cropIdx,drivingDirection);
      
      % % calculate power anomalies
      %--------------------------------------------------------------------
      [PrPar, PrPer, ~, ~] = computePowerAnomalies(synt.shh, synt.svv, ...
      synt.shv, synt.svh, z, dzPowerNorm);
      
      % calculate coherence 
      %--------------------------------------------------------------------
      [chhvv,phaseDer] = computePhaseDerivative(synt.shh, synt.svv, epsilonAverage, deltaEpsilon, centerFrequency, lightSpeed, dzPowerNorm, dZ);
      
      % save to Struct
      % -------------------------------------------------------------------
      S.dPhv{I} = PrPer;
      S.dPhh{I} = PrPar;
      S.chhvv{I} = chhvv;
      S.phasder{I} = phaseDer;
      S.drivingDirection{I} = drivingDirection;
      S.drivingSigma{I} = drivingSigma;
      
    end
    
    % saving
    % ---------------------------------------------------------------------
    save(['output/profile',num2str(linenum),'_',num2str(fname(9:16)),'_',num2str(2*trav),'trace_average.mat'],'S')
  
  end
end
%% function definitions

% compute scattering matrix
% -------------------------------------------------------------------------
function [shhs,svvs,shvs,svhs,drivingDirection,drivingSigma] = computeS(HH,VV,HV,VH,i,trav,drivingDir,I)
    if i>2*trav
      if i>length(HH.r(:,1))-trav
         
         drivingDirection = mean(drivingDir(i-2*trav:i));
         drivingSigma = std(drivingDir(i-2*trav:i));

         shhs = HH.r(i-2*trav:i,:)' + 1i*HH.i(i-2*trav:i,:)';
         svvs = VV.r(i-2*trav:i,:)' + 1i*VV.i(i-2*trav:i,:)';  
         svhs = VH.r(i-2*trav:i,:)' + 1i*VH.i(i-2*trav:i,:)';
         shvs = HV.r(i-2*trav:i,:)' + 1i*HV.i(i-2*trav:i,:)';
      else
        drivingDirection = mean(drivingDir(i-trav:i+trav));  
        drivingSigma = std(drivingDir(i-trav:i+trav));                      
        shhs = HH.r(i-trav:i+trav,:)' + 1i*HH.i(i-trav:i+trav,:)';
        svvs = VV.r(i-trav:i+trav,:)' + 1i*VV.i(i-trav:i+trav,:)';  
        svhs = VH.r(i-trav:i+trav,:)' + 1i*VH.i(i-trav:i+trav,:)';
        shvs = HV.r(i-trav:i+trav,:)' + 1i*HV.i(i-trav:i+trav,:)';
      end

    elseif length(drivingDir)<20
         drivingDirection = mean(drivingDir);
         drivingSigma = std(drivingDir);  
        shhs = HH.r' + 1i*HH.i';
        svvs = VV.r' + 1i*VV.i';  
        svhs = VH.r' + 1i*VH.i';
        shvs = HV.r' + 1i*HV.i';
        disp(['I= ',I, 'i= ',i, ', less than 20 traces averaging']);
    else
      drivingDirection = mean(drivingDir(i:i+2*trav));
      drivingSigma = std(drivingDir(i:i+2*trav));
      shhs = HH.r(1:i+2*trav,:)' + 1i*HH.i(1:i+2*trav,:)';
      svvs = VV.r(1:i+2*trav,:)' + 1i*VV.i(1:i+2*trav,:)';  
      svhs = VH.r(1:i+2*trav,:)' + 1i*VH.i(1:i+2*trav,:)';
      shvs = HV.r(1:i+2*trav,:)' + 1i*HV.i(1:i+2*trav,:)';
    end

    % find artefacts and ignore those traces
    % ---------------------------------------------------------------------
    shhs_dB = 20*log10(abs(shhs));
    svvs_dB = 20*log10(abs(svvs));
    svhs_dB = 20*log10(abs(svhs));
    shvs_dB = 20*log10(abs(shvs));
    
    [~,yindex] = find(shhs_dB>150);
    shhs(:,unique(yindex)) = nan;
    clear yindex
    [~,yindex] = find(svvs_dB>150);
    svvs(:,unique(yindex)) = nan;
    clear yindex
    [~,yindex] = find(svhs_dB>150);
    svhs(:,unique(yindex)) = nan;
    clear yindex
    [~,yindex] = find(shvs_dB>150);
    shvs(:,unique(yindex)) = nan;

end

%% synthesize azimuthal response
% -------------------------------------------------------------------------
function synt= synthesizeResponse(shh,svv,shv,svh,n,drivingDir)

theta = drivingDir+90;
 if theta>360
     theta = theta-360;
 end
 thetanew = zeros(1,360);
 VV = zeros(n,360);
 HH = zeros(n,360);
 VH = zeros(n,360);
 HV = zeros(n,360);
  for N = 1:n
    Sn =  [mean(shh(N,:),'omitnan') mean(svh(N,:),'omitnan');mean(shv(N,:),'omitnan'), mean(svv(N,:),'omitnan')] ;
    t = 1;
    for gamma = 1:360
        thetanew(t) = theta + gamma;
        R = [cosd(gamma) sind(gamma); -sind(gamma), cosd(gamma)];
        Srot = R*Sn*R';
        VV(N,t)= Srot(2,2);
        HH(N,t) = Srot(1,1);
        HV(N,t) = Srot(1,2);
        VH(N,t) = Srot(2,1);
        t=t+1;
    end
  end
  
  % rotate relative to true north
  % -----------------------------------------------------------------------
  [~,B]=find(thetanew>360);
  thetanew(B) = thetanew(B)-360;
  [~,BB] = min(thetanew);
  synt.svv = [VV(:,BB:end),VV(:,1:BB-1)];
  synt.shh = [HH(:,BB:end),HH(:,1:BB-1)];
  synt.shv = [HV(:,BB:end),HV(:,1:BB-1)];
  synt.svh = [VH(:,BB:end),VH(:,1:BB-1)];
end
