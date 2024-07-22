%% clearing
%--------------------------------------------------------------------------
clear all;
close all;
format long;

%% Initializing
% constants
%--------------------------------------------------------------------------
lightSpeed = 299792458;                                                     % Speed of light [m a^-1]

epsilonPerpendicular = 3.15;                                                % Perpendicular component of real part of dielectric permittivity
epsilonParallel = 3.1840;                                                   % Parallel component of real part of dielectric permittivity
epsilonAverage = mean([epsilonParallel, epsilonPerpendicular]);             % mean dielectric permittivity
deltaEpsilon = 0.034;                                                       % Dielectric anisotropy

% other parameters
%--------------------------------------------------------------------------
centerFrequency = 330e6;                                                    % Center frequency [Hz]
azimuthResolution = 360;                                                    % Azimuthal resolution (number of shots)
dzPowerNorm = 50;                                                           % Moving average depth in m
dZ = 0.336;                                                                 % Depth step

% azimuthal increments
%--------------------------------------------------------------------------
x = linspace(0, 2 * pi, azimuthResolution);                                 % azimuthal increments
xdeg = rad2deg(x);                                                          % azimuthal increments in degrees

% coordinate system
%--------------------------------------------------------------------------
wgs84 = wgs84Ellipsoid;                                                     % coordinate system


%%  1) turning circle 
% Read UWB data
%--------------------------------------------------------------------------
data = readUWBData('20220624_concatenated_combined.h5');

% Read antenna GPS data
%--------------------------------------------------------------------------
GPS.port = gpsRead('220624_P_airr_PPP.txt');
GPS.starboard = gpsRead('220624_S_airr_PPP.txt');

%convert timestamp
%--------------------------------------------------------------------------
GPS.port.timestamp = GPS.port.hour/24 + GPS.port.min/60/24 + GPS.port.sec/3600/24;
GPS.starboard.timestamp = GPS.starboard.hour/24 + GPS.starboard.min/60/24 + GPS.starboard.sec/3600/24;

data.HHtimestamp = data.HHdatetime-floor(data.HHdatetime);
data.VVtimestamp = data.VVdatetime-floor(data.VVdatetime);

%interpolate coordinates on timestamp
%--------------------------------------------------------------------------
data.GPSlat_p = interp1(GPS.port.timestamp,GPS.port.latitude,data.HHtimestamp);
data.GPSlon_p = interp1(GPS.port.timestamp,GPS.port.longitude,data.HHtimestamp);

data.GPSlat_s = interp1(GPS.starboard.timestamp,GPS.starboard.latitude,data.HHtimestamp);
data.GPSlon_s = interp1(GPS.starboard.timestamp,GPS.starboard.longitude,data.HHtimestamp);


% indexing turning circle from radargram
%--------------------------------------------------------------------------
idx1 = 2805;
idx2 = 3262;

% time-depth conversion
%--------------------------------------------------------------------------
z = data.HHtime .* 1e6 .* 168 ./ 2;                                         % depth in m

% calculate azimuth of HH polarized traces
%--------------------------------------------------------------------------
az = azimuth(data.GPSlat_p,data.GPSlon_p,data.GPSlat_s,data.GPSlon_s,wgs84);% azimuth of HH polarization (perpendicular to driving direction)
az(az<0) = 360+az(az<0);                                                    % HH azimuth (perpendicular to heading), clockwise from North

azpos =az(idx1:idx2);
azposlin = 1:360;

% Interpolate Radargram
%--------------------------------------------------------------------------
data.HH_chirplin = interpolateRadargram(data.HH_chirp, azpos, azposlin, length(data.HHtime), idx1, idx2);
data.VV_chirplin = interpolateRadargram(data.VV_chirp, azpos, azposlin, length(data.VVtime), idx1, idx2);
data.HV_chirplin = interpolateRadargram(data.HV_chirp, azpos, azposlin, length(data.HHtime), idx1, idx2);
data.VH_chirplin = interpolateRadargram(data.VH_chirp, azpos, azposlin, length(data.VVtime), idx1, idx2);

% Compute Power Anomalies
%--------------------------------------------------------------------------
[circle.shh, circle.svv, circle.shv, circle.svh] = computeS(data.HH_chirplin,data.VV_chirplin,data.HV_chirplin, data.VH_chirplin);
%%
[circle.PrPar, circle.PrPer, circle.Prvv, circle.Prhv] = computePowerAnomalies(circle.shh,...
  circle.svv, circle.shv, circle.svh, z, dzPowerNorm);

% Compute coherence & phase Derivative
%--------------------------------------------------------------------------
[circle.chhvv, circle.phase_der] = computePhaseDerivative(circle.shh, circle.svv, epsilonAverage, deltaEpsilon, centerFrequency, lightSpeed, dzPowerNorm, dZ);

% save time
%--------------------------------------------------------------------------
circle.time = data.HHtime;

save('output/circle.mat','circle')
save('output/data.mat','data')
%% 2) Synthesized Response 
% synthesize azimuthal response from one trace
%--------------------------------------------------------------------------

idx = idx1-100;                                                             % synthesize power response from trace before circle starts
theta = az(idx);                                                            % polarization direction of trace used for synthesizing 
% compute Scattering matrix
% -------------------------------------------------------------------------
[s.shh, s.svv, s.shv, s.svh] = computeS(data.HH_chirp,data.VV_chirp,data.HV_chirp, data.VH_chirp);

% synthesize azimuthal response
% -------------------------------------------------------------------------
synt = synthesizeResponse(s, theta, idx);

% compute power anomalies
%--------------------------------------------------------------------------
[synt.PrPar, synt.PrPer, synt.Prvv, synt.Prhv] = computePowerAnomalies(synt.shh,...
  synt.svv, synt.shv, synt.svh, z(1:11310), dzPowerNorm);

% Compute coherence & phase Derivative
%--------------------------------------------------------------------------
[synt.chhvv, synt.phase_der] = computePhaseDerivative(synt.shh, synt.svv, epsilonAverage, deltaEpsilon, centerFrequency, lightSpeed, dzPowerNorm, dZ);

% save time
%--------------------------------------------------------------------------
synt.time = data.HHtime;

save('output/synt.mat','synt')

%% 3) Fujita Model with EGRIP COF
% read cof data
% -------------------------------------------------------------------------
cof = prepCofInput();

% scattering + birefringence
% -------------------------------------------------------------------------
birsc.depth = cof.depth;
birsc.rxdBs = cof.rxdBs/2;
birsc.rydBs = cof.rydBs/2; 
birsc.exw = cof.exw;
birsc.ex = cof.ex;
birsc.eyw = cof.eyw;
birsc.ey = cof.ey;
birsc.condx = 1e-5*ones(length(birsc.depth),1);
birsc.condy = birsc.condx;

% fujita Model
%--------------------------------------------------------------------------
[model.shh, model.svv, model.shv, model.svh, model.z] = fujitaModel(birsc);

% rotate model domain relative to north:
%--------------------------------------------------------------------------

model.shh = [model.shh(:,end-32:end),model.shh(:,1:end-33)];
model.svv = [model.svv(:,end-32:end),model.svv(:,1:end-33)];
model.shv = [model.shv(:,end-32:end),model.shv(:,1:end-33)];
model.svh = [model.svh(:,end-32:end),model.svh(:,1:end-33)];
%
% calculate power anomalies
%--------------------------------------------------------------------------
[model.PrPar, model.PrPer, model.Prvv, model.Prhv] = computePowerAnomalies(model.shh, model.svv, model.shv, model.svh, model.z(1:end-1), dzPowerNorm);

% Compute coherence & phase Derivative
%--------------------------------------------------------------------------
[model.chhvv, model.phase_der] = computePhaseDerivative(model.shh, model.svv, epsilonAverage, deltaEpsilon, centerFrequency, lightSpeed, dzPowerNorm, dZ);

save('output/model.mat','model')
%% 4) plot Fig. 4

plotFig4(model, synt, circle, xdeg);


%% 5) plot Fig. 5

plotFig5(model, synt, circle, xdeg,z);


%% Function Definitions
% read UWB data
% -------------------------------------------------------------------------
function data = readUWBData(filename)

    % Reads the UWB data from an HDF5 file
    % ---------------------------------------------------------------------
    data.HHdatetime = h5read(['input/',filename], '/concatenated/HH/datetime')/1000/60/60/24;
    data.VVdatetime = h5read(['input/',filename], '/concatenated/VV/datetime')/1000/60/60/24;
    data.HH_chirp = h5read(['input/',filename], '/concatenated/HH/Chirps');
    data.VV_chirp = h5read(['input/',filename], '/concatenated/VV/Chirps');
    data.HV_chirp = h5read(['input/',filename], '/concatenated/HV/Chirps');
    data.VH_chirp = h5read(['input/',filename], '/concatenated/VH/Chirps');
    data.HHlat = h5read(['input/',filename], '/concatenated/HH/lat');
    data.HHlon = h5read(['input/',filename], '/concatenated/HH/lon');
    data.VVlat = h5read(['input/',filename], '/concatenated/VV/lat');
    data.VVlon = h5read(['input/',filename], '/concatenated/VV/lon');
    data.VVdist = h5read(['input/',filename], '/concatenated/VV/distance');
    data.HHdist = h5read(['input/',filename], '/concatenated/HH/distance');
    data.HHtime = h5read(['input/',filename], '/concatenated/HH/_time');
    data.VVtime = h5read(['input/',filename], '/concatenated/VV/_time');
    data.HH_gpstime = h5read(['input/',filename], '/concatenated/HH/_time');
    data.VV_gpstime = h5read(['input/',filename], '/concatenated/VV/_time');
    
end

% read GPS data for turning circle
% -------------------------------------------------------------------------
function gps = gpsRead(gpsname)
   allGps = readmatrix(['input/',gpsname],'NumHeaderLines',17);
   gps.latitude = allGps(:,7);
   gps.longitude = allGps(:,8);
   gps.time = allGps(:,10);
   gps.hour = allGps(:,4);
   gps.min = allGps(:,5);
   gps.sec = allGps(:,6);
end

% interpolate radargram
% -------------------------------------------------------------------------
function chirplin = interpolateRadargram(chirp, azpos, azposlin, timeLength, idx1, idx2)
    chirplin.r = [];
    chirplin.i = [];
    for i = 1:timeLength
        chirplin.r(:, i) = interp1(azpos, chirp.r(idx1:idx2, i), azposlin, 'linear', 'extrap');
        chirplin.i(:, i) = interp1(azpos, chirp.i(idx1:idx2, i), azposlin, 'linear', 'extrap');
    end
end

% compute scattering matrix
% -------------------------------------------------------------------------
function [shh,svv,shv,svh] = computeS(HH,VV,HV,VH)
  shh = HH.r' + 1i*HH.i';
  svv = VV.r' + 1i*VV.i';
  svh = VH.r' + 1i*VH.i';
  shv = HV.r' + 1i*HV.i';
end

% synthesize azimuthal response
% -------------------------------------------------------------------------
function synt= synthesizeResponse(matrix, theta, idx)
  n = 11310;
  VV = zeros(n,360);
  HH = zeros(n,360);
  HV = zeros(n,360);
  VH = zeros(n,360);

  for N = 1:n
    Sn =  [matrix.shh(N,idx) matrix.svh(N,idx);matrix.shv(N,idx), matrix.svv(N,idx)] ;
    thetanew = zeros(1,360);
    for gamma = 1:360
      thetanew(gamma) = theta + gamma;
      R = [cosd((gamma)) sind((gamma)); -sind((gamma)), cosd((gamma))];     % clockwise rotation
      Srot = R*Sn*R';
      VV(N,gamma)= Srot(2,2);
      HH(N,gamma) = Srot(1,1);
      HV(N,gamma) = Srot(1,2);
      VH(N,gamma) = Srot(2,1);
    end
  end
  [~,B]=find(thetanew>360);
  thetanew(B) = thetanew(B)-360;
  [~,BB] = min(thetanew);
  
  synt.svv = [VV(:,BB+1:end),VV(:,1:BB)];
  synt.shh = [HH(:,BB+1:end),HH(:,1:BB)];
  synt.shv = [HV(:,BB+1:end),HV(:,1:BB)];
  synt.svh = [VH(:,BB+1:end),VH(:,1:BB)];
end


% plot Fig. 4
% -------------------------------------------------------------------------
function f1 = plotFig4(model, synt, circle, xdeg)

  f1 = figure();
  
  % Parameters for panel arrangement
  %------------------------------------------------------------------------
  gap = [0.02 0.02];                                                        % [horizontal gap, vertical gap]
  marg_h = [0.12 0.05];                                                     % [bottom margin, top margin]
  marg_w = [0.05 0.05];                                                     % [left margin, right margin], increased left margin

  % Create tight subplots
  % -----------------------------------------------------------------------
  ha = tight_subplot(3, 4, gap, marg_h, marg_w);

  % turning circle
  % -----------------------------------------------------------------------
  cutoff = 11310;                                                           % crop depth
  
  axes(ha(1));
  imagesc(xdeg,circle.time(1:cutoff).*1e6,imgaussfilt(circle.PrPar(1:cutoff,:),10));
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  clim([-10 10])
  xticks([0 90 180 270 360])
  xticklabels([])
  title('dP_{HH}')
  set(gca,'fontsize',14)
  axis([0 360 circle.time(1).*1e6 circle.time(cutoff).*1e6])
  
  axes(ha(2));
  imagesc(xdeg,circle.time(1:cutoff).*1e6,imgaussfilt(circle.PrPer(1:cutoff,:),10));
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  xticklabels([])
  title('dP_{HV}')
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  yticklabels([])
  clim([-10 10])
  axis([0 360 circle.time(1).*1e6 circle.time(cutoff).*1e6])

  axes(ha(3));
  imagesc(xdeg,circle.time(1:cutoff).*1e6,imgaussfilt(angle(circle.chhvv(1:cutoff,:)),10));
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  xticklabels([])
  title('\phi_{HHVV}')
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  yticklabels([])
  axis([0 360 circle.time(1).*1e6 circle.time(11310).*1e6])

  axes(ha(4));
  imagesc(xdeg,circle.time(1:cutoff).*1e6,imgaussfilt(circle.phase_der(1:cutoff,:),10))
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  title('\psi_{HHVV}')
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  clim([-0.2 0.2])
  xticklabels([])
  yticklabels([])
  axis([0 360 circle.time(1).*1e6 circle.time(11310).*1e6])

  % synthesized
  % -----------------------------------------------------------------------
  axes(ha(5));
  imagesc(xdeg,synt.time(1:cutoff).*1e6,imgaussfilt(synt.PrPar,10)) 
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  clim([-10 10])
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  axis([0 360 synt.time(1).*1e6 synt.time(cutoff).*1e6])
  xticklabels([])

  axes(ha(6));
  imagesc(xdeg,synt.time(1:cutoff).*1e6,imgaussfilt(synt.PrPer,10));hold on
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  yticklabels([])
  xticklabels([])
  clim([-10 10])
  axis([0 360 synt.time(1).*1e6 synt.time(cutoff).*1e6])

  axes(ha(7));
  imagesc(xdeg,synt.time(1:cutoff).*1e6,imgaussfilt(angle(synt.chhvv),10));
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  set(gca,'fontsize',14)
  yticklabels([])
  xticklabels([])
  xticks([0 90 180 270 360])
  axis([0 360 synt.time(1).*1e6 synt.time(cutoff).*1e6])
  
  axes(ha(8));
  imagesc(xdeg,synt.time(1:cutoff).*1e6,imgaussfilt(synt.phase_der,10));hold on
  yticks([11.71, 17.63, 23.55, 29.47])
  yticklabels([1000, 1500, 2000, 2500])
  clim([-0.2 0.2])
  yticklabels([])
  xticklabels([])
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  axis([0 360 synt.time(1).*1e6 synt.time(cutoff).*1e6])

  % modeled
  % -----------------------------------------------------------------------
  axes(ha(9));
  imagesc(xdeg,model.z,imgaussfilt(model.PrPar,10)) 
  c=colorbar('location','southoutside');
  c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4)];
  c.Label.String = '[dB]';
  clim([-10 10])
  xticks([0 90 180 270 360])
  set(gca,'fontsize',14)
  axis([0 360 885 model.z(end)]);

  axes(ha(10));
  imagesc(xdeg,model.z,imgaussfilt(model.PrPer,10));
  c= colorbar('location','southoutside');
  c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4)];
  c.Label.String = '[dB]';
  xticks([0 90 180 270 360])
  clim([-10 10])
  yticklabels([])
  set(gca,'fontsize',14)
  axis([0 360 885 model.z(end)]);

  axes(ha(11));
  imagesc(xdeg,model.z,angle(model.chhvv));
  c=colorbar('location','southoutside','Ticks',[-3.14,0,3.14]);
  c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4)];
  clim([-pi pi])
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  yticklabels([])
  axis([0 360 885 model.z(end)]);
  
  axes(ha(12));
  imagesc(xdeg,model.z,(model.phase_der))
  clim([-1 1])
  set(gca,'fontsize',14)
  xticks([0 90 180 270 360])
  axis([0 360 885 model.z(end)]);
  yticklabels([])
  c= colorbar('location','southoutside');
  c.Position = [c.Position(1), 0.05,c.Position(3),c.Position(4)];
  colormap((brewermap([],"RdBu")));

end

% plot Fig. 5
% -------------------------------------------------------------------------
function  f2 = plotFig5(model, synt, circle, xdeg,z)

  f2=figure();
  
  gap = [0.02 0.03];                                                        % [horizontal gap, vertical gap]
  marg_h = [0.05 0.05];                                                     % [bottom margin, top margin]
  marg_w = [0.05 0.04];                                                     % [left margin, right margin], increased left margin

  % Create tight subplots
  % -----------------------------------------------------------------------
  hb = tight_subplot(5, 3, gap, marg_h, marg_w);
  
  j= 1;
  p=1;
  zmod = model.z(1:end-1);

  for i = 500:500:2500
    try
      % modelled dPhh
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      interv = find(zmod<i); interv=find(zmod(interv)>j);
      y = mean(model.PrPar(interv,:));

      axes(hb(3+(p-1)*3));
      plot(xdeg, y,'k.');hold on  

      % curvefitting
      % -------------------------------------------------------------------
      fitcomb = fittype(@(a180, b180, a90, b90, c, xdeg) a180 * sind(2*xdeg + b180) ...
        + a90*sind(4*xdeg + b90)+c, 'independent', 'xdeg');
      [sine_fit_comb, ~] = fit(xdeg', y', fitcomb, 'StartPoint', ...
        [range(y)/2,0, range(y)/2, 0, mean(y)]);

      % combined signal
      % -------------------------------------------------------------------
      curvecomb = sine_fit_comb.a180*sind(2*xdeg+sine_fit_comb.b180)+...
        sine_fit_comb.a90*sind(4*xdeg+sine_fit_comb.b90) + sine_fit_comb.c;
      plot(xdeg,curvecomb,'-','linewidth',4,'color',[0.6,0.6,0.6,0.6])
      
      % 180-periodic signal
      % -------------------------------------------------------------------
      comb180 = sine_fit_comb.a180*sind(2*xdeg+sine_fit_comb.b180)+ sine_fit_comb.c;
      l1=plot(xdeg,comb180, 'r','linewidth',2);hold on
      
      % 90-periodic signal
      % -------------------------------------------------------------------
      comb90 = sine_fit_comb.a90*sind(4*xdeg+sine_fit_comb.b90)+ sine_fit_comb.c;
      l2=plot(xdeg,comb90,'b','linewidth',2);


      legend([l1,l2],['A_{180}=',num2str(round(abs(sine_fit_comb.a180),2))],...
        ['A_{90}=',num2str(round(abs(sine_fit_comb.a90),2))],...
     'location','northoutside','orientation','horizontal')
      legend('boxoff')

      axis([0 360 -9 9])
      xticks([0 90 180 270 360])
      if 3+(p-1)*3 ~=3
        yticklabels([])
      end   

      if 3+(p-1)*3 ~=9
        set(gca,'Xticklabel',[])
      end

      set(gca,'fontsize',14)
      box on
      grid on
    catch
      delete((hb(3+(p-1)*3)));    
    end
  
    try
      % synthetic dPhh
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      intervr  = find(z<i); interv=find(z(interv)>j);
      y = mean(synt.PrPar(intervr,:));

      axes(hb(2+(p-1)*3));
      plot(xdeg, y,'k.');hold on  
        
      % curvefitting
      % -------------------------------------------------------------------
      fitcomb = fittype(@(a180, b180, a90, b90, c, xdeg) a180 * sind(2*xdeg + b180)...
        + a90*sind(4*xdeg + b90)+c, 'independent', 'xdeg');
      [sine_fit_comb, ~] = fit(xdeg', y', fitcomb, 'StartPoint', ...
        [range(y)/2,0, range(y)/2, 0, mean(y)]);

      % combined signal
      % -------------------------------------------------------------------
      curvecomb = sine_fit_comb.a180*sind(2*xdeg+sine_fit_comb.b180)+...
        sine_fit_comb.a90*sind(4*xdeg+sine_fit_comb.b90) + sine_fit_comb.c;
      plot(xdeg,curvecomb,'-','linewidth',4,'color',[0.6,0.6,0.6,0.6])
      
      % 180-periodic signal
      % -------------------------------------------------------------------
      comb180 = sine_fit_comb.a180*sind(2*xdeg+sine_fit_comb.b180)+ sine_fit_comb.c;
      l1=plot(xdeg,comb180, 'r','linewidth',2);hold on

      % 90-periodic signal
      % -------------------------------------------------------------------
      comb90 = sine_fit_comb.a90*sind(4*xdeg+sine_fit_comb.b90)+ sine_fit_comb.c;
      l2=plot(xdeg,comb90,'b','linewidth',2);

      legend([l1,l2],['A_{180}=',num2str(round(abs(sine_fit_comb.a180),2))],...
        ['A_{90}=',num2str(round(abs(sine_fit_comb.a90),2))],...
     'location','northoutside','orientation','horizontal')
      legend('boxoff')

      axis([0 360 -9 9])
      xticks([0 90 180 270 360])
      set(gca,'fontsize',14)
      
      if 2+(p-1)*3 ~=14
        set(gca,'Xticklabel',[])
      end
      yticklabels([])
      box on
      grid on
  
    catch
      delete((hb(2+(p-1)*3)));
    end
  
    try
      % observed dPhh
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      intervr  = find(z<i); interv=find(z(interv)>j);
      y = mean(circle.PrPar(intervr,1:end));

      axes(hb(1+(p-1)*3));
      plot(xdeg, y,'k.');hold on   

      % curve-fitting
      % -------------------------------------------------------------------
      fitcomb = fittype(@(a180, b180, a90, b90, c, xdeg) a180 * sind(2*xdeg + b180)...
        + a90*sind(4*xdeg + b90)+c, 'independent', 'xdeg');
      [sine_fit_comb, ~] = fit(xdeg', y', fitcomb, 'StartPoint', ...
        [range(y)/2,0, range(y)/2, 0, mean(y)]);

      % combined curve
      % -------------------------------------------------------------------
      curvecomb = sine_fit_comb.a180*sind(2*xdeg+sine_fit_comb.b180)+...
        sine_fit_comb.a90*sind(4*xdeg+sine_fit_comb.b90) + sine_fit_comb.c;
      plot(xdeg,curvecomb,'-','linewidth',4,'color',[0.6,0.6,0.6,0.6])

      % 180-periodic signal
      % -------------------------------------------------------------------
      comb180 = sine_fit_comb.a180*sind(2*xdeg+sine_fit_comb.b180)+ sine_fit_comb.c;
      l1=plot(xdeg,comb180, 'r','linewidth',2);hold on

      % 90-periodic signal
      % -------------------------------------------------------------------
      comb90 = sine_fit_comb.a90*sind(4*xdeg+sine_fit_comb.b90)+ sine_fit_comb.c;
      l2=plot(xdeg,comb90,'b','linewidth',2);

      legend([l1,l2],['A_{180}=',num2str(round(abs(sine_fit_comb.a180),2))],...
        ['A_{90}=',num2str(round(abs(sine_fit_comb.a90),2))],...
       'location','northoutside','orientation','horizontal')
      legend('boxoff')
      
      axis([0 360 -9 9])
      xticks([0 90 180 270 360])
      if 1+(p-1)*3 ~=13
        set(gca,'Xticklabel',[])
      end
      set(gca,'fontsize',14)
      box on
      grid on

    catch
      delete((hb(1+(p-1)*3)));
    end
    
  j = i;
  p = p+1;
  end
end