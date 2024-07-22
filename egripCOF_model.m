clear all
close all
format long
%% Initializing
lightSpeed = 299792458;                                                     % Speed of light [m a^-1]

epsilonPerpendicular = 3.15;                                                % Parallel component of real part of dielectric permittivity
epsilonParallel = 3.1840;                                                   % Perpendicular component of real part of dielectric permittivity
epsilonAverage = mean([epsilonParallel, epsilonPerpendicular]);             % mean dielectric permittivity
deltaEpsilon = 0.034;                                                       % Dielectric anisotropy

% other parameters
%---------------------------------------------------------------------------------------------------------------
centerFrequency = 330e6;                                                    % Center frequency [Hz]
azimuthResolution = 360;                                                    % Azimuthal resolution (number of shots)
dzPowerNorm = 50;                                                           % Moving average depth in m
dZ = 0.336;                                                                 % Depth step

% azimuthal increments
%---------------------------------------------------------------------------------------------------------------
x = linspace(0, 2 * pi, azimuthResolution);                                 % azimuthal increments
xdeg = rad2deg(x);                                                          % azimuthal increments in degrees

% coordinate system
%---------------------------------------------------------------------------------------------------------------
wgs84 = wgs84Ellipsoid;                                                     % coordinate system

%% Load and preprocess cof data

cof=prepCofInput();
 
%% Fujita Model
% %---------------------------------------------------------------------------------------------------------------
% birefringence only
%---------------------------------------------------------------------------------------------------------------
bir.depth = cof.depth;
bir.rxdBs = cof.rxdBs;
bir.rydBs = cof.rxdBs;                                                      % assume both are the same
bir.exw = cof.exw;
bir.ex = cof.ex;
bir.eyw = cof.eyw;
bir.ey = cof.ey;
bir.condx = 2e-5*ones(length(bir.depth),1);
bir.condy = bir.condx;
[bir.shh, bir.svv, bir.shv, bir.svh, bir.z] = fujitaModel(bir);
%
% calculate power anomalies
%--------------------------------------------------------------------------
[bir.PrPar, bir.PrPer, bir.Prvv, bir.Prhv] = computePowerAnomalies(bir.shh, bir.svv, bir.shv, bir.svh, bir.z, dzPowerNorm);

%% scattering only
% -------------------------------------------------------------------------
sc.depth = cof.depth;
sc.rxdBs = cof.rxdBs;
sc.rydBs = cof.rydBs; 
sc.exw = cof.exw;
sc.ex = cof.ex;
sc.eyw = cof.exw;
sc.ey = cof.ex;
sc.condx = 2e-5*ones(length(sc.depth),1);
sc.condy = sc.condx;
[sc.shh, sc.svv, sc.shv, sc.svh, sc.z] = fujitaModel(sc);

% calculate power anomalies
%--------------------------------------------------------------------------
[sc.PrPar, sc.PrPer, sc.Prvv, sc.Prhv] = computePowerAnomalies(sc.shh, sc.svv, sc.shv, sc.svh, sc.z, dzPowerNorm);

%% birefringence + scattering
% -------------------------------------------------------------------------
birsc.depth = cof.depth;
birsc.rxdBs = cof.rxdBs;
birsc.rydBs = cof.rydBs; 
birsc.exw = cof.exw;
birsc.ex = cof.ex;
birsc.eyw = cof.eyw;
birsc.ey = cof.ey;
birsc.condx = 2e-5*ones(length(birsc.depth),1);
birsc.condy = birsc.condx;
[birsc.shh, birsc.svv, birsc.shv, birsc.svh, birsc.z] = fujitaModel(birsc);

% calculate power anomalies
%--------------------------------------------------------------------------
[birsc.PrPar, birsc.PrPer, birsc.Prvv, birsc.Prhv] = computePowerAnomalies(birsc.shh, birsc.svv, birsc.shv, birsc.svh, birsc.z, dzPowerNorm);


%% Plotting
close all

plotFig3_part1(cof, bir, sc, birsc)

plotFig3_part2(bir, sc, birsc)



%% function definitions

function plotFig3_part1(cof, bir, sc, birsc)
  figure()
  clf
  
  % Parameters for panel arrangement
  % -----------------------------------------------------------------------
  n = 5;                                                                    % Number of subplots
  gap = [0.02 0.02];                                                        % [horizontal gap, vertical gap]
  marg_h = [0.1 0.05];                                                      % [bottom margin, top margin]
  marg_w = [0.1 0.05];                                                      % [left margin, right margin], increased left margin

  % Create tight subplots
  % -----------------------------------------------------------------------
  ha = tight_subplot(1, n, gap, marg_h, marg_w);

  % cof
  % -----------------------------------------------------------------------
  axes(ha(1));
  plot(cof.ex,cof.depth,'.','color',[0.6,0.6,0.8]);hold on
  plot(cof.exw,cof.depth,'d','markersize',5,'color',[0,0,0],...
     'markerfacecolor',[0.6,0.6,0.8]);hold on
   
  plot(cof.ey,cof.depth,'.','color',[0.9,0.7,0.2]);hold on
  plot(cof.eyw,cof.depth,'v','markersize',5,'color',[0,0,0],...
     'markerfacecolor',[0.9,0.7,0.2]);hold on
   
   plot(cof.ez,cof.depth,'.','color',[0.7,0.8,0.2]);hold on
  plot(cof.ezw,cof.depth,'s','markersize',5,'color',[0,0,0],...
     'markerfacecolor',[0.7,0.8,0.2]);hold on 
  axis ij 
  set(gca,'fontsize', 14)
  ylabel('depth [m]')
  grid on
  axis([0 0.8 0 cof.depth(end)])
  xlabel('eigenvalues')

  
  % reflection ratio
  % -----------------------------------------------------------------------
  axes(ha(2));
  plot((cof.rydB-cof.rxdB),cof.depth,'o','markersize',5,'color',[0,0,0],...
    'markerfacecolor',[1,0.8,0.8]);hold on
  axis ij
  plot((cof.rydBs-cof.rxdBs),cof.depth,'color',[1,0.4,0.5],'linewidth',3);
  axis ij
  grid on
  set(gca,'fontsize',14)
  xlabel('reflection ratio [dB]')
  yticks(0:200:1600)
  yticklabels([])
  axis([-15 15 0 cof.depth(end)])

  % birefringence only
  % -----------------------------------------------------------------------
  axes(ha(3));
  imagesc(1:360,cof.depth,bir.PrPar)
  set(gca,'fontsize',14)
  colormap(brewermap([],"RdBu"));
  xlabel('azimuth angle \theta [^o]')
  xticks([0 180 360]);
  yticklabels([])
  axis([0 360 0 cof.depth(end)])
  clim([-15 15])
  grid on

  
  % scattering only
  % -----------------------------------------------------------------------
  axes(ha(4));
  imagesc(1:360,cof.depth,sc.PrPar)
  set(gca,'fontsize',14)
  colormap(brewermap([],"RdBu"));
  xlabel('azimuth angle \theta [^o]')
  yticklabels([])
  axis([0 360 0 cof.depth(end)])
  xticks([0 180 360]);
  clim([-15 15])
  grid on
  
  % combined birefringence + scattering
  % -----------------------------------------------------------------------
  axes(ha(5));
  imagesc(1:360,cof.depth,birsc.PrPar)
  set(gca,'fontsize',14)
  colormap(brewermap([],"RdBu"));
  xlabel('azimuth angle \theta [^o]')
  yticklabels([])
  axis([0 360 0 cof.depth(end)])
  xticks([0 180 360]);
  grid on
  clim([-15 15])

end


function plotFig3_part2(bir, sc, birsc)

  figure()
  n = 5;                                                                    % Number of subplots
  gap = [0.01 0.02];                                                        % [horizontal gap, vertical gap]
  marg_h = [0.1 0.05];                                                      % [bottom margin, top margin]
  marg_w = [0.1 0.05];                                                      % [left margin, right margin], increased left margin

  % Create tight subplots
  % -----------------------------------------------------------------------
  hb = tight_subplot(n, 1, gap, marg_h, marg_w);

  axes(hb(1))
  plot(1:360,bir.PrPar(272,:),'--','color',[0.6,0.6,0.6],'linewidth',2); hold on
  plot(1:360,sc.PrPar(272,:),':','color',[0.4,0.4,0.4],'linewidth',2); hold on
  plot(1:360,birsc.PrPar(272,:),'-','color',[0.2,0.2,0.2],'linewidth',2); hold on
  axis([1,360,-15 15])
  xticks([0 90 180 270 360])
  xticklabels([])
  yticks([-10 0 10])
  grid on
  set(gca,'fontsize',14)
  
  axes(hb(2))
  plot(1:360,bir.PrPar(391,:),'--','color',[0.6,0.6,0.6],'linewidth',2); hold on
  plot(1:360,sc.PrPar(391,:),':','color',[0.4,0.4,0.4],'linewidth',2); hold on
  plot(1:360,birsc.PrPar(391,:),'-','color',[0.2,0.2,0.2],'linewidth',2); hold on
  axis([1,360,-15 15])
  xticks([0 90 180 270 360])
  xticklabels([])
  yticks([-10 0 10])
  grid on
  set(gca,'fontsize',14)
  
  axes(hb(3))
  plot(1:360,bir.PrPar(506,:),'--','color',[0.6,0.6,0.6],'linewidth',2); hold on
  plot(1:360,sc.PrPar(506,:),':','color',[0.4,0.4,0.4],'linewidth',2); hold on
  plot(1:360,birsc.PrPar(506,:),'-','color',[0.2,0.2,0.2],'linewidth',2); hold on
  axis([1,360,-15 15])
  xticks([0 90 180 270 360])
  xticklabels([])
  yticks([-10 0 10])
  grid on
  set(gca,'fontsize',14)
  
  axes(hb(4))
  plot(1:360,bir.PrPar(605,:),'--','color',[0.6,0.6,0.6],'linewidth',2); hold on
  plot(1:360,sc.PrPar(605,:),':','color',[0.4,0.4,0.4],'linewidth',2); hold on
  plot(1:360,birsc.PrPar(605,:),'-','color',[0.2,0.2,0.2],'linewidth',2); hold on
  axis([1,360,-15 15])
  xticks([0 90 180 270 360])
  xticklabels([])
  yticks([-10 0 10])
  grid on
  set(gca,'fontsize',14)
  
  axes(hb(5))
  plot(1:360,bir.PrPar(731,:),'--','color',[0.6,0.6,0.6],'linewidth',2); hold on
  plot(1:360,sc.PrPar(731,:),':','color',[0.4,0.4,0.4],'linewidth',2); hold on
  plot(1:360,birsc.PrPar(731,:),'-','color',[0.2,0.2,0.2],'linewidth',2); hold on
  axis([1,360,-15 15])
  xticks([0 90 180 270 360])
  yticks([-10 0 10])
  grid on
  set(gca,'fontsize',14)

end
