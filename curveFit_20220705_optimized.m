%% data loading and set up matrices
%--------------------------------------------------------------------------
clear all
close all

%--------------------------------------------------------------------------
% Load data
files = {'prof20220705_103058_20trace_av.mat', ...
         'prof20220705_110610_20trace_av.mat', ...
         'prof20220705_135130_20trace_av.mat', ...
         'prof20220705_180134_20trace_av.mat'};
       
data = cell(1, numel(files));
for i = 1:numel(files)
    data{i} = load(files{i});
end
%--------------------------------------------------------------------------

% Extract latitudes and longitudes
latvect = cellfun(@(x) x.S.lat, data, 'UniformOutput', false);
lonvect = cellfun(@(x) x.S.lon, data, 'UniformOutput', false);
latvect = cat(1, latvect{:});
lonvect = cat(1, lonvect{:});
%--------------------------------------------------------------------------

% Calculate distances and times
wgs84 = wgs84Ellipsoid;
dist = distance(latvect(1:end-1), lonvect(1:end-1), latvect(2:end), lonvect(2:end), wgs84);
time = 2 * data{1}.S.z / 168;

%--------------------------------------------------------------------------
% % Plot coordinates
% figure()
% plot(lonvect, latvect, '*')
%--------------------------------------------------------------------------

% Initialize arrays
numAngles = length(latvect);
numSamples = 360 * numAngles;
dPHH = zeros(size(time, 1), numSamples);
dPHV = zeros(size(time, 1), numSamples);
CHHVV = zeros(size(time, 1), numSamples);
Phaseder = zeros(size(time, 1), numSamples);
%--------------------------------------------------------------------------

% Populate arrays
offset = 0;
for i = 1:numel(data)
    numAngles_i = numel(data{i}.S.dPhh);
    for j = 1:numAngles_i
        dPHH(:, offset + (1:360)) = data{i}.S.dPhh{j};
        dPHV(:, offset + (1:360)) = data{i}.S.dPhv{j};
        CHHVV(:, offset + (1:360)) = angle(data{i}.S.chhvv{j});
        Phaseder(:, offset + (1:360)) = data{i}.S.phasder{j};
        offset = offset + 360;
    end
end
%--------------------------------------------------------------------------

% Driving direction
drivingDir= [data{1,1}.S.drivingDirection,data{1,2}.S.drivingDirection,...
  data{1,3}.S.drivingDirection,data{1,4}.S.drivingDirection];


%% curve Fitting
%--------------------------------------------------------------------------
% Initialize structure
Sav = struct('lat', latvect, 'lon', lonvect, 'drivingDir', drivingDir);

Q = size(dPHH);
N = Q(2) / 360;
Sav.i = 1:N;
Sav.dPHH = cell(1, N);
Sav.dPHV = cell(1, N);
Sav.Phaseder = cell(1, N);
Sav.CHHVV = cell(1, N);
Sav.curvefit180 = cell(1, N);
Sav.curvefit90 = cell(1, N);
Sav.maxAmp180 = cell(1, N);
Sav.maxAmp90 = cell(1, N);
Sav.maxAmpDir180 = cell(1, N);
Sav.maxAmpDir90 = cell(1, N);
Sav.time = time;
Sav.depth = data{1}.S.z;
%--------------------------------------------------------------------------

% Loop over coordinate points 
for i = 1:N
   % find dPHH and correct for driving direction
    matd=dPHH(:,(i-1)*360+1:i*360);
    azN = Sav.drivingDir(i);
    deg = 1:360;
    Hdir = round(Sav.drivingDir(i)+90);
    if Hdir>360
      Hdir = Hdir-360;
    end
    absOrientation = [Hdir:360,1:Hdir];

    [~,helpidx] = min(absOrientation);
    mat = [matd(:,helpidx:end),matd(:,1:helpidx-1)]; %oriented towards north and rotating clockwise. 

    if i ==13
      bla=0
    end
    Sav.dPHH{1,i} = mat;
    %--------------------------------------------------------------------------

    % curve fitting
    dpi = 1;
    % loop over 6 depths
    % save dp points
    Sav.dp = [1095,2286,5262,8238];
  %--------------------------------------------------------------------------
    for dp=Sav.dp  %[500,1000,2000,5000,8000,10000]

        % Sample data
        x = 1:360; % Assuming data spans 360 degrees
        y = mean(mat(dp-50:dp+50,:));% Your dataset here

        % Fit a sine curve with 180-degree periodicity
        fit180 = fittype(@(a, b, c, x) a * sind(2*x + b) + c, 'independent', 'x');
        [sine_fit_180, gof180] = fit(x', y', fit180, 'StartPoint', [range(y)/2, 0, mean(y)]);


        % Fit a sine curve with 90-degree periodicity
        fit90 = fittype(@(a, b, c, x) a * sind(4*x + b) + c, 'independent', 'x');
        [sine_fit_90, gof90] = fit(x', y', fit90, 'StartPoint', [range(y)/2, 0, mean(y)]);

        % Calculate x-location of maximum amplitude
        curve180 = sine_fit_180.a*sind(2*x+sine_fit_180.b)+sine_fit_180.c;
        [maxamp180(dpi),maxampidx180(dpi)] = max(curve180);
        maxampDir180(dpi) = x(maxampidx180(dpi));

        curve90 = sine_fit_90.a*sind(4*x+sine_fit_90.b)+sine_fit_90.c;
        [maxamp90(dpi),maxampidx90(dpi)] = max(curve90);
        maxampDir90(dpi) = x(maxampidx90(dpi));
        
        rsquare90(dpi) = gof90.rsquare;
        rsquare180(dpi) = gof180.rsquare;
        
        sse90(dpi) = gof90.sse;
        sse180(dpi) = gof180.sse;
        
        rmse90(dpi) = gof90.rmse;
        rmse180(dpi) = gof180.rmse;

        Sav.curvefit180{dpi,i} = curve180;
        Sav.curvefit90{dpi,i} = curve90;
    
        dpi = dpi+1;
    end
    %--------------------------------------------------------------------------

  % saving  
  Sav.maxAmp180{1,i} = maxamp180;
  Sav.maxAmp90{1,i} = maxamp90;
  Sav.maxAmpDir180{1,i} = maxampDir180;
  Sav.maxAmpDir90{1,i} = maxampDir90;
  
  Sav.rmse90{1,i} = rmse90;
  Sav.rmse180{1,i} = rmse180;
  Sav.sse90{1,i} = sse90;
  Sav.sse180{1,i} = sse180;
  Sav.rsquare90{1,i} = rsquare90;
  Sav.rsquare180{1,i} = rsquare180;
  %--------------------------------------------------------------------------

  % find and save dPHV
  matd= dPHV(:,(i-1)*360+1:i*360);
  mat = [matd(:,helpidx:end),matd(:,1:helpidx-1)]; %oriented towards north and rotating clockwise. 
  Sav.dPHV{1,i} = mat;
  %--------------------------------------------------------------------------
  
  % find and save CHHVV
  matd = CHHVV(:,(i-1)*360+1:i*360);
  mat = [matd(:,helpidx:end),matd(:,1:helpidx-1)]; %oriented towards north and rotating clockwise. 

  Sav.CHHVV{1,i} = mat;
  %--------------------------------------------------------------------------
  
  % find and save Phase der
  matd = Phaseder(:,(i-1)*360+1:i*360);
  mat = [matd(:,helpidx:end),matd(:,1:helpidx-1)]; %oriented towards north and rotating clockwise. 
  Sav.Phaseder{1,i} = mat;
  %--------------------------------------------------------------------------
end

% saving
save('curveFit20220705_20trace_av.mat','Sav', '-v7.3')
%--------------------------------------------------------------------------

%% plot profile figure 
%--------------------------------------------------------------------------
clear all
close all

load('curveFit20220705_20trace_av.mat')
%

figure()
set(gcf, 'color', 'none'); 
% Loop over coordinate points 
ii =1;
N = length(Sav.lat);
for i = 1:N
  if mod(i,2)~=0
    continue
  end
  
  % plot dPHH;
  subplot(3,floor(N/2),ii)
  M=Sav.dPHH{1,i};
  imagesc(1:360,Sav.time(1:11214),M(1:11214,:));hold on 
  axis ij
  colormap((brewermap([],"RdBu")));hold on
  % plot fitted curves
  for dpi = 1:length(Sav.dp)
    plot(1:360,-Sav.curvefit180{dpi,i}*0.5+Sav.time(Sav.dp(dpi)),'k-', 'linewidth',2);hold on
    plot(1:360,-Sav.curvefit90{dpi,i}*0.5+Sav.time(Sav.dp(dpi)),'k:','linewidth',2);hold on
  end

  if i>2
    yticklabels('')
    yticks([11.71, 17.63,23.55,29.45])
    xticklabels('')
    xticks([90 180 270 360])
  else
    set(gca,'fontsize',14)
    yticks([11.71, 17.63,23.55,29.45])
    yticklabels([1000 1500 2000 2500])
    xticklabels({'$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'})
    set(gca,'TickLabelInterpreter', 'latex')
    xticks([90 180 270 360])
  end
  box on
  up =max(max(Sav.dPHH{1,i}));
  lo = min(min(Sav.dPHH{1,i}));
  clbnd = min([abs(up),abs(lo)]);
  %caxis([-clbnd clbnd])
  caxis([-4 4])

  %--------------------------------------------------------------------------
  % plot dPHV;  
  subplot(3,floor(N/2),floor(N/2)+ii)
  M = Sav.dPHV{1,i};
  imagesc(1:360,Sav.time(1:11214),M(1:11214,:));hold on 
  colormap((brewermap([],"RdBu")));hold on
  if i>2
    yticklabels('')
    yticks([11.71, 17.63,23.55,29.45])
    xticklabels('')
    xticks([90 180 270 360])
  else
    set(gca,'fontsize',14)
    yticks([11.71, 17.63,23.55,29.45])
    yticklabels([1000 1500 2000 2500])  
    xticklabels({'$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'})
    set(gca,'TickLabelInterpreter', 'latex')
    xticks([90 180 270 360])
  end
 
  box on
  up =max(max(Sav.dPHV{1,i}));
  lo = min(min(Sav.dPHV{1,i}));
  clbnd = min([abs(up),abs(lo)]);
  %caxis([-clbnd clbnd])
  caxis([-4 4])

  %-------------------------------------------------------------------------- 
  
  % plot CHHVV;
  subplot(3,floor(N/2),2*floor(N/2)+ii)
  M = Sav.CHHVV{1,i};
  imagesc(1:360,Sav.time(1:11214),M(1:11214,:));hold on 
  colormap((brewermap([],"RdBu")));hold on
   if i>2
    yticklabels('')
    yticks([11.71, 17.63,23.55,29.45])
    xticklabels('')
    xticks([90 180 270 360])
  else
    set(gca,'fontsize',14)
    yticks([11.71, 17.63,23.55,29.45])
    yticklabels([1000 1500 2000 2500])
    xticklabels({'$\frac{\pi}{2}$', '$\pi$', '$\frac{3\pi}{2}$', '$2\pi$'})
    set(gca,'TickLabelInterpreter', 'latex')
    xticks([90 180 270 360])
  end
  box on
  up =max(max(Sav.CHHVV{1,i}));
  lo = min(min(Sav.CHHVV{1,i}));
  clbnd = min([abs(up),abs(lo)]);
  %caxis([-clbnd clbnd])
  caxis([-4 4])
  %--------------------------------------------------------------------------  
  
  ii = ii+1;

end


export_fig('20220705_transparent','-nocrop', '-transparent', '-png')
%--------------------------------------------------------------------------
%% plot scattering orientation
%--------------------------------------------------------------------------
  % Create a figure
  figure;
% plotting
for i = Sav.i
  
  % Define the angle (in degrees) from north
  angle_from_north180 = Sav.maxAmpDir180{1,i};  
  angle_from_north90 = Sav.maxAmpDir90{1,i}; 
  %--------------------------------------------------------------------------
  
  % Convert angle to radians
  angle_rad180 = deg2rad(90-angle_from_north180);
  angle_rad90 = deg2rad(90-angle_from_north90);
  %--------------------------------------------------------------------------
  
  % Define the length of the array
  array_length180 = Sav.maxAmp180{1,i}*0.01;
  array_length90 = Sav.maxAmp90{1,i}*0.01;
  %--------------------------------------------------------------------------
  
  % Define the coordinates of the starting point
  x_start = Sav.lon(i);
  y_start = Sav.lat(i);
  %--------------------------------------------------------------------------
  
  cmap = jet(length(Sav.dp));
  for depthIdx = 1:length(Sav.dp)
   
    % Calculate the components of the array vector
    array_dx180 = array_length180(depthIdx) * cos(angle_rad180(depthIdx));
    array_dy180 = array_length180(depthIdx) * sin(angle_rad180(depthIdx));
    %--------------------------------------------------------------------------

    % Define the coordinates of the ending point (tip of the arrow)
    x_end_positive = array_dx180;
    y_end_positive = array_dy180;
    x_end_negative = -array_dx180;
    y_end_negative = -array_dy180;
    %--------------------------------------------------------------------------
    
    % Plot the positive arrow
    quiver(x_start, y_start, x_end_positive, y_end_positive, 0, 'Color', cmap(depthIdx,:));  
    hold on; % Keep the same plot to add another arrow
    %--------------------------------------------------------------------------

    % Plot the negative arrow
    quiver(x_start, y_start, x_end_negative, y_end_negative, 0, 'Color',cmap(depthIdx,:));
    %--------------------------------------------------------------------------

  end
end
 
% Set axis equal and add grid
axis equal;
grid on;

% Add labels and title
xlabel('X');
ylabel('Y');
%--------------------------------------------------------------------------

%% save to CSV file

% separate by depth
for depthIdx = 1:length(Sav.dp)

  % Define coordinates of the starting points (as a map)
  x_starts = Sav.lon;  % Example x coordinates
  y_starts = Sav.lat;  % Example y coordinates
  %--------------------------------------------------------------------------

  for j=1:length(x_starts)
    
    % Define the length of the array
    n = Sav.maxAmp180{1,j};
    array_length180(j) = n(depthIdx);
    m = Sav.maxAmp90{1,j};
    array_length90(j) = m(depthIdx);
    %--------------------------------------------------------------------------

    %Define the angle (in degrees) from north
    nn = Sav.maxAmpDir180{1,j}; 
    angles_from_north180(j) = nn(depthIdx);
    mm = Sav.maxAmpDir90{1,j};  
    angles_from_north90(j) = mm(depthIdx);
    %--------------------------------------------------------------------------
     % curve-fitting statistics
    
    rmse90n = Sav.rmse90{1,j};rmse90(j)=rmse90n(depthIdx);
    rmse180n = Sav.rmse180{1,j};rmse180(j)=rmse180n(depthIdx);
    sse90n = Sav.sse90{1,j};sse90(j)=sse90n(depthIdx);
    sse180n = Sav.sse180{1,j};sse180(j)=sse180n(depthIdx);
    rsquare90n = Sav.rsquare90{1,j};rsquare90(j)=rsquare90n(depthIdx);
    rsquare180n = Sav.rsquare180{1,j};rsquare180(j)=rsquare180n(depthIdx);
  end
  
  % POSITIVE vectors
  % Open a file for writing
  filename = ['20220705_dir_pos',num2str(Sav.dp(depthIdx)),'.csv'];
  fileID = fopen(filename, 'w');
  %--------------------------------------------------------------------------
  
  % Write header
  fprintf(fileID, 'x_start,y_start,amp180,angleFromNorth180,drivingDir,constAmp,amp90,angleFromNorth90,rmse90,rmse180,rsquare90,rsquare180,sse90,sse180\n');
  %--------------------------------------------------------------------------
  
  % Loop over each arrow
  for i = 1:numel(angles_from_north180)
      % Convert angle to radians and adjust for compass direction
      angle_rad180 = deg2rad(90 - angles_from_north180(i));
      angle_rad90 = deg2rad(90 - angles_from_north90(i));
      %--------------------------------------------------------------------------

      % Calculate the components of the array vector
      array_dx180 = array_length180 * cos(angle_rad180);
      array_dy180 = array_length180 * sin(angle_rad180);
      
      array_dx90 = array_length90 * cos(angle_rad90);
      array_dy90 = array_length90 * sin(angle_rad90);
      %--------------------------------------------------------------------------

      % Define the coordinates of the ending point (tip of the arrow)
      x_end180 = x_starts(i) + array_dx180;
      y_end180 = y_starts(i) + array_dy180;
      
      x_end90 = x_starts(i) + array_dx90;
      y_end90 = y_starts(i) + array_dy90;
      %--------------------------------------------------------------------------
      
      % define amplitude and angle from north
      amp180 = array_length180(i);
      angleFromNorth180 = angles_from_north180(i);
      
      amp90 = array_length90(i);
      angleFromNorth90 = angles_from_north90(i);
      %--------------------------------------------------------------------------

      drivingDir = Sav.drivingDir(i);
      % Write coordinates to the file
      fprintf(fileID, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', x_starts(i), y_starts(i), ...
        amp180*10, angleFromNorth180,drivingDir,amp180*0+1,amp90*10,angleFromNorth90,...
        rmse90(i),rmse180(i),rsquare90(i),rsquare180(i),sse90(i),sse180(i));
  end
  %--------------------------------------------------------------------------

  % Close the file
  fclose(fileID);
  %--------------------------------------------------------------------------
  %--------------------------------------------------------------------------
  
  % NEGATIVE vectors
  % Open a file for writing
  filename = ['20220705_dir_neg',num2str(Sav.dp(depthIdx)),'.csv'];
  fileID = fopen(filename, 'w');
  %--------------------------------------------------------------------------
  
  % Write header
  fprintf(fileID, 'x_start,y_start,amp180,angleFromNorth180,drivingDir,constAmp,amp90,angleFromNorth90\n');
  %--------------------------------------------------------------------------

  % Loop over each arrow
  for i = 1:numel(angles_from_north180)
      % Convert angle to radians and adjust for compass direction
      angle_rad180 = deg2rad(90 - angles_from_north180(i));
      angle_rad90 = deg2rad(90 - angles_from_north90(i));
      %--------------------------------------------------------------------------

      % Calculate the components of the array vector
      array_dx180 = array_length180 * cos(angle_rad180);
      array_dy180 = array_length180 * sin(angle_rad180);
      
      array_dx90 = array_length90 * cos(angle_rad90);
      array_dy90 = array_length90 * sin(angle_rad90);
      %--------------------------------------------------------------------------

      % Define the coordinates of the ending point (tip of the arrow)
      x_end180 = x_starts(i) + array_dx180;
      y_end180 = y_starts(i) + array_dy180;
      
      x_end90 = x_starts(i) + array_dx90;
      y_end90 = y_starts(i) + array_dy90;
      %--------------------------------------------------------------------------
      
      % define amplitude and angle from north
      amp180 = array_length180(i);
      angleFromNorth180 = angles_from_north180(i)-180;
      
      amp90 = array_length90(i);
      angleFromNorth90 = angles_from_north90(i)-180;
      %--------------------------------------------------------------------------
      drivingDir = Sav.drivingDir(i);
      
      % Write coordinates to the file
      fprintf(fileID, '%f,%f,%f,%f,%f,%f,%f,%f\n', x_starts(i), y_starts(i), ...
        amp180*10, angleFromNorth180,drivingDir,amp180*0+1,amp90*10,angleFromNorth90);
      %--------------------------------------------------------------------------
  end

  % Close the file
  fclose(fileID);
  %--------------------------------------------------------------------------

  disp(['Arrows data saved to ', filename]);
end
%--------------------------------------------------------------------------
