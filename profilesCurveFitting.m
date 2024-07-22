
%--------------------------------------------------------------------------
function profilesCurveFitting(linenum,mode)
  % linenum = profile number
  % mode == 0: for fitting only 4 depths 800m, 1000m, 1500m, and 2000m
  % mode == 1: for fitting high resolution, at depth interval of 20m

  path2file = 'output/';

  files = dir([path2file,'profile',num2str(linenum),'*trace_average.mat']);
  numSegments = length(files);

  % data loading and set up matrices
  %--------------------------------------------------------------------------

  data = cell(1, numSegments);
  for i = 1:numSegments
      data{i} = load([path2file,files(i,1).name]);
  end

  % Extract latitudes and longitudes
  %--------------------------------------------------------------------------
    
  latvect = cellfun(@(x) x.S.lat, data, 'UniformOutput', false);
  lonvect = cellfun(@(x) x.S.lon, data, 'UniformOutput', false);
  latvect = cat(1, latvect{:});
  lonvect = cat(1, lonvect{:});
    
  % Calculate distances and times
  %--------------------------------------------------------------------------
  time = 2 * data{1}.S.z / 168;
  matHeight = 14189;
   
  % Initialize arrays
  %--------------------------------------------------------------------------
  numAngles = length(latvect);
  numSamples = 360 * numAngles;
  dPHH = zeros(matHeight, numSamples);
  dPHV = zeros(matHeight, numSamples);
  CHHVV = zeros(matHeight, numSamples);
  Phaseder = zeros(matHeight, numSamples);
     
  % Populate arrays
  %--------------------------------------------------------------------------
  offset = 0;
  drivingDir =[];
  drivingSigma =[];
  for i = 1:numel(data)
    numAngles_i = numel(data{i}.S.dPhh);
    for j = 1:numAngles_i
        d = data{i}.S.dPhh{j};
        dPHH(:, offset + (1:360)) = d(1:matHeight,:);
        d = data{i}.S.dPhv{j};
        dPHV(:, offset + (1:360)) = d(1:matHeight,:);
        d = angle(data{i}.S.chhvv{j});
        CHHVV(:, offset + (1:360)) = d(1:matHeight,:);
        d = data{i}.S.phasder{j};
        Phaseder(:, offset + (1:360)) =  d(1:matHeight,:);
        offset = offset + 360;
    end
     drivingDir= [drivingDir, data{1,i}.S.drivingDirection];
     drivingSigma= [drivingSigma, data{1,i}.S.drivingSigma];
  end
    
  %% curve Fitting
  %--------------------------------------------------------------------------
  % Initialize structure
  %--------------------------------------------------------------------------
  Sav = struct('lat', latvect, 'lon', lonvect);
    
  for i=1:length(latvect)
    Sav.drivingDir(i) = drivingDir{1,i};
    Sav.drivingSigma(i) = drivingSigma{1,i};
  end
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

  % Loop over coordinate points 
  %--------------------------------------------------------------------------
  for i = 1:N
     % find dPHH and correct for driving direction
     Mat = dPHH(:,(i-1)*360+1:i*360);
     Sav.dPHH{1,i} = Mat;
    
     % curve fitting
     %----------------------------------------------------------------------
     dpi = 1;
     % loop over depths
     % save dp points
     %----------------------------------------------------------------------
     if mode==0
         Sav.dp = [1095,2286,5262,8238];
     elseif mode==1
         [~,Sav.dp] = min(abs(Sav.depth-(650:20:2550)));
     end
     %----------------------------------------------------------------------
     maxamp180 = zeros(1,length(Sav.dp));
     maxampidx180 = zeros(1,length(Sav.dp));
     maxampDir180 = zeros(1,length(Sav.dp));
     maxamp90 = zeros(1,length(Sav.dp));
     maxampidx90 = zeros(1,length(Sav.dp));
     maxampDir90 = zeros(1,length(Sav.dp));
     rsquare90 = zeros(1,length(Sav.dp));
     rsquare180 = zeros(1,length(Sav.dp));
     sse90 = zeros(1,length(Sav.dp));
     sse180 = zeros(1,length(Sav.dp));
     rmse90 = zeros(1,length(Sav.dp));
     rmse180 = zeros(1,length(Sav.dp));
     for dp=Sav.dp  
         % Sample data
        %------------------------------------------------------------------
        x = 1:360; 
        y = mean(Mat(dp-50:dp+50,:),'omitnan');

        fitcomb = fittype(@(a180, b180, a90, b90, c, x) a180 * sind(2*x + b180) ...
        + a90*sind(4*x + b90)+c, 'independent', 'x');
        [sine_fit_comb, ~] = fit(x', y', fitcomb, 'StartPoint', ...
        [range(y)/2,0, range(y)/2, 0, mean(y)]);

        % combined signal
        % -------------------------------------------------------------------
        curvecomb = sine_fit_comb.a180*sind(2*x+sine_fit_comb.b180)+...
        sine_fit_comb.a90*sind(4*x+sine_fit_comb.b90) + sine_fit_comb.c;

       
        % Fit a sine curve with 180-degree periodicity
        %--------------------------------------------------------------------------
        fit180 = fittype(@(a, b, c, x) a * sind(2*x + b) + c, 'independent', 'x');
        [sine_fit_180, gof180] = fit(x', y', fit180, 'StartPoint', [range(y)/2, 0, mean(y,'omitnan')]);


        % Fit a sine curve with 90-degree periodicity
        %--------------------------------------------------------------------------
        fit90 = fittype(@(a, b, c, x) a * sind(4*x + b) + c, 'independent', 'x');
        [sine_fit_90, gof90] = fit(x', y', fit90, 'StartPoint', [range(y)/2, 0, mean(y,'omitnan')]);

        % Calculate x-location of maximum amplitude
        %--------------------------------------------------------------------------
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
        
    
     % saving  
     % --------------------------------------------------------------------------
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
      
    
     % find and save dPHV
     %--------------------------------------------------------------------------
     Sav.dPHV{1,i} = dPHV(:,(i-1)*360+1:i*360);
      
     % find and save CHHVV
     %--------------------------------------------------------------------------
     Sav.CHHVV{1,i} = CHHVV(:,(i-1)*360+1:i*360);
      
     % find and save Phase der
     %--------------------------------------------------------------------------
     Sav.Phaseder{1,i} = Phaseder(:,(i-1)*360+1:i*360);
      
  end


  % Save the results
  %--------------------------------------------------------------------------
  if mode == 0
      save(['output/curveFit_',num2str(linenum),'_20trace_average.mat'], 'Sav', '-v7.3');
  elseif mode == 1
      save(['output/curveFit_',num2str(linenum),'_20trace_average_highres.mat'], 'Sav', '-v7.3');
  end

end