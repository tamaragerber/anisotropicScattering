function saveCSVfile(Sav,linenum)

  % separate by depth
  % -----------------------------------------------------------------------
  for depthIdx = 1:length(Sav.dp)

    % Define coordinates of the starting points (as a map)
    % ---------------------------------------------------------------------
    x_starts = Sav.lon;  
    y_starts = Sav.lat;  

    for j=1:length(x_starts)
      % Define the length of the array
      % -------------------------------------------------------------------
      n = Sav.maxAmp180{1,j};
      array_length180(j) = n(depthIdx);
      m = Sav.maxAmp90{1,j};
      array_length90(j) = m(depthIdx);

      %Define the angle (in degrees) from north
      %--------------------------------------------------------------------
      nn = Sav.maxAmpDir180{1,j}; 
      angles_from_north180(j) = nn(depthIdx);
      mm = Sav.maxAmpDir90{1,j};  
      angles_from_north90(j) = mm(depthIdx);
      
      % curve-fitting statistics
      %--------------------------------------------------------------------
      rmse90n = Sav.rmse90{1,j};rmse90(j)=rmse90n(depthIdx);
      rmse180n = Sav.rmse180{1,j};rmse180(j)=rmse180n(depthIdx);
      sse90n = Sav.sse90{1,j};sse90(j)=sse90n(depthIdx);
      sse180n = Sav.sse180{1,j};sse180(j)=sse180n(depthIdx);
      rsquare90n = Sav.rsquare90{1,j};rsquare90(j)=rsquare90n(depthIdx);
      rsquare180n = Sav.rsquare180{1,j};rsquare180(j)=rsquare180n(depthIdx);
    end

    % Open a file for writing
    %----------------------------------------------------------------------
    filename = ['output/',num2str(linenum),'_dir',num2str(Sav.depth(Sav.dp(depthIdx))),'.csv'];
    fileID = fopen(filename, 'w');
    
    % Write header
    %----------------------------------------------------------------------
    fprintf(fileID, 'x_start,y_start,amp180,angleFromNorth180,drivingDir,drivingSigma,constAmp,amp90,angleFromNorth90,rmse90,rmse180,rsquare90,rsquare180,sse90,sse180\n');
    
    % Loop over each arrow
    %----------------------------------------------------------------------
    for i = 1:numel(angles_from_north180)
        
        % define amplitude and angle from north
        %------------------------------------------------------------------
        amp180 = array_length180(i);
        angleFromNorth180 = angles_from_north180(i);

        amp90 = array_length90(i);
        angleFromNorth90 = angles_from_north90(i);

        drivingDir = Sav.drivingDir(i);
        drivingSigma = Sav.drivingSigma(i);

        % Write coordinates to the file
        %------------------------------------------------------------------
        fprintf(fileID, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', x_starts(i), y_starts(i), ...
          amp180*10, angleFromNorth180,drivingDir,drivingSigma,amp180*0+1,amp90*10,angleFromNorth90,...
          rmse90(i),rmse180(i),rsquare90(i),rsquare180(i),sse90(i),sse180(i));
    end

    % Close the file
    %----------------------------------------------------------------------
    fclose(fileID);
    disp(['Orientation data saved to ', filename]);

  end
end