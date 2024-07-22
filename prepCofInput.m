function cof = prepCofInput()

  fabric = readmatrix('EGRIP_cAxis.csv');

  % Sort by depth
  %------------------------------------------------------------------------
  cof.depth = fabric(:,5);
  [cof.depth, s] = sort(cof.depth);
  fabric = fabric(s, :);
  [cof.depth, u] = unique(cof.depth);
  fabric = fabric(u, :);

  cof.lambda1 = fabric(:,15);
  cof.lambda1_w = fabric(:,16);
  cof.lambda2 = fabric(:,17);
  cof.lambda2_w = fabric(:,18);
  cof.lambda3 = fabric(:,19);
  cof.lambda3_w = fabric(:,20);
  cof.lambda2 = fabric(:,17);
  cof.lambda2_w = fabric(:,18);
  cof.lambda3 = fabric(:,19);
  cof.lambda3_w = fabric(:,20);
  clear fabric  

  % define conductivity
  %------------------------------------------------------------------------
  cof.condx = 2e-5;
  cof.condy = 2e-5;

  
  % Calculate differences and plot eigenvalues
  % convert lambda1,2,3 to lambdax,y,z
  % -----------------------------------------------------------------------
  cof = preprocessEigenvalues(cof);

  % Calculate and plot power ratio
  % -----------------------------------------------------------------------
  cof = calculatePowerRatio(cof);

end


%% additional function

% preprocessing of COF data
% -------------------------------------------------------------------------
function cof = preprocessEigenvalues(cof)
    
    cof.chidx = 140;                                                        % index of flipping direction
    cof.exw = cof.lambda1_w;
    cof.eyw = [cof.lambda2_w(1:cof.chidx); cof.lambda3_w(cof.chidx+1:end)];
    cof.ezw = [cof.lambda3_w(1:cof.chidx); cof.lambda2_w(cof.chidx+1:end)];
    cof.ex = cof.lambda1;
    cof.ey = [cof.lambda2(1:cof.chidx); cof.lambda3(cof.chidx+1:end)];
    cof.ez = [cof.lambda3(1:cof.chidx); cof.lambda2(cof.chidx+1:end)];

end