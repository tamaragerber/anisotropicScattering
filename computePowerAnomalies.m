% power anomalies
% -------------------------------------------------------------------------
function [PrPar, PrPer, Prvv, Prhv] = computePowerAnomalies(shh, svv, ...
  shv, svh, z, dzPowerNorm)
    
    PrPar = computeAnomaly(shh, z, dzPowerNorm);
    PrPer = computeAnomaly(svh, z, dzPowerNorm);
    Prvv = computeAnomaly(svv, z, dzPowerNorm);
    Prhv = computeAnomaly(shv, z, dzPowerNorm);
end

%% Helper function to compute anomalies
% -------------------------------------------------------------------------
function anomaly = computeAnomaly(data, z, dzPowerNorm)
    anomaly = 20 * log10(abs(data));
    anomaly = detrend(anomaly', 'constant')';
    anomaly = AverageDepth(anomaly, z(2:end), dzPowerNorm);
end