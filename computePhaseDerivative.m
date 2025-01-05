% coherence phase and phase derivative
% -------------------------------------------------------------------------
function [chhvv, phase_der] = computePhaseDerivative(shh, svv, ei, dei, fc, co, dzPowerNorm, dZ)
    % Computes phase derivative
    chhvv1 = movsum(shh .* conj(svv), [dzPowerNorm/dZ], 1, 'omitnan');
    chhvv2 = sqrt(movsum(abs(shh).^2, [dzPowerNorm/dZ], 1, 'omitnan'));
    chhvv3 = sqrt(movsum(abs(svv).^2, [dzPowerNorm/dZ], 1, 'omitnan'));
    chhvv = chhvv1 ./ (chhvv2 .* chhvv3);
    
    [~, fy] = gradient(imgaussfilt(angle(chhvv), 1));
    phase_der = 2 * co * sqrt(ei) / (4 * pi * fc * dei) .* fy;
end