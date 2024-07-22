% Function for averaging depth between two layers
function y=AverageDepth(x,z,dz)
 % Averages data over depth
  %  y = movmean(x, dz, 2, 'Endpoints', 'discard');
  for i=1:length(z)
      iok=find(z>z(i)-dz/2&z<z(i)+dz/2);
      y(i,:)=mean(x(iok,:),1);
      
  end
end