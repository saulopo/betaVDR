function [xg,yg,zg] = LoadXYZ(FileName);

% Loads 3-column xyz file and converts to gridded data

M = dlmread(FileName);%,',',1,0);
x = M(:,1); y = M(:,2); z = M(:,3);

% detect grid dimensions

if(y(2)==y(1))     % x is the fast index
  
  n = 0;           % number of points on direction y
  yo = y(1);
  while y(n+1)==yo
    n = n + 1;
  end

  m = length(x)/n; % number of points on direction x
  
  xg = reshape(x,n,m);
  yg = reshape(y,n,m); 
  zg = reshape(z,n,m); 
  xg = xg'; yg = yg'; zg = zg';

else % y is the fast index

  m = 0;           % number of points on direction x
  xo = x(1);
  while x(m+1)==xo
    m = m + 1;
  end
  
  n = length(x)/m; % number of points on direction y
  
  xg = reshape(x,m,n);
  yg = reshape(y,m,n); 
  zg = reshape(z,m,n);

end
