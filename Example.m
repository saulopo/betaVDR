clc;close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthetic example that illustrates the method betaVDR
% from the manuscript "A Stable Finite Difference Method Based
% on Upward Continuation to Evaluate Vertical Derivatives of
% Potential Field Data", published in Pure and Applied Geophysics
%
% https://doi.org/10.1007/s00024-022-03164-z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the experiment

Method = {'Freq','BFD','TN','\beta VDR'};

[x,y,Mo]=LoadXYZ('input.dat');
dx=1; dy=1;

NoiseLevel = 2/100; % 2% noise (change as needed)

noise = randn(size(Mo)); noise=noise/max(abs(noise(:)));
M = Mo + NoiseLevel*max(Mo(:))*noise;

[ny,nx] = size(M);
[Mx,My] = gradient(M,dx,dy);
Mh = sqrt(Mx.^2+My.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute vertical derivatives

Mzv = zeros(ny,nx,4); % save derivatives in a multi-index array

beta = 50;
h = 2*dx;

Mzv(:,:,1) = dzFreq(M,dx,dy,1);     % 1. in frequency domain
Mzv(:,:,2) = BDF(M,dx,dy,h);        % 2. Backward finite diff
Mzv(:,:,3) = TN4th(M,dx,dy,0,h);    % 3. Tran & Nguyen 4th order
Mzv(:,:,4) = betaVDR(M,dx,dy,beta); % 4. beta VDR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute filters

TDX = 0*Mzv; ED = TDX;
for i = 1:4
	Mz = Mzv(:,:,i);
	TDX(:,:,i) = atan(Mh./abs(Mz));
	[Axz,Ayz]  = gradient((sign(Mz)+1)/2,1,1);
	ED(:,:,i)  = sqrt(Axz.^2+Ayz.^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures

FontAxes = 14; FontColorbar = 12;
load ColorMap.mat

figure % Vertical derivative
for i = 1:4 
	subplot(2,2,i)
	pcolor(x,y,Mzv(:,:,i)); shading interp; axis image;
	colormap(ColorMap); caxis([min(Mzv(:)),max(Mzv(:))]);
	set(gca,'FontSize',FontAxes); 
	cb=colorbar; set(cb,'FontSize',FontColorbar);
	tit = get(cb,'title'); 	set(tit,'String','nT/km');
	xlabel('Easting (km)');	ylabel('Northing (km)');
	%grid off; box on; set(gca,'layer','top') % use this line on matlab
	title(['TDX/',char(Method(i))])
end
set(gcf,'position',[33    68   980   593]);

figure % TDX
for i = 1:4 
	subplot(2,2,i)
	pcolor(x,y,TDX(:,:,i)); shading interp; axis image;
	colormap(ColorMap); caxis([0,pi/2]);
	set(gca,'FontSize',FontAxes);
	cb=colorbar; set(cb,'FontSize',FontColorbar);
	tit = get(cb,'title'); 	set(tit,'String','rad');
	xlabel('Easting (km)');	ylabel('Northing (km)');
	%grid off; box on; set(gca,'layer','top') % use this line on matlab
	title(['TDX/',char(Method(i))])
end
set(gcf,'position',[33    68   980   593]);

figure; % ED
for i = 1:4
	subplot(2,2,i)
	pcolor(x,y,ED(:,:,i)); shading interp; axis image; 
	colormap(ColorMap); caxis([0,1])
	set(gca,'FontSize',FontAxes);
	cb=colorbar; set(cb,'FontSize',FontColorbar);
	xlabel('Easting (km)');	ylabel('Northing (km)');
	%grid off; box on; set(gca,'layer','top') % use this line on matlab
	title(['ED/',char(Method(i))])
end
set(gcf,'position',[33    68   980   593]);
