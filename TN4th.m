function Tz = TN4th(Trans_data,dx,dy,beta,dz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth-order Tran-Nguyen method to compute the
% vertical derivative
%
% Adapted from the code from that accompanies the article by
% Arisoy & Dikmen (2011), "Potensoft: MATLAB-based software for
% potential field data processing, modeling and mapping."
% Computers & Geosciences 37, 935–942


if(nargin<5)
	dz = 0.1*dx;
end

gridsize = size(Trans_data); 
nx = gridsize(1);
ny = gridsize(2);
nmax=max([nx ny]); npts=2^nextpow2(nmax);  
data = Trans_data;
cdiff=floor((npts-ny)/2); rdiff=floor((npts-nx)/2);
if cdiff == 0 && rdiff == 0
nmax = nmax+1; npts=2^nextpow2(nmax);
cdiff=floor((npts-ny)/2); rdiff=floor((npts-nx)/2);
end
% data=data-mean(data(:));
data1=taper2d(data,npts,nx,ny,rdiff,cdiff);
dkx= 2.*pi./((npts-1).*dx);
dky= 2.*pi./((npts-1).*dy);
nyqx= (npts/2)+1;
nyqy= (npts/2)+1; 
for j=1:npts
       if j <= nyqx
         kx(j)=(j-1)*dkx;
       else
         kx(j)=(j-(npts+1))*dkx;
       end 
for i=1:npts
       if i <= nyqy
         ky(i)=(i-1)*dky;
       else
         ky(i)=(i-(npts+1))*dky;
       end 
       k(i,j)=sqrt(kx(j).^2+ky(i).^2);
       if(k(i,j)==0) k(i,j)=0.00000001; end;
end
end
freq=k;
%%
h  = beta*dz;

f=fft2(data1);
UC=exp(-freq*(h)).*f;
UCinv=ifft2(UC);
UC1=real(UCinv);
%
UC=exp(-freq*(h+dz)).*f;
UCinv=ifft2(UC);
UC2=real(UCinv);
%
UC=exp(-freq*(h+2*dz)).*f;
UCinv=ifft2(UC);
UC3=real(UCinv);
%
UC=exp(-freq*(h+3*dz)).*f;
UCinv=ifft2(UC);
UC4=real(UCinv);
%
UC=exp(-freq*(h+4*dz)).*f;
UCinv=ifft2(UC);
UC5=real(UCinv);
%

coef=[25,-48,36,-16,3];

Tz=(coef(1)*UC1+coef(2)*UC2+coef(3)*UC3+coef(4)*UC4+coef(5)*UC5)./(12*dz);
Tz=Tz(1+rdiff:nx+rdiff,1+cdiff:ny+cdiff);


