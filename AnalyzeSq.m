tic
% Input/Output files
   %param=load('param');
   %run = param(1);
   %ntime = param(2);
   run = 2189;
   ntime = 10;

% Real space grid inputs
   L = 19.4456;                 % half of sim box size 50000
   %L = 20.664;                  % half of sim box size 60000
   %L = 38.8911;                 % half of sim box size 400000
   Ndiv=121;                     % number of bins (assumed will be odd)

for time = ntime:ntime
%for time = 1:ntime
tic
% Load xyz data into 3D arrays
   loadfile = sprintf('%d_%d.xyz', run, time);
   savefile = sprintf('%d_%d.sk', run, time);
   data=load(loadfile);
   ndata=length(data);
   j=1;
   k=1;
   for i = 1:ndata
      X(i,1:3) = data(i,2:4);
      if data(i,1) == 1
     XA(j,1:3) = data(i,2:4);
         j=j+1;
      else
         XB(k,1:3) = data(i,2:4);
         k=k+1;
      end
   end
toc
tic
% Grid real space for density arrays
   dx =L*2/(Ndiv);                        % size of bin
   x = -L:dx:L;                   % bins
   y = -L+dx/2:dx:L-dx/2;             % midpoints
   [g edges mid loc] = histcn(X,x,x,x);   % n-d histogram
   f = g(1:Ndiv,1:Ndiv,1:Ndiv);
   [gA edges mid loc] = histcn(XA,x,x,x); % n-d histogram A type
   fA = gA(1:Ndiv,1:Ndiv,1:Ndiv);
   [gB edges mid loc] = histcn(XB,x,x,x); % n-d histogram B type
   fB = gB(1:Ndiv,1:Ndiv,1:Ndiv);
toc
tic
% Fourier Transform
   L2=abs(y(Ndiv)-y(1));
   omega=(2*pi/(L2+dx))*[1:1:Ndiv];           % setting correct scale
   omega=omega-omega(ceil(Ndiv/2));           % centering around zero assumes odd N
% Do Fourier transforms with correct centering
   ftk=(fftshift(fftn(fftshift(f))*dx));      % All
   ftkA=(fftshift(fftn(fftshift(fA))*dx));    % A type
   ftkB=(fftshift(fftn(fftshift(fB))*dx));    % B type
   ftkM=(fftshift(fftn(fftshift(fA-fB))*dx)); % Steinmueller
% Calculate structure factures
   sk =  1/ndata * abs(ftk).^2;
   skA = 1/ndata * abs(ftkA).^2;
   skB = 1/ndata * abs(ftkB).^2;
   skM = 1/ndata * abs(ftkM).^2;

   [k1 k2 k3] = meshgrid(omega);         % k-space grid
toc
tic
% Radial Average
   normk=sqrt(k1.*k1+k2.*k2+k3.*k3);
   Nbins=Ndiv*1.5;               % Choose bin size so at least 1 grid point in each bin
   kmax = sqrt(3)*omega(Ndiv);               % Exclude points where full sphere of radius normk not within box
   dk=kmax/Nbins;
   % Bin in k-space
   C=zeros(Nbins+1,6);
   for i=1:Ndiv
      for j=1:Ndiv
         for k=1:Ndiv
        kval=normk(i,j,k);
            bindex = floor(kval/dk) + 1;
        C(bindex,1:6) = C(bindex,1:6) + [1 kval sk(i,j,k) skA(i,j,k) skB(i,j,k) skM(i,j,k)];
         end
      end
   end
toc
% Output radial averaged Sk
   kstr=floor(omega(Ndiv)/dk)
   %D(:,2:6) = C(1:kstr,2:6)./C(1:kstr,1);   % Calculate aveages by dividing through by counts
   D(:,2) = C(1:kstr,2)./C(1:kstr,1);   % Calculate aveages by dividing through by counts
   D(:,3) = C(1:kstr,3)./C(1:kstr,1);   % Calculate aveages by dividing through by counts
   D(:,4) = C(1:kstr,4)./C(1:kstr,1);   % Calculate aveages by dividing through by counts
   D(:,5) = C(1:kstr,5)./C(1:kstr,1);   % Calculate aveages by dividing through by counts
   D(:,6) = C(1:kstr,6)./C(1:kstr,1);   % Calculate aveages by dividing through by counts
   D(:,1) = 0:kstr-1;
   D(:,1) = D(:,1)*dk+dk/2;
   % Output array containing | midpoint of k-space bin | average k for bin | Sk | SkA | SkB | SkM |
   dlmwrite(savefile,D,'\t');
end
toc
% Check if FT looks sensible
% figure
% hold on
% imagesc(omega,omega,ftk(:,:,ceil(Ndiv/3)))
% colorbar('vert')
% hold off

