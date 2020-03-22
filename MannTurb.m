function [PHI,k2,k3,k11,k2_log,k3_log] = MannTurb(alphaEps,GAMMA,L,varargin)
% [PHI,k2,k3,k11,k2_log,k3_log] = MannTurb(alphaEps,GAMMA,L,varargin)
% compute the two sided Mann spectral tensor for wind turbulence.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alphaEps = 1st constant of Mann spectral tensor [1 x 1]
% GAMMA = shear constant (2nd constant) [1 x 1]
% L = Integral length scale [1 x 1]
% varargin: It can be;
%  - N1: Number of points in the along-wind direction [1 x 1]
%  - N2: Number of points in the across-wind direction [1 x 1]
%  - N3: Number of points in the vertical-wind direction [1 x 1]
%  - k1min: min value of wavenumber for k1 [1 x 1]
%  - k2min: min value of wavenumber for k2 [1 x 1]
%  - k3min: min value of wavenumber for k3 [1 x 1]
%  - k1max = max value of wavenumber for k1 [1 x 1]
%  - k2max: max value of wavenumber for k2 [1 x 1]
%  - k3max: max value of wavenumber for k3 [1 x 1]
%  - Ninterp: Number of interpolation points for 2F1 approximation [ 1 x 1]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHI : 5D Spectral tensor for the three wind components [2N1 x 2N2 x 2N3 x 3 x 3]
% k2 : 2-sided wavenumber [2 x N1, 2 x N2, 2 x N3]
% k3 : 2-sided wavenumber [2 x N1, 2 x N2, 2 x N3]
% k11: : Single-sided wavenumber [1,1 x N1]
% k2_log : two sided wavenumber [1, 2 x N2]
% k3_log : two sided wavenumber [1, 2 x N3]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTHAX:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PHI,k2,k3,k11,k2_log,k3_log] = MannTurb(alphaEps,GAMMA,L)
% [PHI,k2,k3,k11,k2_log,k3_log] = MannTurb(alphaEps,GAMMA,L,'N1',50)
% [PHI,k2,k3,k11,k2_log,k3_log] = MannTurb(alphaEps,GAMMA,L,'Ninterp',100)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Etienne Cheynet - last modified : 16/04/2015
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see also MannCoherence.m fitMannTensor.m

%% INPUT parser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('N1',100); % 
p.addOptional('N2',50); % 
p.addOptional('N3',50); % 
p.addOptional('k1min',-4); % 
p.addOptional('k2min',-4); % 
p.addOptional('k3min',-4); % 
p.addOptional('k1max',log10(10)); % 
p.addOptional('k2max',log10(10)); % max value of wavenumber for k2
p.addOptional('k3max',log10(10)); % 
p.addOptional('Ninterp',100); % 
p.parse(varargin{:});
% check number of input: Number of outputs must be >=3 and <=6.
nargoutchk(3,6)
% shorthen the variables name
N1 = p.Results.N1 ;
N2 = p.Results.N2 ;
N3 = p.Results.N3 ;
k1min = p.Results.k1min ;
k2min = p.Results.k2min ;
k3min = p.Results.k3min ;
k1max = p.Results.k1max ;
k2max = p.Results.k2max ;
k3max = p.Results.k3max ;
Ninterp = p.Results.Ninterp ;
% cosntruct wavenumber vectors
k1_log=[-fliplr(logspace(k1min,k1max,N1)),logspace(k1min,k1max,N1)];
k2_log=[-fliplr(logspace(k2min,k2max,N2)),logspace(k2min,k2max,N2)];
k3_log=[-fliplr(logspace(k3min,k3max,N3)),logspace(k3min,k3max,N3)];

%%%%%%%%%%%%%%%%%%
% PREALLOCATIONS
%%%%%%%%%%%%%%%%%%
PHI = zeros(2*N1,2*N2,2*N3,3,3); % spectral tensor
% 3D box where the spectral tensor is built
[k2,k1,k3]=meshgrid(k2_log,k1_log,k3_log);
% definition of k
k=sqrt(k1.^2+k2.^2+k3.^2); % [2Nk1 x 2Nk2 x 2Nk3] matrix

if GAMMA<1e-3,
    warning('GAMMA <1e-3 is set to 0')
    ES= VonKarmanIsoTensor(alphaEps,L,k);
    % Diagonal terms
    PHI(:,:,:,1,1)= ES./(4*pi.*k.^4).*(k.^2-k1.^2);
    PHI(:,:,:,2,2)= ES./(4*pi.*k.^4).*(k.^2-k2.^2);
    PHI(:,:,:,3,3)= ES./(4*pi.*k.^4).*(k1.^2+k2.^2);
    % off-diagonal terms
    PHI(:,:,:,1,2)= ES./(4*pi.*k.^4).*(-k1.*k2);
    PHI(:,:,:,2,1)=PHI(:,:,:,1,2);
    PHI(:,:,:,1,3)= ES./(4*pi.*k.^2.*k.^2).*(-k1.*k3);
    PHI(:,:,:,3,1)=PHI(:,:,:,1,3);
    PHI(:,:,:,2,3)= ES./(4*pi.*k.^2.*k.^2).*(-k2.*k3);
    PHI(:,:,:,3,2)=PHI(:,:,:,2,3);
else
    %%%%%%%%%%%%%%%%%%
    % HYPERGEOM TRICK
    %%%%%%%%%%%%%%%%%%
    % Linear interpolation for Hyeprgeom 2F1
    ValToInterp = -(k(:).*L).^(-2);
    x = sort(ValToInterp);
    x = x(1:round(numel(ValToInterp)/Ninterp):end);
    F = griddedInterpolant(x,hypergeom([1/3,17/6],4/3,x));
    HYPERGEOM = F(ValToInterp);
    %     HYPERGEOM = interp1(x,hypergeom([1/3,17/6],4/3,x),ValToInterp);
    be = GAMMA.*(k(:).*L).^(-2/3).*(HYPERGEOM).^(-1/2);
    %      be = GAMMA.*(k(:).*L).^(-2/3);
    be = reshape(be,2*N1,2*N2,2*N3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % k30,k0,C1,C2,xi1 and xi2
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    k30=k3+be.*k1;% definition of k30
    k0=sqrt(k1.^2+k2.^2+k30.^2); % definition of k0
    % CALCULATION OF C1
    A = be.*k1.^2.*(k0.^2-2.*k30.^2+be.*k1.*k30);
    B = k.^2.*(k1.^2+k2.^2);
    C1=A./B;
    % CALCULATION OF C2
    arg1 = be.*k1.*sqrt(k1.^2+k2.^2);
    arg2 = (k0.^2-k30.*k1.*be);
    C2=k2.*k0.^2./((k1.^2+k2.^2).^(3/2)).*atan2(arg1,arg2);
    % CALCULATION OF xi1 and xi2
    xi1= C1-k2./k1.*C2;
    xi2= k2./k1.*C1+C2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECTRAL TENSOR CALCULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % isotropic tensor with k0
    ES= VonKarmanIsoTensor(alphaEps,L,k0);
    % Diagonal terms
    PHI(:,:,:,1,1)= ES./(4*pi.*k0.^4).*(k0.^2-k1.^2-2*k1.*k30.*xi1+(k1.^2+k2.^2).*xi1.^2);
    PHI(:,:,:,2,2)= ES./(4*pi.*k0.^4).*(k0.^2-k2.^2-2*k2.*k30.*xi2+(k1.^2+k2.^2).*xi2.^2);
    PHI(:,:,:,3,3)= ES./(4*pi.*k.^4).*(k1.^2+k2.^2);
    % off-diagonal terms
    PHI(:,:,:,1,2)= ES./(4*pi.*k0.^4).*(-k1.*k2-k1.*k30.*xi2-k2.*k30.*xi1+(k1.^2+k2.^2).*xi1.*xi2);
    PHI(:,:,:,2,1)=PHI(:,:,:,1,2);
    PHI(:,:,:,1,3)= ES./(4*pi.*k.^2.*k0.^2).*(-k1.*k30+(k1.^2+k2.^2).*xi1);
    PHI(:,:,:,3,1)=PHI(:,:,:,1,3);
    PHI(:,:,:,2,3)= ES./(4*pi.*k.^2.*k0.^2).*(-k2.*k30+(k1.^2+k2.^2).*xi2);
    PHI(:,:,:,3,2)=PHI(:,:,:,2,3);
end
% single sided wavenumber for PHI
k11 = k1_log(end-N1+1:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
% FUN 1: VON KARMAN ISOTROPIC SPECTRAL TENSOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Von-Karman spectral tensor (1948)
    function Ek = VonKarmanIsoTensor(alphaEps,L,k)
        Ek =  alphaEps.*L.^(5/3).*(L.*k).^4./((1+(L.*k).^2).^(17/6));
    end
end

