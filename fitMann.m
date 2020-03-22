function [GAMMA,L,alphaEps] = fitMann(k11,Su,Sv,Sw,Suw,varargin)
% [Gamma,L,alphaEps] =
% fitMann(k11,Su,Sv,Sw,Suw,alphaEps,guess) fits the Mann
% spectral tensor to measured 2-sided single and cross-spectra. The fitting
% procedure can be done using 3 or 2- unknown parameters. If alphaEps is
% known, then L and GAMMA are found using the fitting procedure. If
% alphaEps is unknown, then every prameter is found using the
% fitting procedure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k11: single sided wavenumber in the along wind direction.
% Su, Sv, Sw, Suw : [1 x N1] single-point, single-sided non-normalized 
% wind spectra.
% varargin: it can be;
%  - N2: Number of points in the across-wind direction [1 x 1]
%  - N3: Number of points in the vertical-wind direction [1 x 1]
%  - k2min: min value of wavenumber for k2 [1 x 1]
%  - k3min: min value of wavenumber for k3 [1 x 1]
%  - k2max: max value of wavenumber for k2 [1 x 1]
%  - k3max: max value of wavenumber for k3 [1 x 1]
%  - tolX: tolerance for fitting procedure
%  - tolFun: tolerance for fitting procedure
%  - Ninterp: Number of interpolation points for 2F1 approximation [ 1 x 1]
%  - guess: [GAMMA,L] or [GAMMA,L,alphaEps] is the first guess of the 
% coefficients to be fitted. By default it is a [ 1 x 3 ] vector
%  - alphaEps: if unknown, it is empty by default. Otherwise, it is user's
%  defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alphaEps = 1st constant of Mann spectral tensor [1 x 1]
% GAMMA = shear constant (2nd constant) [1 x 1]
% L = Integral length scale [1 x 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author and version:
%  E Cheynet - UiB - last modified:  27/02/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see also MannTurb.m MannCoherence.m

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT parser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser();
p.CaseSensitive = false;
p.addOptional('alphaEps',[]);
p.addOptional('guess',[3,50,0.5]);
p.addOptional('N1',[]);
p.addOptional('N2',100);
p.addOptional('N3',100);
p.addOptional('k2min',-5); % 
p.addOptional('k3min',-5); % 
p.addOptional('k2max',log10(50)); % max value of wavenumber for k2
p.addOptional('k3max',log10(50)); 
p.addOptional('Ninterp',100);
p.addOptional('tolX',1e-3);
p.addOptional('tolFun',1e-3);
p.parse(varargin{:});
% check number of input: Number of outputs must be >=5 and <=17.
% shorthen the variables name
alphaEps = p.Results.alphaEps;
guess = p.Results.guess ;
N1 = p.Results.N1 ;
N2 = p.Results.N2 ;
N3 = p.Results.N3 ;
k2min = p.Results.k2min ;
k3min = p.Results.k3min ;
k2max = p.Results.k2max ;
k3max = p.Results.k3max ;
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Ninterp = p.Results.Ninterp ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE NUMBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double sided k2 and k3 with logarithmically spaced interval points
k2_log=[-fliplr(logspace(k2min,k2max,N2)),logspace(k2min,k2max,N2)];
k3_log=[-fliplr(logspace(k3min,k3max,N3)),logspace(k3min,k3max,N3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a vector kTot used for nlinfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  kTot is k11 repeated 9 times and stored in a [Ndk1 x 3 x 3] matrix

if isempty(N1),
    N1 = numel(k11);
end
if N1<100,
    warning([' N1 contains only ',num2str(N1),' data points. It may not be enough to provide an accurate fit']);
end
if N1~=numel(k11),
    dummyK = k11;
    k11 = logspace(log10(dummyK(1)),log10(dummyK(end)),N1);
    Su = interp1(dummyK,Su,k11);
    Sv = interp1(dummyK,Sv,k11);
    Sw = interp1(dummyK,Sw,k11);
    Suw = interp1(dummyK,Suw,k11);
end






kTot = zeros(N1,3,3);
for ii=1:3,
    for jj=1:3,
        kTot(:,ii,jj) = k11;
    end
end
kTot = reshape(kTot,[],1); % is [9 x Ndk1,1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D to 1D transformation (numerical trick)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformation of 3D data into 1D
clear S
S(:,1,1)= k11(:).*Su(:);
S(:,2,2)= k11(:).*Sv(:);
S(:,3,3)= k11(:).*Sw(:);
S(:,1,2)= 0;
S(:,1,3)= k11(:).*real(Suw(:));
S(:,2,3)= 0;
S(:,2,1)=S(:,1,2);
S(:,3,1)=S(:,1,3);
S(:,3,2)=S(:,2,3);
S = reshape(S,[],1); % is [9 x Ndk1,1]
S(isnan(S)) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isinf(1./max(abs(Suw(:)))),
    warning('Suw is detected as unknown or negligible')
    options=optimset('Display','iter'); % increase the fitting precision if Suw = 0
else
    options=optimset('TolX',tolX,'TolFun',tolFun,'Display','iter');
end
% Coeff = nlinfit(kTot,S,modelFun,guess,options); % fit Mann turbulence at every step
if ~isempty(alphaEps),
    guess = guess(1:2);
    modelFun2 = @MannTurb2; % transform a nested function into anonymous function
    Coeff = lsqcurvefit(@(para,kTot) modelFun2(para,kTot,alphaEps),guess,kTot,S,[0,1],[6,500],options);
    GAMMA = abs(Coeff(1));
    L = abs(Coeff(2));
else
    modelFun1 = @MannTurb1; % transform a nested function into anonymous function
    if numel(guess)<3,
        error('error: guess must contains [GAMMA,L,alphaEps]');
    else
        Coeff = lsqcurvefit(@(para,kTot) modelFun1(para,kTot),guess,kTot,S,[0,1,0],[6,100,2],options);
    end
    GAMMA = abs(Coeff(1));
    L = abs(Coeff(2));
    alphaEps = abs(Coeff(3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN 1: VON KARMAN ISOTROPIC SPECTRAL TENSOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Von-Karman spectral tensor (1948)
    function Ek = VonKarmanIsoTensor(alphaEps,L,k)
        Ek =  alphaEps.*L.^(5/3).*(L.*k).^4./((1+(L.*k).^2).^(17/6));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN 2: MANN I SPECTRAL TENSOR FOR 3 PARA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etienne model -- should be modified by Lene with her model
    function [FM] = MannTurb1(para,kTot)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GOAL---------------------------------------------
        %  Compute the Mann spectral tensor
        % INPUT---------------------------------------------
        % para = [alphaEps,GAMMA,L]
        % alphaEps = 1st constant of Mann spectral tensor [1 x 1]
        % GAMMA = shear constant (2nd constant) [1 x 1]
        % L = Integral length scale [1 x 1]
        % Ktot: vector rpeviously defined. is used here only because nlinfit requires it. actually is useless.  -> trick
        % OUTPUT---------------------------------------------
        % PHI = 5D Spectral tensor for the three wind components  [2Ndk1 x 2Ndk2 x 2Ndk3 x 3 x 3]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        kTot; % trick: kTot is useless here
        GAMMA0 = para(1);
        L0 = para(2);
        PHI = zeros(2*N1,2*N2,2*N3,3,3); % spectral tensor
        % 3D box where the spectral tensor is built
        k1_log=[-fliplr(k11),k11]; % double sided k1 with logarithmically spaced interval points
        [k2,k1,k3]=meshgrid(k2_log,k1_log,k3_log);
        % definition of k
        k=sqrt(k1.^2+k2.^2+k3.^2); % [2Nk1 x 2Nk2 x 2Nk3] matrix
        %%%%%%%%%%%%%%%%%%
        % HYPERGEOM TRICK
        %%%%%%%%%%%%%%%%%%,
        ValToInterp = -(k(:).*L0).^(-2);
        x = sort(ValToInterp);
        x = x(1:round(numel(ValToInterp)/Ninterp):end);
        F = griddedInterpolant(x,hypergeom([1/3,17/6],4/3,x));
        HYPERG = F(ValToInterp);
        be = GAMMA0.*(k(:).*L0).^(-2/3).*(HYPERG).^(-1/2);
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
        % isotropic tensor with k0
        ES= VonKarmanIsoTensor(para(3),L0,k0);
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
        % ratio
        FM= squeeze(trapz(k3_log,trapz(k2_log,PHI,2),3));
        FM = FM(end-N1+1:end,:,:);
        FM(:,1,1) = k11(:).*FM(:,1,1);
        FM(:,2,2) = k11(:).*FM(:,2,2);
        FM(:,3,3) = k11(:).*FM(:,3,3);
        if isinf(1./max(abs(Suw(:))))
            FM(:,1,3) = 0;
            FM(:,3,1) = 0;
        else
            FM(:,1,3) = k11(:).*FM(:,1,3);
            FM(:,3,1) = k11(:).*FM(:,3,1);
        end
        FM(:,1,2) = 0;
        FM(:,2,1) = 0;
        FM(:,2,3) = 0;
        FM(:,3,2) = 0;
        FM = reshape(FM,[],1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN 3: MANN I SPECTRAL TENSOR FOR 2 PARA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [FM] = MannTurb2(para,kTot,alphaEps)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GOAL---------------------------------------------
        %  Compute the Mann spectral tensor
        % INPUT---------------------------------------------
        % para = [alphaEps,GAMMA,L]
        % alphaEps = 1st constant of Mann spectral tensor [1 x 1]
        % GAMMA = shear constant (2nd constant) [1 x 1]
        % L = Integral length scale [1 x 1]
        % Ktot: vector rpeviously defined. is used here only because nlinfit requires it. actually is useless.  -> trick
        % OUTPUT---------------------------------------------
        % PHI = 5D Spectral tensor for the three wind components  [2Ndk1 x 2Ndk2 x 2Ndk3 x 3 x 3]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        kTot; % trick: Ktot is useless here
        GAMMA0 = para(1);
        L0 = para(2);
        PHI = zeros(2*N1,2*N2,2*N3,3,3); % spectral tensor
        % 3D box where the spectral tensor is built
        k1_log=[-fliplr(k11),k11]; % double sided k1 with logarithmically spaced interval points
        [k2,k1,k3]=meshgrid(k2_log,k1_log,k3_log);
        % definition of k
        k=sqrt(k1.^2+k2.^2+k3.^2); % [2Nk1 x 2Nk2 x 2Nk3] matrix
        
        %%%%%%%%%%%%%%%%%%
        % HYPERGEOM TRICK
        %%%%%%%%%%%%%%%%%%,
        ValToInterp = -(k(:).*L0).^(-2);
        x = sort(ValToInterp);
        x = x(1:round(numel(ValToInterp)/Ninterp):end);
        F = griddedInterpolant(x,hypergeom([1/3,17/6],4/3,x));
        HYPERG = F(ValToInterp);
        be = GAMMA0.*(k(:).*L0).^(-2/3).*(HYPERG).^(-1/2);
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
        % isotropic tensor with k0
        ES= VonKarmanIsoTensor(alphaEps,L0,k0);
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
        % ratio = abs(squeeze((trapz(k1_log,trapz(k3_log,trapz(k2_log,PHI,2),3),1))));
        FM= squeeze(trapz(k3_log,trapz(k2_log,PHI,2),3));
        FM = FM(end-N1+1:end,:,:);
        FM(:,1,1) = k11'.*FM(:,1,1);
        FM(:,2,2) = k11'.*FM(:,2,2);
        FM(:,3,3) = k11'.*FM(:,3,3);
        if isinf(1./max(abs(Suw(:))))
            FM(:,1,3) = 0;
            FM(:,3,1) = 0;
        else
            FM(:,1,3) = k11'.*FM(:,1,3);
            FM(:,3,1) = k11'.*FM(:,3,1);
        end
        FM(:,1,2) = 0;
        FM(:,2,1) = 0;
        FM(:,2,3) = 0;
        FM(:,3,2) = 0;
        FM = reshape(FM,[],1);
    end
end