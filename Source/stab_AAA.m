function [r,om,ff,w,errvec,pol] = stab_AAA(F,Om,tol,mmax,conType)
% aaa rational approximation of data F on set 1j*Om
% [r,om,fu,w,errvec,pol] = stab_AAA(F,Om,tol,mmax,conType) %
% Input: F = vector of data values
%        Om = vector of sample points (nonnegative angular frequencies)
%        tol = relative tolerance tol, set to 1e-13 if omitted
%        mmax: Maximum number of iterations set to 100 if omitted
%        conType: Equivalent constraints formulations for stability
%        enforcement based on Finsler's Lemma. Arguments can be strings of values
%        - 'Scalar': constraints include a scalar multiplier (faster but
%        sometimes numerically less reliable
%        - 'Vector': constraints include an instrumental vector of
%        variables (numerically more robust)
%        - 'Nullspace': stability constraints are projected over a smaller
%        subspace. No instrumental variables are needed (very robust. Yet
%        solver is observed to handle this less efficiently).
%
% Output: r = AAA approximant to F (function handle)
%         r = model response
%         om = barycentric nodes
%         fu = target function values at nodes
%         w = barycentric weights
%         errvec = Least squares residual error
%         pol = model poles

% number of sample points
M = length(Om);

% default relative tol 1e-13
if nargin<3
    tol = 1e-3;
end

% default max type (99,99)
if nargin<4
    mmax = 100;
end

if nargin<5
    conType='Nullspace';
end

% work with column vectors
Om=Om(:);
F=F(:);

% left scaling matrices
SFr = spdiags(real(F),0,M,M);
SFi = spdiags(imag(F),0,M,M);

% initializations
J = 1:M; % indices of non-interpolating points
om = []; % support points
fr = []; fi = []; % data values at support points (real and imag)
Cp = []; Cm = []; % Cauchy matrices
errvec = []; % error evolution through iterations
R = mean(F); % start with mean of F as initial approximant

% main loop: add one point at the time
for m = 1:mmax
    
    % select next support point: where model vs data error is maximum
    [~,j] = max(abs(F-R));
    
    %%% Must ensure here that Om(j) is not zero %%%
    
    % update support points and data values
    om = [om; Om(j)]; %#ok<*AGROW>
    fr = [fr; real(F(j))];
    fi = [fi; imag(F(j))];
    
    % update index vector: this removes current sample index j from J,
    % equivalent to J = setdiff(J,j)
    J(J==j) = [];
    
    % add column related to current support point to Cauchy matrices
    Cp = [Cp 1./(Om-Om(j))+1./(Om+Om(j))];
    Cm = [Cm 1./(Om-Om(j))-1./(Om+Om(j))];
    
    % right scaling matrices
    Sfr = diag(fr);
    Sfi = diag(fi);
    
    % Loewner matrix
    A=[ SFi(J,J)*Cp(J,:)-Cm(J,:)*Sfi, SFr(J,J)*Cm(J,:)-Cm(J,:)*Sfr;
        -SFr(J,J)*Cp(J,:)+Cp(J,:)*Sfr, SFi(J,J)*Cm(J,:)-Cp(J,:)*Sfi];
    
    % Shape coeffs alternating real and imag (not efficient but easier)
    Ar=A(:,1:m);
    Ai=A(:,m+1:end);
    A2=[];
    for ii=1:size(Ar,2)
        A2=[A2,Ar(:,ii),Ai(:,ii)];
    end

    % Compute homogeneous Least squares optimal solution
    [~,~,V2] = svd(A2,0);    
    Cden=V2(:,end);
    
    % Collect real (wr) and imag (wi) weight components
    wr = Cden(1:2:end-1);
    wi = Cden(2:2:end);
    
    % Shape them
    wt=[wr;wi];

    % Compute the numerator and denominator of the barycentric interpolant
    N = 1j*Cp*(wi.*fi - wr.*fr) + Cm*(wr.*fi + wi.*fr);
    D = [-1j*Cp, Cm]*wt;
    
    % rational approximation
    R = F; % Initialize by data, so that the support points are correct
    R(J) = N(J)./D(J); % replace only the LS-computed points
    
    % error
    err = norm(F-R,inf);
    
    % max error at sample points: add to cumulative error history
    errvec = [errvec; err];
    
    
    
    
    % At last iteration, enforce stability
    if err <= tol*norm(F,inf)
        
            % Barycentric nodes
            cp=1i*om;
          
            % Stability Enforcement
            [Cden] = stabAAA_sdp(cp,A2,conType);
            
            % Make column and normalize
            Cden=Cden'/norm(Cden);
          
            % real and imaginary parts of the weights
            wr = Cden(1:2:end-1);
            wi = Cden(2:2:end);
            wt=[wr;wi];

            % numerator and denominator
            N = 1j*Cp*(wi.*fi - wr.*fr) + Cm*(wr.*fi + wi.*fr);
            D = [-1j*Cp, Cm]*wt;
            
            % rational approximation
            R = F; % Initialize by data, so that the support points are correct
            R(J) = N(J)./D(J); % replace only the LS-computed points
            
            % error
            err = norm(F-R,inf);
            
            % max error at sample points: add to cumulative error history
            errvec = [errvec; err];
            break;

    end
end

% Output arguments
w = [wr wi];

% Construction of interpolant from barycentric form (manual...)
CCp = -1j./bsxfun(@minus,Om,om.');
CCm = -1j./bsxfun(@minus,Om,-om.');
ww = wr+1j*wi;
ff = fr+1j*fi;
r = ( CCp*(ww.*ff) + CCm*conj(ww.*ff) )./( CCp*(ww) + CCm*conj(ww) );

% Exact interpolation at support points.
indices=isnan(r);
r(indices)=F(indices);

% poles, residues, and zeros
B = eye(2*m+1);
B(end,end) = 0;
A = [zeros(m,m), diag(om), 2*ones(m,1);
    -diag(om), zeros(m,m+1);
    w(:).' 0];

% poles 
pol = eig(A,B);
pol = pol(~isinf(pol));


