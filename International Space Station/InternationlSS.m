%% Testing script for AAA
clear 
close all

addpath(genpath('D:\TB\idempar'));

% Add to path Mosek solver (or equivalent)
addpath(genpath('C:\Program Files\Mosek\10.0'));

% Add to path Yalmip
addpath(genpath('D:\TB\YALMIP-master'));

% Load Data
load DataISS

% Normalize Freq Axis
f = fpoints/fpoints(end);

% Normlize response
Fvect = FF/norm(FF,inf);

% Define angular frequency
Om = 2*pi*f;

% Define Error Tolerance
tol =0.0001;

% Define Constraint representation 
conType = 'Vector';


%%%%%%%%%%%%%%%%%%%%%%%% proposed %%%%%%%%%%%%%%%%%%%%%%%
[r,om,fu,w,errvec,pol] = stab_AAA(Fvect,Om,tol,100,conType);


% Visual Comparison
figure
loglog(f,abs(Fvect),'linewidth',5,'color','r')
hold on
loglog(f,abs(r),'-.','linewidth',3,'color','b')
title('Model-Data comparison')
xlabel('Normalized frequency')
ylabel('Magnitude');
legend('Data','stabAAA')
set(gca,'FontSize',24)
grid on
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%% scatter poles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(pol,'ro','MarkerSize',24,'Marker','*');
title('Poles of the stabAAA model')
xlabel('Real Part (normalized)')
ylabel('Imaginary Part (normalized)')
grid on
legend('stabAAA poles');
set(gca,'FontSize',24)

% Compute error metrics
InfNorm_stabAAA=max(abs(Fvect-r));
Norm_stabAAA=norm((Fvect-r));

