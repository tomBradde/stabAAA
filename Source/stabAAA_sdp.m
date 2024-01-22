function [Cm] = stabAAA_sdp(cp,A,conType)
% returns a set of barycentric weights that renders the AAA model
% structurally stable

% Compute SVD to perform change of variables
[~,S,V]=svd(A,0);

% Optimal Unconstrained Solution
x=V(:,end);

% Optimal unconstrained solution under image of V' (normalized)
xBar=[zeros(length(x)-1,1);1e-2];

% Denominator input to state matrix
Bm=kron(ones(length(cp),1),[2;0]);

% Denominator A matrix
Am=[];
for ii=1:length(cp)
    blockA=[0,imag(cp(ii));-imag(cp(ii)),0];
    Am=blkdiag(Am,blockA);
end
n=size(A,2);

% Compute state space transformation
invS=diag(1./diag(S));
T=V*S;

% Transformed Denominator state space
Bm=invS*V'*Bm;
Am=invS*V'*Am*T;

% Initialize solver settings (needs Mosek solver)
ops = sdpsettings('solver','mosek','verbose',0);

% normalize Bm
Bm=Bm/norm(Bm);


% Storage function matrix
Y=sdpvar(n);

% Objective for epighraphic representation
r=sdpvar(1);
obj=r;

% Define the selected constraint formulation
if strcmp(conType,'Vector')
    p=sdpvar(n,1);
    X=(Y*Am'+Am*Y) + (Bm*p'+p*Bm')*norm(Am);
    % Define problem constraints
    conditions=[
        Y>=eye(n)*0;
        X<=-eye(n)*0;
        [r,Bm'+xBar'*Y;(Bm'+xBar'*Y)',Y]>=eye(n+1)*0;
        r>=0;
        ];

elseif strcmp(conType,'Nullspace')
    Bnull=null(Bm');
    X=Bnull'*(Y*Am'+Am*Y)*Bnull;

    % Define problem constraints
    conditions=[
        Y>=eye(n)*0;
        X<=-eye(n-1)*0;
        [r,Bm'+xBar'*Y;(Bm'+xBar'*Y)',Y]>=eye(n+1)*0;
        r>=0;
        ];

elseif strcmp(conType,'Scalar')
    k=sdpvar(1);
    X=(Y*Am'+Am*Y)-2*k*(Bm*Bm');

    % Define problem constraints
    conditions=[
        Y>=eye(n)*0;
        X<=-eye(n)*0;
        [r,Bm'+xBar'*Y;(Bm'+xBar'*Y)',Y]>=eye(n+1)*0;
        r>=0;
        ];

end


%Solve the problem
sol1=optimize(conditions,obj,ops);
solverProblems=sol1.problem;
[pFeas,~]=check(conditions);

% print if the solver had problems
if solverProblems~=0
    fprintf('The solver had some problems in computing the optimal solution');
end

% print if the solver violated some constraints
if isempty(find(pFeas(1:2)<=0,1))
else
    fprintf('Some constraints are violated by the retrieved solution');
end

% Retrieve storage matrix
Y=double(Y);

% Compute denominator state to output matrix
Cm=Y\Bm;

% Change representation to obtain real and imaginary parts of weights
Cm=Cm'*invS*V';

% Plot weights of stable model against weights of unconstrained model
figure
plot(x);
hold on
plot(Cm/norm(Cm))
legend('Constrained','Unconstrained')
title('AAA model weights')

end

