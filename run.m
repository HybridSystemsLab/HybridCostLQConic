%This script allows one to reproduce Example 1 in Ferrante Sanfelice CDC
%2018

%This code requires YALMIP, a parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. Any SDP solver can be used.    


%Definition of the data
Ac=[0 1;
   -1 0] ;

Ad=expm(Ac-eye(2));

np=max(size(Ac));

Qc=eye(np);
M=[-1, 0.5;0.5 0];

%%Solution to the SDP problem to get the matrix P defining the
%%Lyapunov-like inequality

tau1=sdpvar(1,1,'full');
tau2=sdpvar(1,1,'full');
lambda=sdpvar(1,1,'full');
P=sdpvar(np,np,'symmetric');

flow=[Ac'*P+P*Ac+Qc+tau1*M<=0, P>=eye(np)*1e-8];
jump=[Ad'*P*Ad-P+Qc-tau2*M<=0];
problem=flow+jump+[tau1>=0,tau2>=0, P-lambda*eye(np)<=0];
options=sdpsettings('solver','mosek','verbose',2);
solution=solvesdp(problem,lambda,options);
P=double(P);

%% Simulation of the hybrid system including computation of the cost along the solution
global  Ac Ad M

%initial condition
xp0=[2,6];
% simulation horizon
TSPAN=[0  50];
JSPAN = [0 400];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',0.01);

maxStepCoefficient = .1;  % set the maximum step length. At each run of the
                   % integrator the option 'MaxStep' is set to 
                   % (time length of last integration)*maxStepCoefficient.
                   %  Default value = 0.1

% simulate
x0 = [xp0';0];
[t, x, j] = HyEQsolver(@f,@g,@C,@D,x0,TSPAN,JSPAN,rule,options,maxStepCoefficient);

cost=x(:,3);
%% Plot
figure(1);

hold on;
v2 = [0 0; 0 8; 8 8; 9 9];
hold on;
grid on;

v3 = [0 0; 0 -8; -8 -8; -9 -9];

patch('Vertices',v2,'FaceColor',[0.94 0.94 0.94])

patch('Vertices',v3,'FaceColor',[0.94 0.94 0.94])

grid on;
box on;

plotHarc(x(:,1),j,x(:,2));

figure(2);
plotHarcColor(t,j,cost);





