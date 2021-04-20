clc;
clear;

global mC mH mO;

%withAshMoist =[0.4546 0.0621 0.4833];
%sum(withAshMoist)

CHO = [50.38    6.19    43.43]; %C H O
CH = [CHO(1)/(CHO(1)+CHO(2)) CHO(2)/(CHO(1)+CHO(2))]
CHO = CHO./100;
alpha = 0.68; beta = 0.95; gamma = 0.39;
delta = 0.8; epsilon = 0.73;


%% involve OPT
%% OPT1
% alphaOpt1 = [-0.586 2.255 0];
% betaOpt1 = [0.995 -0.012 0.162];
% gammaOpt1 = [1.015 -0.045 -0.005];
% deltaOpt1 = [0.294 0.986 0.002];
% epsilonOpt1 = [0.734 -0.372 0.021];
% 
% alpha = alphaOpt1(1) + alphaOpt1(2) * CHO(1) + alphaOpt1(3) * CHO(2)
% beta = betaOpt1(1) + betaOpt1(2) * CHO(1) + betaOpt1(3) * CHO(2)
% gamma = gammaOpt1(1) + gammaOpt1(2) * CHO(1) + gammaOpt1(3) * CHO(2)
% delta = deltaOpt1(1) + deltaOpt1(2) * CHO(1) + deltaOpt1(3) * CHO(2)
% epsilon = epsilonOpt1(1) + epsilonOpt1(2) * CHO(1) + epsilonOpt1(3) * CHO(2)

%% OPT2
alphaOpt2 = [1.503 -0.037 -13.807];
betaOpt2 = [2.079 -2.160 -0.207];
gammaOpt2 = [12.697 -25.284 13.422];
deltaOpt2 = [-1.75 3.428 13.422];
epsilonOpt2 = [-2.339 1.303 41.335];

alpha = alphaOpt2(1) + alphaOpt2(2) * CHO(1) + alphaOpt2(3) * CHO(2)
beta = betaOpt2(1) + betaOpt2(2) * CHO(1) + betaOpt2(3) * CHO(2)
gamma = gammaOpt2(1) + gammaOpt2(2) * CHO(1) + gammaOpt2(3) * CHO(2)
delta = deltaOpt2(1) + deltaOpt2(2) * CHO(1) + deltaOpt2(3) * CHO(2)
epsilon = epsilonOpt2(1) + epsilonOpt2(2) * CHO(1) + epsilonOpt2(3) * CHO(2)

if alpha > 1
    alpha = 1;
end

if beta > 1
    beta = 1;
end

if gamma >1
    gamma =1;
end

if delta >1
    delta =1;
end

if epsilon >1
    epsilon =1;
end

mC = 12;
mH = 1;
mO = 16;

CELL = [6 10 5];
XYHW = [5 8 4];
LIGH = [22 28 9];
LIGC = [15 14 4];
LIGO = [20 22 10];
TGL = [57 100 7];
TANN = [15 12 7];


SP1 = alpha .* CELL + (1-alpha) .* XYHW;
SP2 = delta * beta .* LIGH + delta * (1-beta) .* LIGC + (1-delta) .* TGL;
SP3 = epsilon * gamma .* LIGO + epsilon * (1-gamma) .* LIGC ...
    + (1-epsilon) .* TANN;
SP = [SP1; SP2; SP3];

CHratio = CHO(1)/CHO(2);
COratio = CHO(1)/CHO(3);

A = zeros(3,3);
for i=1:3
    A(1,i)= mC *SP(i,1) - mH * CHratio * SP(i,2);
    A(2,i)= mC *SP(i,1) - mO * COratio * SP(i,3);
    A(3,i) = 1;
end
b = [0.0 0.0 1.0]';

% x = A\b

%conditions
% A * x <= b

fracs=[1.0 ,1.0, 1];
maxFrac=1;

%Aeq *x = beq
Aeq=[1.0, 1.0, 1.0];
beq=1.0;
%lb <= x <= ub
%lower boundary
lb = [0, 0.0, 0.0];
%upper boundary
ub = [1, 1.0, 1.0];

x = lsqlin(A,b,A,fracs,[],[],lb,ub)

totalMM = 0;
for i = 1:3
    totalMM = totalMM + x(i) * molarMass(SP(i,:));
end

YRM1 =  x(1)*molarMass(SP1)/totalMM;

YCELLmol = x(1) * alpha;
YXYHWmol = x(1) * (1-alpha);
YCELL = YCELLmol*molarMass(CELL)/totalMM
YXYHW = YXYHWmol*molarMass(XYHW)/totalMM

YLIGHmol = x(2) * delta * beta;
YLIGC1mol = x(2)* delta *(1-beta);
YLIGOmol = x(3) * epsilon * gamma;
YLIGC2mol = x(3) * epsilon * (1-gamma);
YTGLmol = x(2) * (1-delta);
YTANNmol = x(3) * (1-epsilon);

YLIGH = YLIGHmol*molarMass(LIGH)/totalMM
YLIGO = YLIGOmol*molarMass(LIGO)/totalMM
YLIGC1 = YLIGC1mol*molarMass(LIGC)/totalMM;
YLIGC2 = YLIGC2mol*molarMass(LIGC)/totalMM;
YLIGC = YLIGC1 + YLIGC2
YTGL = YTGLmol*molarMass(TGL)/totalMM
YTANN = YTANNmol*molarMass(TANN)/totalMM

YCELL + YXYHW + YLIGH + YLIGO + YLIGC +YTGL + YTANN