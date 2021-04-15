clc;
clear;

global mC mH mO;

withAshMoist =[0.4546 0.0621 0.4833];
sum(withAshMoist)

CHO = [0.513105051786092 0.0563305854998943 0.430564362714014]; %C H O
alpha = 0.65; beta = 0.8; gamma = 0.8;
delta = 1; epsilon = 0.82;
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
ub = [1, 0.5, 0.5];

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