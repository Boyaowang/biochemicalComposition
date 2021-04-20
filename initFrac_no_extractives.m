clc;
clear;

global mC mH mO;

withAshMoist =[0.4972    0.0608    0.4420];
sum(withAshMoist)

CHO = [0.534 0.060 0.406]; %C H O
alpha = 0.6; beta = 0.8; gamma = 0.8; %delta = 0.8; epsilon = 1;
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
SP2 = beta .* LIGH + (1-beta) .* LIGC;
SP3 = gamma .* LIGO + (1-gamma) .* LIGC;
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
lb = [0.0, 0.0, 0.0];
%upper boundary
ub = [1.0, 1.0, 1.0];

x = lsqlin(A,b,A,fracs,[],[],lb,ub)

totalMM = 0;
for i = 1:3
    totalMM = totalMM + x(i) * molarMass(SP(i,:));
end

YRM1 =  x(1)*molarMass(SP1)/totalMM

YCELLmol = x(1) * alpha;
YXYHWmol = x(1) * (1-alpha);
YCELL = YCELLmol*molarMass(CELL)/totalMM
YXYHW = YXYHWmol*molarMass(XYHW)/totalMM

YLIGHmol = x(2) * beta;
YLIGC1mol = x(2) * (1-beta);
YLIGOmol = x(3) * gamma;
YLIGC2mol = x(3) * (1-gamma);

YLIGH = YLIGHmol*molarMass(LIGH)/totalMM
YLIGO = YLIGOmol*molarMass(LIGO)/totalMM
YLIGC1 = YLIGC1mol*molarMass(LIGC)/totalMM
YLIGC2 = YLIGC2mol*molarMass(LIGC)/totalMM
YLIGC = YLIGC1 + YLIGC2

%% save the first result for later processing
Ymass1 = [YCELL YXYHW YLIGH YLIGO YLIGC]
for i=1:5
   if Ymass1(i)<1e-5
       Ymass(i)=0;
   end     
end

%% update the results based on new mass fractions
for i=1:5
Ymass2(i) = Ymass1(i)/sum(Ymass1);
end
Ymass2
