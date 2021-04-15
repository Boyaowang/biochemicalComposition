CELL = [6 10 5];
XYHW = [5 8 4];
LIGH = [22 28 9];
LIGO = [20 22 10];
LIGC = [15 14 4];
TGL = [57 100 7];
TANN = [15 12 7];

a =[0.3697	0.1622	0.0019	0.2427	0.0374	0.0000	0.1862];
sum(a)

moleTot = a(1)/molarMass(CELL) + a(2)/molarMass(XYHW) ...
        + a(3)/molarMass(LIGH) + a(4)/molarMass(LIGO) ...
        + a(5)/molarMass(LIGC) + a(6)/molarMass(TGL)...
        + a(7)/molarMass(TANN)
    
%% expected mole concentration
YmolCELLExp = a(1)/molarMass(CELL) / moleTot
YmolXYHWExp = a(2)/molarMass(XYHW) / moleTot

SP1Exp = YmolCELLExp + YmolXYHWExp

