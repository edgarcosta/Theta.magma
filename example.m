//SetDebugOnError(true);
//AttachSpec("spec");

SetVerbose("ThetaFlint", 2);
chars := [ Matrix(2*g,1, [c[i] : i in [1..2*g]]) : c in CartesianPower({0,1/2},2*g) ] where g:=2;
R<x> := PolynomialRing(RationalsExtra(500));
C := HyperellipticCurve(R![-2, 2, -4, 2, -1], R![1, 1, 0, 1]);
tau := SmallPeriodMatrix(PeriodMatrix(C));
time tflint := [ThetaFlint(c, ZeroMatrix(Integers(), 2,1), tau) : c in chars];
time t := [Theta(c, ZeroMatrix(Integers(), 2,1), tau) : c in chars];
RealField(8)!Max([Abs(elt) : elt in Eltseq(Vector(tflint) - Vector(t))]);
