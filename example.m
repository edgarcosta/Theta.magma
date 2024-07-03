//SetDebugOnError(true);
AttachSpec("spec");

chars := [ KMatrixSpace(RationalField(), 4, 1) |
[0, 0, 0, 0],
[ 0, 0, 0, 1/2 ],
[ 0, 0, 1/2, 0 ],
[ 0, 0, 1/2, 1/2 ],
[ 0, 1/2, 0, 0 ],
[ 0, 1/2, 0, 1/2 ],
[ 0, 1/2, 1/2, 0 ],
[ 0, 1/2, 1/2, 1/2 ],
[ 1/2, 0, 0, 0 ],
[ 1/2, 0, 0, 1/2 ],
[ 1/2, 0, 1/2, 0 ],
[ 1/2, 0, 1/2, 1/2 ],
[ 1/2, 1/2, 0, 0 ],
[ 1/2, 1/2, 0, 1/2 ],
[ 1/2, 1/2, 1/2, 0 ],
[ 1/2, 1/2, 1/2, 1/2 ]
];
SetVerbose("ThetaFlint", 2);
R<x> := PolynomialRing(RationalsExtra(500));
C := HyperellipticCurve(R![-2, 2, -4, 2, -1], R![1, 1, 0, 1]);
tau := SmallPeriodMatrix(PeriodMatrix(C));
time tflint := [ThetaFlint(c, ZeroMatrix(Integers(), 2,1), tau) : c in chars];
time t := [Theta(c, ZeroMatrix(Integers(), 2,1), tau) : c in chars];
RealField(8)!Max([Abs(elt) : elt in Eltseq(Vector(tflint) - Vector(t))]);
