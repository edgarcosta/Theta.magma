function IntegerToHalfCharacteristic(n, g)
  s := IntegerToSequence(n, 2);
  s cat:= [0 : _ in [1..(g - #s)]];
  return Transpose(Matrix(Rationals(), [Reverse(s)]))/2;
end function;

function PairToCharacteristic(pair, g)
  assert #pair eq 2;
  return VerticalJoin([IntegerToHalfCharacteristic(n, g) : n in pair]);
end function;

function HalfCharacteristicToInteger(char)
  assert Ncols(char) eq 1;
  char := ChangeRing(2*char, Integers());
  return SequenceToInteger(Reverse(Eltseq(char)), 2);
end function;

function CharacteristicToPair(char)
  assert Ncols(char) eq 1;
  assert Nrows(char) mod 2 eq 0;
  g := Nrows(char) div 2;
    return <HalfCharacteristicToInteger(Submatrix(char, 1 + g*i, 1, g, 1)) : i in [0,1]>;
end function;

function CharacteristicToInteger(char)
  g := Nrows(char) div 2;
  a, b := Explode(CharacteristicToPair(char));
  return 2^g*a + b;
end function;

function IntegerToCharacteristic(n, g)
  a := n div 2^g;
  b := n mod 2^g;
  return PairToCharacteristic(<a,b>, g);
end function;

function IsPositiveDefiniteImproved(M : prec := 0);
  /* To avoid stupid numerical behavior of IsPositiveDefinite */
  if prec eq 0 then
    RR := BaseRing(M); prec := Precision(RR) div 2;
  end if;
  RRSmall := RealField(prec);
  /* Deal with zero entries that are represented by 10^(-N) */
  for i in [1..#Rows(M)] do
    for j in [1..#Rows(M)] do
      if Abs(M[i,j]) lt 10^(-prec) then
        M[i,j] := 0;
      end if;
    end for;
  end for;
  if Abs(Determinant(M)) lt 10^(-prec) or not IsSymmetricImproved(M) then
    return false;
  end if;
return IsPositiveDefinite(ChangeRing(M, RRSmall));
end function;


intrinsic ThetaFlint(char::Mtrx, z::Mtrx, tau::Mtrx[FldCom]) -> SeqEnum
{ Computes the multidimensional theta function with characteristic char (a (2g)x1 matrix) at z (a gx1 matrix) and tau (a symmetric gxg matrix with positive definite imaginary part).}
  g := NumberOfRows(tau);
  require NumberOfRows(char) eq 2*g:
    "The first argument must have 2g rows.";
  require NumberOfColumns(char) eq 1:
    "The first argument must have 1 column.";
  require NumberOfRows(z) eq g:
    "The second argument must have g rows.";
  require NumberOfColumns(z) eq 1:
    "The second argument must have 1 column.";
  // we force it to be symmetric
  tau := (tau + Transpose(tau))/2;
  Imtau := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
  require IsPositiveDefiniteImproved(Imtau) :
    "The third argument should be a symmetric complex matrix with positive definite imaginary part.";

  n := CharacteristicToInteger(char);
  v := ThetaFlint(z, tau);
  return v[n];
end intrinsic;
