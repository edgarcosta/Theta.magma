

theta_flint_cache := NewStore();


intrinsic CacheClearThetaFlint()
{Clear the internal cache for GeometricHomomorphismRepresentationCC}
  // We need to save and restore the id, otherwise horrific things might
  // happen
  StoreClear(theta_flint_cache);
  StoreSet(theta_flint_cache, "cache", AssociativeArray());
end intrinsic;

intrinsic CacheThetaFlint() -> Assoc
{Clear the internal cache for GeometricHomomorphismRepresentationCC}
  // We need to save and restore the id, otherwise horrific things might
  // happen
  if not StoreIsDefined(theta_flint_cache, "cache") then
    StoreSet(theta_flint_cache, "cache", AssociativeArray());
  end if;
  return StoreGet(theta_flint_cache, "cache");
end intrinsic;

function CacheKeys(tau, z)
  g := Nrows(tau);
  assert Ncols(z) eq 1;
  assert Nrows(z) eq g;
  assert Ncols(tau) eq g;
  try
    precz := Precision(BaseRing(z));
  catch e
    precz := -1;
    z := ChangeRing(z, Rationals());
  end try;
  k1 := <g, Precision(BaseRing(tau)), precz>;
  k2 := <tau, z>;
  return k1, k2;
end function;
procedure SetCache(tau, z, v)
  bool, cache := StoreIsDefined(theta_flint_cache, "cache");
  if not bool then
    cache := AssociativeArray();
  end if;
  k1, k2 := CacheKeys(tau, z);
  if not IsDefined(cache, k1) then
    cache[k1] := AssociativeArray();
  end if;
  cache[k1][k2] := v;
  StoreSet(theta_flint_cache, "cache", cache);
end procedure;

function GetCache(tau, z)
  k1, k2 := CacheKeys(tau, z);
  bool, cache := StoreIsDefined(theta_flint_cache, "cache");
  if not bool then
    cache := AssociativeArray();
    return false, _;
  end if;
  bool, cache1 := IsDefined(cache, k1);
  if not bool then
    return false, _;
  end if;
  return IsDefined(cache1[k2]);
end function;

forward call_python_flint;
intrinsic ThetaFlint(z::Mtrx, tau::Mtrx[FldCom]) -> FldComElt
{ Computes the multidimensional theta function for all characteristics at z (a gx1 matrix) and tau (a symmetric gxg matrix with positive definite imaginary part).}
  bool, v := GetCache(tau, z);
  if bool then
    return v;
  end if;
  v := call_python_flint(tau, z);
  SetCache(tau, z, v);
  return v;
end intrinsic;

to_arb := func<elt | ReplaceCharacter(Sprintf("arb('%o +/- %.*o')", elt, 3, Max(Abs(elt)*10^(1-Precision(Parent(elt))), 10^(1-Precision(Parent(elt)))) ), "E", "e")>;
to_acb := func<elt | Sprintf("acb(%o, %o)", to_arb(Real(elt)), to_arb(Imaginary(elt)))>;
to_acb_list := func<elt | Sprintf("[%o]", Join([to_acb(x) : x in elt], ", "))> ;
to_acb_matrix := func<elt | Sprintf("acb_mat([%o])", Join([to_acb_list(Eltseq(x)) : x in Rows(elt)], ", ")) >;


procedure install_python_flint()
cmd := "
nl = '\\n'
import sys
import subprocess

cmd = [sys.executable, '-m', 'pip', 'install', '--pre', '--upgrade', 'python-flint']

if sys.prefix == sys.base_prefix:
    cmd.append('--user')
cmdfull = ' '.join(cmd)
try:
    from flint import acb
except ImportError:
    _ = sys.stderr.write('Trying to install python-flint' + nl)
    _ = sys.stderr.write(f'Running: {cmdfull}' + nl)
    subprocess.run(cmd, check=True, stdout=sys.stderr.buffer)
";
  Pipe("python", cmd);
end procedure;

function call_python_flint(tau, z)
  cmd := "
try:
    from flint import acb_mat, acb, arb, ctx
except ImportError:
    print('ImportError')
    exit(0)

finally:
    from flint import acb_mat, acb, arb, ctx

log_10_2 = arb.const_log2()/arb.const_log10()
def arb_to_magma(x):
    digits = (max(x.rel_accuracy_bits(), x.rel_one_accuracy_bits())*log_10_2).floor().unique_fmpz()
    digits_str = f'p{digits}'
    return x.str(digits, more=True, radius=False) + digits_str, digits
def acb_to_magma(x):
    real_str, real_digits = arb_to_magma(x.real)
    imag_str, imag_digits = arb_to_magma(x.imag)
    return f'[{real_str}, {imag_str}]', min(real_digits, imag_digits)
def acb_entries_to_magma(m):
    pairs = [acb_to_magma(x) for x in m]
    digits = min(x[1] for x in pairs)
    return f'[ ComplexField({digits}) | ' + ', '.join([x[0] for x in pairs])  + ']'

ctx.dps = %o
tau = %o
z = %o
theta = tau.theta(z)
print(acb_entries_to_magma(theta))
";
  digits := Precision(BaseRing(tau));
  g := Nrows(tau);
  assert g eq Ncols(tau);
  assert g eq Nrows(z);
  assert 1 eq Ncols(z);
  if IsZero(z) then
    acb_z := Sprintf("acb_mat([[0] for _ in range(%o)])", g);
  else
    acb_z := to_acb_matrix(z);
    try
      digits := Max(digits, Precision(BaseRing(z)));
    catch e
      _ := true;
    end try;
  end if;
  working_digits := Ceiling(digits*1.1 + 10);
  acb_tau := to_acb_matrix(tau);
  out := Pipe("python", Sprintf(cmd, working_digits, acb_tau, acb_z));
  if out eq "ImportError\n" then
    install_python_flint();
    out := Pipe("python", Sprintf(cmd, working_digits, acb_tau, acb_z));
    if out eq "ImportError\n" then
      error "Could not install python-flint";
    end if;
  end if;
  return eval out;
end function;



