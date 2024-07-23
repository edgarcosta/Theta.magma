python_flint_version := "0.7.0a4";

declare verbose ThetaFlint, 2;
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
  return IsDefined(cache1, k2);
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


function GetPaths()
  filenames := GetFilenames(CacheClearThetaFlint);
  assert #filenames eq 1;
  package_path := "/" cat Join(s[1..(#s - 1)], "/") where s := Split(filenames[1,1],"/");
  venv_path := package_path cat "/flint";
  python_path := "'" cat venv_path cat "/bin/python'";
  return [venv_path, python_path];
end function;

procedure create_venv(venv_path)
  try
    _ := Pipe(Sprintf("test -d %o", venv_path), "");
  catch e
    vprintf ThetaFlint: "creating venv python-flint...";
    cmd := "python3 -m venv '%o' --without-pip";
    vtime ThetaFlint:
    _ := Pipe("sh", Sprintf(cmd, venv_path));
    _ := Pipe(Sprintf("test -d %o", venv_path), "");
  end try;
end procedure;

function python_version()
  version_info := Pipe("python3 --version", ""); // Python 3.X.Y
  return Join(Split(Split(version_info, " ")[2], ".")[1..2], ".");
end function;

procedure call_pip(venv_path)
  version := python_version();
  sites_path := Sprintf("%o/lib/python%o/site-packages", venv, version);
  package_path := Sprintf("%o/python_flint-%o.dist-info", sites_path, python_flint_version);
  // check if the package is there
  try
    _ := Pipe(Sprintf("test -d %o", package_path), "");
  catch e // if not installs it
    vprintf ThetaFlint: "installing python-flint...";
    vtime ThetaFlint:
      cmd := "python3 -m pip install python-flint==%o --no-input --disable-pip-version-check --force-reinstall --pre --target='%o'";
    _ := Pipe(Sprintf(cmd, python_flint_version, sites_path), "");
    _ := Pipe(Sprintf("test -d %o", package_path), "");
  end try;
end procedure;

procedure install_python_flint()
  venv := GetPaths()[1];
  create_venv(venv);
  call_pip(venv);
end procedure;

function call_python_flint(tau, z)
  cmd := "
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
  cmd := Sprintf(cmd, working_digits, acb_tau, acb_z);
  python_path := GetPaths()[2];
  try
    _ := Pipe(Sprintf("test -f %o", python_path), "");
  catch e
    install_python_flint();
    _ := Pipe(Sprintf("test -f %o", python_path), "");
  end try;

  vprintf ThetaFlint: "Calling python...";
  vtime ThetaFlint:
  out := Pipe(python_path, cmd);
  return eval out;
end function;



