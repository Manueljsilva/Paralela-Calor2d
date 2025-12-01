{ pkgs ? import <nixpkgs> {} }:

let
  python = pkgs.python312;
  pythonEnv = python.withPackages (ps: with ps; [
    pip
    numpy
    pandas
    matplotlib
    scikit-learn
  ]);
in

pkgs.mkShell {
  name = "develop";

  buildInputs = [
    pkgs.mpich
    pkgs.llvmPackages_13.openmp
    pkgs.gcc
    pkgs.gnumake
    pythonEnv
  ];

  shellHook = ''
    export PYTHONNOUSERSITE=1
    printf "Entered nix-shell: python=%s\n" "${python.interpreter}"
  '';
}
