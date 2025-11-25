{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.mpich
    pkgs.llvmPackages_13.openmp
    pkgs.gcc
    pkgs.gnumake
  ];
}
