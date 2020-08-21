{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.fuse pkgs.cargo pkgs.pkg-config pkgs.rustfmt
  ];
}
