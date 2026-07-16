{ pkgs, lib, config, inputs, ... }:

{
  cachix.enable = false;

  packages = [ pkgs.git pkgs.fuse pkgs.pkg-config pkgs.libiconv ];

  languages.rust.enable = true;
}
