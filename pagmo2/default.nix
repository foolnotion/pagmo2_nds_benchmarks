let
  pkgs = import ~/nix/nixpkgs {};
in
  pkgs.gcc11Stdenv.mkDerivation {
    name = "pagmo2-env";
    hardeningDisable = [ "all" ]; 
    impureUseNativeOptimizations = true;
    nativeBuildInputs = with pkgs; [ cmake clang_12 clang-tools ];

    buildInputs = with pkgs; [ cmake boost eigen tbb ];
  }
