# seals-2
Make (or select) a specific directory %libdir% for a local Perl library.

Place the YIW/ directory with its content in %libdir% (e.g. ~/perl5/lib/YIW).

Setup your shell to include PERL_LIB_YIW environment variable with the value of %libdir%.

E.g. with bash:

export PERL_LIB_YIW=~/perl5/lib

Place the scripts in bin/ directories into a location, included in PATH environment variable.

If your system does not have Perl at /usr/bin/perl, edit the path in the shebang line of all scripts to point to the available Perl binary.
