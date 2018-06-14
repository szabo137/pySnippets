here=$(pwd)/sfTrident/phaseInt/pulseLib/cosPulse
gcc -shared -fPIC -o $here/cIntegrands.so $here/cIntegrands.c
echo "compile ctype internal integrals: done"
