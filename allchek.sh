solver="./lean-proof.sh"
options=""
options=""
time=""
ulimit="ulimit -S -t 30"
# file="qflia-liageneric-nocutlemmas.txt"
# file="qflia-liageneric.txt"

#### Get all files containing "unsat" but not "(xor" nor "(distinct", which we do not and will not support
# grep -L "(xor\|(distinct\|QF_.*\(N\|BV\)\|--incremental" $(grep -l -nr "unsat" test/regress/cli/regress0/)


# file="leanUfRegressions.txt"
file="leanUflraRegressions-smaller.txt"

echo "Options: $options"
echo

while read name; do
    echo "$solver $name";
    $ulimit; $time $solver $name $options;
    echo "=====================================";
done < $file