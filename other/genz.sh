pgm=../src/testpack_digitalnet
gnz=$1
sdim=$2
declare -a dnnames=("NX" "Sobol" "NXLW" "SobolLW")
declare -a dnindex=(0 1 3 4)
declare -a gnznames=("Oscillatory" "ProductPeak" "CornerPeak" "Gaussian" \
                    "C0Function" "Discontinuous")

for ((i=0;i<4;i++)) do
#    output="${gnznames[${gnz}-1]}${dnnames[$i]}RMSE.$sdim.txt"
    output="${gnznames[${gnz}-1]}${dnnames[$i]}.$sdim.txt"
#    $pgm -s $sdim -m 10 -M 18 -S 1 -g $gnz -d ${dnindex[$i]} -r100 > $output
    $pgm -s $sdim -m 10 -M 18 -S 1 -g $gnz -d ${dnindex[$i]} > $output
done
