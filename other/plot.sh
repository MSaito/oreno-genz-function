gnz=$1
sdim=$2
declare -a gnznames=("Oscillatory" "ProductPeak" "CornerPeak" "Gaussian" \
                    "C0Function" "Discontinuous")
gn=${gnznames[${gnz}-1]}
gnuplot <<EOF
set terminal postscript eps color
set output "0${gnz}-${gn}-s${sdim}.eps"
set datafile separator ","
set title "${gn} (s = ${sdim})"
plot "${gn}SobolRMSE.${sdim}.txt" using 1:3 title "Sobol" with line lw 2 dt 2 lc rgb "red", \
"${gn}SobolLWRMSE.${sdim}.txt" using 1:3 title "Sobol Low WAFOM" with line lw 2 lc rgb "red", \
"${gn}NXRMSE.${sdim}.txt" using 1:3 title "NX" with line lw 2 dt 2 lc rgb "dark-green", \
"${gn}NXLWRMSE.${sdim}.txt" using 1:3 title "NX Low WAFOM" with line lw 2 lc rgb "dark-green"
EOF
