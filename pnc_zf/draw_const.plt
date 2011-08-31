set style line 1 lt 1 lc rgb "black" lw 1
set style line 2 lt 3 lc rgb "black" lw 1 
set style line 3 lt 2 pt 6 lc rgb "black" lw 1 
set style line 4 lt 2 pt 7 lc rgb "black" lw 1
set style line 5 lt 2 pt 8 lc rgb "black" lw 1  
set style line 6 lt 2 pt 9 lc rgb "black" lw 1  

#set term postscript eps enhanced color blacktext "Helvetica" 12
#set output "ser.eps"
set term png enhanced
set output "const.png"
#set xrange [-2:2]
#set yrange [-2:2]

plot "sp_constellation.txt" using 1:2 with points pt 1 title "label 00", \
	 "sp_constellation.txt" using 3:4 with points pt 2 title "label 01", \
	 "sp_constellation.txt" using 5:6 with points pt 3 title "label 10", \
	 "sp_constellation.txt" using 7:8 with points pt 4 title "label 11"




