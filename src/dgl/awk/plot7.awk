# 26.02.2005 (C) Ingo Eble
#
# plot7.awk - Schrittweite und Durchmesser über der Zeit.
#
# Zwei Dateien werden angegeben. Die erste wird in BEGIN verarbeitet, die zweite Normal.
#

BEGIN { 
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  print $0 > DF;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 

  BatchString = "set title 'TITEL'\nset xlabel 'XLABEL'\nset grid\nset border 4096 lw 0.5\nset data style lines\nset logscale y\nset lmargin 4\nset rmargin 4\nplot '"DF"' using 5:3 title 'h' with lines -1";

  for(i=1;i<=((NF-5)/2);i=i+1)#
    {
      BatchString = BatchString ", '' using 5:($"(2*i+1+4)"-$"(2*i+4)") title 'd$([y_"i"])$' w l lt "i" lw 2";
    }
  
  BatchString = BatchString "\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nreplot\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nreplot\nquit\n";

  print BatchString > PF;
} 
