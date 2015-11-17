# 17.02.2005 (C) Ingo Eble
#
# plot3.awk - Durchmesser über der Zeit.
#

BEGIN { 
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  print $0 > DF;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 

  BatchString = "set title 'TITEL'\nset xlabel 'XLABEL'\nset grid\nset border 4096 lw 0.5\nset data style lines\nset logscale y\nset lmargin 4\nset rmargin 4\nplot '"DF"'";

  for(i=1;i<=((NF-1)/2);i=i+1)
    {
      if( i > 1 ) BatchString = BatchString ", ''";
      BatchString = BatchString " using 2:($"(2*i+2)"-$"(2*i+1)") title 'diam([$y_"i"$])' w l lt "i" lw 2";
    }
  
  BatchString = BatchString "\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nreplot\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nreplot\nquit\n";

  print BatchString > PF;
} 
