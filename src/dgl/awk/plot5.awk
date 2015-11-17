# 17.02.2005 (C) Ingo Eble
#
# plot5.awk - Phasenplot.
#

BEGIN { 

  print NVAR1, NVAR2;

  V1[1] = 2*NVAR1+1; # Column of lower bound of variable 1.
  V1[2] = 2*NVAR1+2; # Column of upper bound of variable 1.
  V2[1] = 2*NVAR2+1; # Column of lower bound of variable 2.
  V2[2] = 2*NVAR2+2; # Column of upper bound of variable 2.
  V3[1] = 2*NVAR3+1; # Column of lower bound of variable 3.
  V3[2] = 2*NVAR3+2; # Column of upper bound of variable 3.
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  print $0 > DF;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 

  BatchString = "set title 'TITEL'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nset zlabel '$y_3$'\nset data style lines\nset grid\n" "set size 1,1\nset border 4096 lw 0.5\n";

  if( V1[2] <= NF && V2[2] <= NF )#&& V3[2] <= NF ) #:($"V3[1]"+0.5*($"V3[2]"-$"V3[1]"))
    {
      BatchString = BatchString "plot '"DF"' using ($"V1[1]"+0.5*($"V1[2]"-$"V1[1]")):($"V2[1]"+0.5*($"V2[2]"-$"V2[1]")) notitle w l lt 1 lw 2\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nreplot\n";
    }
  else
    {
      print "*** Datei enthält nur Daten von den Variablen 1 bis " ((NF-2)/2) " ***";
    }

  BatchString = BatchString "quit\n";

  print BatchString > PF;
} 
