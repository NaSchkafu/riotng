# 17.02.2005 (C) Ingo Eble
#
# plot1.awk - Plot der Loesung als Intervall ohne Zeitbeschriftung.
#

BEGIN { 
  #FS = "[ \t]*(\\[|,|\\]|\\]\\[)[ \t]*"; # Setzt die Trennzeichen: ' [  ,  ]  ][ ';

  print NVAR1, NVAR2;

  V1[1] = 2*NVAR1+1; # Column of lower bound of variable 1.
  V1[2] = 2*NVAR1+2; # Column of upper bound of variable 1.
  V2[1] = 2*NVAR2+1; # Column of lower bound of variable 2.
  V2[2] = 2*NVAR2+2; # Column of upper bound of variable 2.
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
$1 ~ "*" { 

  if( V1[2] > NF || V2[2] > NF ) 
    {
      print "*** Datei enthält nur Daten von den Variablen 1 bis " ((NF-2)/2) " ***";
      exit;
    }

  print $(V1[1]), $(V2[1]) "\n" $(V1[2]), $(V2[1]) "\n" $(V1[2]), $(V2[2]) "\n" $(V1[1]), $(V2[2]) "\n" $(V1[1]), $(V2[1]) "\n" > DF;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 

  BatchString = "set title 'TITEL'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nset data style lines\nset grid\n" "set size 1,1\nset size ratio -1\nset border 4096 lw 0.5\n";

  if( V1[2] <= NF && V2[2] <= NF ) 
    {
      BatchString = BatchString "plot '"DF"' notitle with lines 3\npause -1\nset terminal pslatex color solid norotate\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nset size 1,1\nset size ratio -1\nset border 4096 lw 0.5\nreplot\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nset size 1,1\nset size ratio -1\nset border 4096 lw 0.5\nreplot\n";
    }

  BatchString = BatchString "quit\n";

  print BatchString > PF;
} 
