# 18.02.2005 (C) Ingo Eble
#
# plot2.awk - Schrittweite über der Zeit
#

BEGIN { 
  Zeile = 1;
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  if( NR > 1 && NF > 3 )
    {
      if( $1 == Zeile )
	{
	  Inhalt[0] = $2;
	  Inhalt[1] = $3;
	}
      else
	{
	  print Inhalt[0], Inhalt[1] > DF;
	  
	  Inhalt[0] = $2;
	  Inhalt[1] = $3;
	  
	  Zeile = $1;
	}
    }
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 

  # Letzte Zeile auch noch 
  print Inhalt[0], Inhalt[1] > DF;

  BatchString = "set title 'TITEL'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nset grid\nset border 4096 lw 0.5\nset data style lines\nset logscale y\nset lmargin 4\nset rmargin 4\nplot '"DF"' notitle with lines 3\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nreplot\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nreplot\nquit\n";

  print BatchString > PF;
} 
