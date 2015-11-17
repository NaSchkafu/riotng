# 17.02.2005 (C) Ingo Eble
#
# plot4.awk - Lösungskurven über der Zeit.
#

BEGIN { 
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  print $0 > DF;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 
  
  BatchString = "set title 'TITEL'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nset data style lines\nset grid\nset size 1,1\nset border 4096 lw 0.5\nset lmargin 4\nset rmargin 4\n";

  if( NVAR ~ /alle/ )
    {
      BatchString = BatchString "plot '" DF "'";

      for(i=1;i<=((NF-2)/2);i=i+1)#
	{
	  if( i > 1 ) BatchString = BatchString ", ''";
	  BatchString = BatchString " using 2 : " (2*i+1) " title '$y_"i"$' w l lt "i" lw 2, '' using 2 : " (2*i+2) " notitle w l lt "i" lw 1";
	}

      BatchString = BatchString "\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nreplot\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nreplot\nquit";
    }
  else
    {
      if( NVAR > (NF-2)/2 )
	{
	  print "*** Datei enthält nur Daten der Variablen 1 bis " ((NF-2)/2) " ***";
	  BatchString = BatchString "quit\n";
	}
      else
	{
	  BatchString = BatchString "plot '" DF "' using 2 : " (2*NVAR+1) " title '$y_"NVAR"$' w l lt 1 lw 2, '' using 2 : " (2*NVAR+2) " notitle w l lt 1 lw 1\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE.tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nreplot\nset terminal postscript eps enhanced color solid\nset output 'EPSFILEDIRTITELFILE.eps'\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\nreplot\nquit\n";
	}
    }

  print BatchString > PF;
} 
