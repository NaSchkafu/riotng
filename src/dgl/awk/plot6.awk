# 17.02.2005 (C) Ingo Eble
#
# plot6.awk - Lösungsmengen.
#

BEGIN { 

  print NVAR1, NVAR2;

  print "WARNUNG: *** Funktionert momentan nur bei 2 Variablen, die dann auch noch 'a' und 'b' heißen müssen ***";

  BatchString = "";

  print "set xlabel 'XLABEL'\nset ylabel 'YLABEL'\nset data style lines\nset grid\nset border 4096 lw 0.5\nset lmargin 8\nset rmargin 4\nset parametric\n" > PF;

  n = 1; # Variablenzaehler.
  m = 0; # Ausgabenzaehler.
  s = 0; # Aufrufzaehler (zaehlt die 'Steps').
  t = 1; # Zaehlt wie oft "Polynomial" zwischen zwei "Step" vorkommt.
}

# Auf Zeilen anwenden, deren erstes Element mit "Step" übereinstimmt
$1 ~ /Step/{ 

  if( s == 1 ) w = n-1;

  T = $6; # Sichere Zeit fuer Ausgabe.

  s = s+1;
  n = 1;
}

# Auf Zeilen anwenden, deren erstes Element mit "Polynomial" übereinstimmt
$1 ~ /Polynomial/{ 
  
  if( n == NVAR1 || n == NVAR2 )
    {
      gsub(/\^/,"**");

      if( n == NVAR1 ) gsub(/Polynomial/,"p" s "(a,b)");
      else             gsub(/Polynomial/,"q" s "(a,b)");

      BatchString = BatchString $0 "\n";
      m = m+1;
    }

  if( m > 1 ) 
    {
#      if( s == 4 || s == 8 || s == 12 || s == 16 )
#	{
      BatchString = BatchString "set title 'Lösungsmenge zur Zeit $t="T"$'\nplot [-1:1] p"s"(t,1),q"s"(t,1) notitle w l 3, p"s"(t,-1),q"s"(t,-1) notitle w l 3, p"s"(1,t),q"s"(1,t) notitle w l 3, p"s"(-1,t),q"s"(-1,t) notitle w l 3\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE_"s".tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nreplot\nset terminal x11\nset xlabel 'XLABEL'\nset ylabel 'YLABEL'\n";
#	}

      m = 0;
      v = 1;
    }

  n = n+1;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben.
END { 

  if( v == 0 ) # Keine zwei Ausgaben kamen zustande, weil nur ein Datensatz vorhanden ist.
    {
      print "*** Datei enthält nur Daten von den Variablen 1 bis " n " ***";
    }
  else if( NVAR2 > w )
    {
      print "*** Datei enthält nur Daten von den Variablen 1 bis " w " ***";
    }
  else
    {#set xrange [-2:2]\nset yrange [0:3.5]\nset grid\n
      BatchString = BatchString "set title ''\nset size (3./5.),1\nset size ratio -1 {1.3,1.3}\nset border 4096 lw 0.5\nset lmargin 8\nset rmargin 4\nplot [-1:1]"; #

      for(i=1; i<s; i++)
	{
	  #if( i >= 46 && i <= 49 ) 
          if( i <= 68 || i == 72 || i == 78 || i == 93 || i == 95 ) #i >= 1 && i <= 3 )
	    {
	  if( i > 1 ) BatchString = BatchString ", ";
	  BatchString = BatchString "p"i"(t,1),q"i"(t,1) notitle w l 3, p"i"(t,-1),q"i"(t,-1) notitle w l 3, p"i"(1,t),q"i"(1,t) notitle w l 3, p"i"(-1,t),q"i"(-1,t) notitle w l 3";
	}
	}
      BatchString = BatchString "\npause -1\nset terminal pslatex color solid norotate auxfile\nset output 'TEXFILEDIRTITELFILE_choice.tex'\nset xlabel 'TLABELX'\nset ylabel 'TLABELY'\nreplot\nquit";

      print BatchString > PF;
    }

  print "\nquit\n" > PF;
} 
