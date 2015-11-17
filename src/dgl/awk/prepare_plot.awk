# 15.02.2005 (C) Ingo Eble
#
# prepare_plot.awk - Bereitet die Daten aus 'Title'_00.txt, 'Title'_02.txt und 'Title'_04.txt fuer Gnuplot auf.
#

BEGIN { 

  print "Titel der DGL?";
  getline TitelDGL < "-";

  print "Was soll geplottet werden?";
  print "(1) Lösungsmenge als Intervalle ohne Zeitbeschriftung";
  print "(2) Lösungsmenge als Intervalle mit Zeitbeschriftung";
  print "(3) Schrittweite über der Zeit";
  print "(4) Durchmesser über der Zeit";
  print "(5) Lösungskurven über der Zeit";
  getline What < "-";

  # Wo befindet sich die Information.
  if     ( What == 1 || What == 2 ) File = TitelDGL "_02.txt";
  else if( What == 3              ) File = TitelDGL "_00.txt";
  else if( What == 4 || What == 5 ) File = TitelDGL "_04.txt";

  print "Geben Sie einen Titel für den Plot ein: ";
  getline TitelPlot < "-";
  
  Gnuplot_Dir = "gnuplot/";
  Texfile_Dir = Gnuplot_Dir "tex/";
  Epsfile_Dir = Gnuplot_Dir "eps/";
  Tmpfile_Dir = Gnuplot_Dir "tmp/";

  DataF1 = Tmpfile_Dir "gnuplot.data.1";
  DataF2 = Tmpfile_Dir "gnuplot.data.2";
  PrgF   = Tmpfile_Dir "gnuplot.1";

  n = 0;
  Zeile = 1;

  if( What == 1 || What == 2 )
    FS = "[ \t]*(\\[|,|\\]|\\]\\[)[ \t]*"; # Setzt die Trennzeichen: ' [  ,  ]  ][ ';

  # Lese die entsprechende Datei.
  ARGV[1] = File;
  ARGC    = 2;

#   err = getline $0 < File; 

#   if( err == -1 ) 
#     {
#       print ERRNO;
#       exit;
#     }

#   print $0;
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  if( What == 1 || What == 2 )
    {
      if( NF != 6 ) 
	{
	  print "Datei muss 6 Felder enthalten!";
	  exit;
	}

      print $2, $4 "\n" $3, $4 "\n" $3, $5 "\n" $2, $5 "\n" $2, $4 "\n" > DataF1; 

      if( What == 2 ) # Koordinaten fuer die Beschriftung sammeln.
	{
	  x[n] = $3;
	  y[n] = $4+0.5*($5-$4);
	  t[n] = $6;
	  
	  n = n+1;
	}
    }
  else if( What == 3 )
    {
      if( NR > 1 )
	{
	  if( NF == 13 ) 
	    {
	      if( $1 == Zeile )
		{
		  Inhalt[0] = $2;
		  Inhalt[1] = $3;
		}
	      else
		{
		  print Inhalt[0], Inhalt[1] > DataF2;
		  
		  Inhalt[0] = $2;
		  Inhalt[1] = $3;
		  
		  Zeile = $1;
		}
	    }
	}
    }
  else if( What == 4 ) print $0 > DataF1;
  else if( What == 5 ) print $0 > DataF1;
}

# In dieser Funktion wird der Batch-File fuer Gnuplot geschrieben und Gnuplot
# gestartet. The comments corresponds to the Gnuplot commands.
END { 

  print "set title '" TitelPlot "'" > PrgF;

  if( What == 1 )
    {
      print "set xlabel 'y_1'" > PrgF;
      print "set ylabel 'y_2'" > PrgF;
    }
  else if( What == 2 )
    {
      print "set xlabel 'y_1'" > PrgF;
      print "set ylabel 'y_2'" > PrgF;

      for(i=0;i<n;i=i+1)
	{
	  print "set label", i+1, "'t=" t[i] "' at", x[i] "," y[i] > PrgF;
	}
    }
  else if( What == 3 )
    {
      print "set xlabel 't'" > PrgF;
      print "set ylabel 'h'" > PrgF;
    }
  else if( What == 4 )
    {
      print "set xlabel 't'" > PrgF;
      print "set ylabel 'diam([y])'" > PrgF;
    }
  else if( What == 5 )
    {
      print "set xlabel 't'" > PrgF;
      print "set ylabel 'y'" > PrgF;
    }

  print "set data style lines" > PrgF;

  if( What == 1 )
    {
      print "set grid"             > PrgF;
      print "set size 1,1"   > PrgF;
      print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;

      print "plot '" DataF1 "' notitle with lines 3" > PrgF;
     #print "plot '" DataF1 "' notitle, '" DataF2 "' notitle with lines 3" > PrgF;
    }
  else if( What == 2 )
    {
      print "set grid"             > PrgF;
      print "set size 1,1"   > PrgF;
      print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;

      print "plot '" DataF1 "' notitle with lines 3" > PrgF;
    }
  else if( What == 3 )
    {
      print "set logscale y" > PrgF;
      print "plot '" DataF2 "' notitle with lines 3" > PrgF;
    }
  else if( What == 4 )
    {
      print "set logscale y" > PrgF;
      String = "plot '" DataF1 "' using " NF ":($2-$1) notitle with lines 3";

      for(i=1;i<((NF-1)/2);i=i+1)
	String = String ", '' using " NF ":($"(2*i+2)"-$"(2*i+1)") notitle with lines "(3+i);

      print String > PrgF;
    }
  else if( What == 5 )
    {
      print "set grid"             > PrgF;
      print "set size 1,1"   > PrgF;
      #print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;

      String = "plot '" DataF1 "' using " NF ":1 notitle w l lt 1 lw 2";
      String = String ", '' using " NF ":2 notitle w l lt 1 lw 1";

#       for(i=3;i<NF;i=i+2)
# 	{
# 	  String = String ", '' using " NF ":"i    " notitle w l lt "(1+i)" lw 2";
# 	  String = String ", '' using " NF ":"(i+1)" notitle w l lt "(1+i);
# 	}

      print String > PrgF;
    }

  print "pause -1" > PrgF;                                  # Waits for user input e.g. press any key

  #
  # Save the plot as tex-file.
  #
  
  print "set terminal pslatex color solid norotate" > PrgF; # Changes output characteristics
  print "set output '" Texfile_Dir TitelPlot ".tex'" > PrgF;    # Directs output to external file

  # Change labels for Latex.

  if( What == 1 )
    {
      print "set xlabel '$y_1$'" > PrgF;
      print "set ylabel '$y_2$'" > PrgF;
      print "set size 1,1"   > PrgF;
      print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;
    }
  else if( What == 2 )
    {
      print "set xlabel '$y_1$'" > PrgF;
      print "set ylabel '$y_2$'" > PrgF;
      print "set size 1,1"   > PrgF;
      print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;
      
      for(i=0;i<n;i=i+1)
	{
	  print "set label", i+1, "'$t=" t[i] "$' at", x[i] "," y[i] > PrgF;
	}
    }
  else if( What == 3 )
    {
      print "set xlabel '$t$'" > PrgF;
      print "set ylabel '$h$'" > PrgF;
    }
  else if( What == 4 )
    {
      print "set xlabel '$t$'" > PrgF;
      print "set ylabel 'diam($[y]$)'" > PrgF;
    }
  else if( What == 5 )
    {
      print "set xlabel '$t$'" > PrgF;
      print "set ylabel '$y$'" > PrgF;
    }


  # Re-plot the last plot to file

  print "replot" > PrgF;

  #
  # Save the plot as eps-file.
  #

  print "set terminal postscript eps enhanced color solid" > PrgF; # Changes output characteristics
  print "set output '" Epsfile_Dir TitelPlot ".eps'" > PrgF;           # Directs output to external file

  # Change labels.
  
  if( What == 1 )
    {
      print "set xlabel 'y_1'" > PrgF;
      print "set ylabel 'y_2'" > PrgF;
      print "set size 1,1"   > PrgF;
      print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;
    }
  else if( What == 2 )
    {
      print "set xlabel 'y_1'" > PrgF;
      print "set ylabel 'y_2'" > PrgF;
      print "set size 1,1"   > PrgF;
      print "set size ratio -1"    > PrgF;
      print "set border 4096 lw 0.5" > PrgF;

      for(i=0;i<n;i=i+1)
	{
	  print "set label", i+1, "'t=" t[i] "' at", x[i] "," y[i] > PrgF;
	}
    }
  else if( What == 3 )
    {
      print "set xlabel 't'" > PrgF;
      print "set ylabel 'h'" > PrgF;
    }
  else if( What == 4 )
    {
      print "set xlabel 't'" > PrgF;
      print "set ylabel 'diam([y])'" > PrgF;
    }
  else if( What == 5 )
    {
      print "set xlabel 't'" > PrgF;
      print "set ylabel 'y'" > PrgF;
    }


  # Re-plot the last plot to file

  print "replot" > PrgF;

  #
  # Quit Gnuplot
  #

  print "quit" > PrgF;
  
  #
  # Now start Gnuplot with the above written batch file.
  #

  #system("/home/ae08/gnuplot/bin/gnuplot '" PrgF "'");
  system("gnuplot '" PrgF "'");
} 
