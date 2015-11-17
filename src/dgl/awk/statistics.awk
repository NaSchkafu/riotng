# 20.01.2005 (C) Ingo Eble
#
# statistics.awk - Berechnet ein paar statistische Kenngrößen aus den Ausgabedateien.
#

BEGIN { 
  
  #print "Funktioniert leider nicht aufgrund des d Zahlenformats";
  #exit 1;

  CONVFMT = "%15.9E";
  OFMT = "%15.9E";

  max_h = 0;
  min_h = 100;
}

# Die folgende in '{...}' gefasste Aktion wird auf jede Zeile angewendet.
{ 
  if( NR > 1 ) # Erste Zeile übergehen.
    {
      if( FT == 1 ) # File-Typ
	{
	  if( NF != 13 ) 
	    {
	      print "Datei muss 13 Felder enthalten!";
	      exit;
	    }
      
	  if( max_h < $3 ) max_h = $3; # Maximale Schrittweite.
	  if( min_h > $3 ) min_h = $3; # Minimale Schrittweite.
	  m_h = m_h + $3;              # Fuer Mittelwert.
      
	  max_iter = max_iter + $5; # Max. gemachte Zahl an Iterationen

# 	  if( $8 == "no" ) 
# 	    {
# 	      errtol  = errtol  +   1; # Wie oft Fehlertoleranz nicht erreicht?
# 	      errtime = errtime + $10; # Zusätzlich beanspruchte Zeit.
# 	    }
	  if( $8 == "yes" ) decreased = decreased + 1; # Wie oft h verringert weil keine Inklusion erzielt wurde?
	  
	  timestep = timestep + $9;  # Gesamtzeit addieren.
	  timeTM   = timeTM   + $10;  # Zeit fuer TM addieren.
	  timeI0   = timeI0   + $11; # Zeit fuer Inneneinschließung.
	}
      else if( FT == 2 )
	{

	}
      else 
	{
	  print "Falsche Eingabe";
	  exit;
	}
    }
}

# In dieser Funktion werden die Kenngroessen ausgerechnet und ausgegeben.
END { 
  
  if( NR > 1 )
    {
      if( FT == 1 )
	{
	  m_timestep = timestep / (NR-1); # NR-1, da erste Zeile nicht gilt.
	  m_timeTM   = timeTM   / (NR-1);
	  m_timeI0   = timeI0   / (NR-1);
	  m_h        = m_h      / (NR-1);
	}
      else if( FT == 2 )
	{

	}
    }

  printf("Statistische Auswertung der Daten ergab:\n");
  printf("========================================\n");

  if( FT == 1 )
    {
#      printf("Rechenzeit insgesamt : %6.2f sek. (total), %6.2f sek. (im Mittel)\n", timestep, m_timestep);
#      printf("Rechenzeit TM        : %6.2f sek. (total), %6.2f sek. (im Mittel), %2.2f % (bzgl. Gesamtzeit)\n", timeTM, m_timeTM, timeTM/timestep*100);
#      printf("Rechenzeit Innenein. : %6.2f sek. (total), %6.2f sek. (im Mittel), %2.2f % (bzgl. Gesamtzeit)\n", timeI0, m_timeI0, timeI0/timestep*100);
#      printf("Rechenzeit Loes+Verb.: %6.2f sek. (total), %6.2f sek. (im Mittel), %2.2f % (bzgl. Gesamtzeit)\n\n", (timestep-timeTM-timeI0), (m_timestep-m_timeTM-m_timeI0), (timestep-timeTM-timeI0)/timestep*100);
      
#       printf("LocalError > Tol.    : %6i mal\n", errtol);
#       printf("zusätzliche Zeit     : %6.2f sek. (total), %2.2f % (bzgl. Gesamtzeit)\n\n", errtime, errtime/timestep*100);

      printf("h decreased            : %3i mal\n", decreased);
      printf("maximale Schrittweite  : %13.13E\n", max_h);
      printf("minimale Schrittweite  : %13.13E\n", min_h);
      printf("im Mittel              : %13.13E\n", m_h);
      printf("Max. Iterationen       : %6i\n", max_iter);
    }
  else if( FT == 2 )
    {

    }

} 
