# 26.02.2005 (C) Ingo Eble
#
# separate_column - Separiert eine Spalte aus übergebenem File.
#
# Stimmt nicht ganz: Im Moment wird nur Zeit und Schrittweite 
# aus 'File'_00.txt extrahiert.
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
	  Inhalt[0] = $C1;
	  Inhalt[1] = $C2;
	}
      else
	{
	  print Zeile, Inhalt[0], Inhalt[1];
	  
	  Inhalt[0] = $C1;
	  Inhalt[1] = $C2;
	  
	  Zeile = $1;
	}
    }
}

END { 

  # Letzte Zeile auch noch.

  print Zeile, Inhalt[0], Inhalt[1];

} 
