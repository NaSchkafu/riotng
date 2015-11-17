# 18.02.2005 (C) Ingo Eble
#
# auslesen.awk 
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
	  Inhalt[0] = $1;
	  Inhalt[1] = $2;
	}
      else
	{
	  print Inhalt[0]-1, Inhalt[1];
	  
	  Inhalt[0] = $1;
	  Inhalt[1] = $2;
	  
	  Zeile = $1;
	}
    }
}

#
END { 

  # Letzte Zeile auch noch 
  print Inhalt[0], Inhalt[1];

} 
