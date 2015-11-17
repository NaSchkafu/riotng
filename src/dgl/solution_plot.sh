#!/bin/bash

while [ -z ${TitelDGL:=""} ]
do
    echo -n "Titel der DGL eingeben: ";
    read TitelDGL;
done

echo -n "1) Plot oder 2) statistische Auswertung? [1] ";
read Wert1;

if [ ${Wert1:=1} -lt 1 ]
then 
   Wert1=1;
elif [ $Wert1 -gt 2 ]
then
   Wert1=2;
fi

echo $Wert1; # Gewählten Wert ausgeben.

if [ $Wert1 -eq 1 ]
then
    while [ ${Wert:=0} -eq 0 ] 
    do
	echo "Was soll geplottet werden?";
	echo "(1) Lösungsmenge als Intervalle ohne Zeitbeschriftung";
	echo "(2) Schrittweite über der Zeit";
	echo "(3) Durchmesser über der Zeit";
	echo "(4) Lösungskurven über der Zeit";
	echo "(5) Phasenplot";
	echo "(6) Lösungsmengen";
	echo "(7) Durchmesser und Schrittweite über der Zeit";
	echo -n "Ihre Wahl: ";
	read Wert;
    done

    echo -n "Titel des Plots eingeben: ";
    read Titel;

    #
    # Setzen des Dateinamens.
    #
    case ${Wert:?"Eingabe fehlgeschlagen"}
    in
	    1) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_03.txt";
	;;
	    2) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_00.txt";
	;;
	[3-5]) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_03.txt";
	;;
	    6) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_04.txt";
	;;
	    7) Datei1=${TitelDGL:?"Eingabe fehlgeschlagen"}"_00.txt";
               Datei2=${TitelDGL:?"Eingabe fehlgeschlagen"}"_03.txt";
	;;
	    *) echo "Falsche Eingabe"
    esac

    #
    # Testen ob die Datei existiert.
    #
    if [ -e ${Datei:=""} -o -e ${Datei1:=""} -a -e ${Datei2:=", "} ]
    then
	#
	# Die benoetigten Pfade und Dateinamen setzen.
	#
	GNUplotDir="gnuplot";
	TexfileDir="tex";
	EpsfileDir="eps";
	TmpfileDir="tmp";

	DataFile1=$GNUplotDir"/"$TmpfileDir"/gnuplot.data.1";
	DataFile2=$GNUplotDir"/"$TmpfileDir"/gnuplot.data.2";
	DataFile3=$GNUplotDir"/"$TmpfileDir"/gnuplot.data.3";
	PrgFile1=$GNUplotDir"/"$TmpfileDir"/gnuplot.1";
	PrgFile2=$GNUplotDir"/"$TmpfileDir"/gnuplot.2";
	#
	# Je nach angegebenem Wert das passende Awk-Skript auf die Datei anwenden.
	#
	case $Wert
	in
		1) echo -n "Welche Variablen [1 2]? "; read NVar1 NVar2; 
		   if [ ${NVar1:=0} -eq 0 ] # Prüfe ob NVar1 definiert und ungleich 0.
		   then 
		    NVar1=1;
		   fi
		   if [ ${NVar2:=0} -eq 0 ] 
		   then 
		    let NVar2=NVar1+1;
		   fi
		   awk -f awk/plot1.awk -v NVAR1=$NVar1 -v NVAR2=$NVar2 -v DF=$DataFile1 -v PF=$PrgFile1 $Datei 
		   XLabel="y_"$NVar1; YLabel="y_"$NVar2;
		   TitelDGL=$TitelDGL"_SolIV_"$NVar1"_"$NVar2;
	    ;;
		2) awk -f awk/plot2.awk -v DF=$DataFile1 -v PF=$PrgFile1 $Datei
		   XLabel="t"; YLabel="h";
                   TitelDGL=$TitelDGL"_h";
	    ;;
		3) awk -f awk/plot3.awk -v DF=$DataFile1 -v PF=$PrgFile1 $Datei
		   XLabel="t"; YLabel="y";
		   TitelDGL=$TitelDGL"_d";
	    ;;
		4) echo -n "Welche Funktion [alle]? "; read NVar;
		   awk -f awk/plot4.awk -v NVAR=${NVar:="alle"} -v DF=$DataFile1 -v PF=$PrgFile1 $Datei
		   XLabel="t"; YLabel="y";
  		   TitelDGL=$TitelDGL"_Sol_"$NVar1;
	    ;;
		5) echo -n "Welche Variablen [1 2]? "; read NVar1 NVar2; 
		   if [ ${NVar1:=0} -eq 0 ] # Prüfe ob NVar1 definiert und ungleich 0.
		   then 
		    NVar1=1;
		   fi
		   if [ ${NVar2:=0} -eq 0 ] 
		   then 
		    let NVar2=NVar1+1;
		   fi #-v NVAR3=3
		   awk -f awk/plot5.awk -v NVAR1=$NVar1 -v NVAR2=$NVar2 -v DF=$DataFile1 -v PF=$PrgFile1 $Datei
		   XLabel="y_"$NVar1; YLabel="y_"$NVar2;
   		   TitelDGL=$TitelDGL"_Phase_"$NVar1"_"$NVar2;
	    ;;
		6) echo -n "Welche Variablen [1 2]? "; read NVar1 NVar2; 
		   if [ ${NVar1:=0} -eq 0 ] # Prüfe ob NVar1 definiert und ungleich 0.
		   then 
		    NVar1=1;
		   fi
		   if [ ${NVar2:=0} -eq 0 ] 
		   then 
		    let NVar2=NVar1+1;
		   fi
		   awk -f awk/plot6.awk -v NVAR1=$NVar1 -v NVAR2=$NVar2 -v PF=$PrgFile1 $Datei
		   XLabel="y_"$NVar1; YLabel="y_"$NVar2;
  		   TitelDGL=$TitelDGL"_SolSet_"$NVar1"_"$NVar2;
	    ;;
	        7) awk -f awk/separate_column.awk -v C1=2 -v C2=3 $Datei1 > $DataFile1
                   awk 'NR==1 || ($1 ~ /-/){print $0}' $Datei2 | nl > $DataFile2
		   join $DataFile1 $DataFile2 > $DataFile3
                   awk -f awk/plot7.awk -v DF=$DataFile1 -v PF=$PrgFile1 $DataFile3
		   XLabel="t";
		   TitelDGL=$TitelDGL"_d_h";
            ;;
	esac
	#
	# Abschließendes bearbeiten der Gnuplot-Batch Datei.
	#
	sed -e's/TITELFILE/'$TitelDGL'/g' $PrgFile1 > $PrgFile2
	sed -e's/TITEL/'$Titel'/g' $PrgFile2 > $PrgFile1
	sed -e's/TEXFILEDIR/'$GNUplotDir'\/'$TexfileDir'\//g' $PrgFile1 > $PrgFile2
	sed -e's/EPSFILEDIR/'$GNUplotDir'\/'$EpsfileDir'\//g' $PrgFile2 > $PrgFile1
	sed -e's/XLABEL/'$XLabel'/g' $PrgFile1 > $PrgFile2
	sed -e's/YLABEL/'$YLabel'/g' $PrgFile2 > $PrgFile1
	sed -e's/TLABELX/$'$XLabel'$/g' $PrgFile1 > $PrgFile2
	sed -e's/TLABELY/$'$YLabel'$/g' $PrgFile2 > $PrgFile1
    	
	#
	# Starten von Gnuplot.
	#
	#$GNUplotDir"/"gnuplot4 $PrgFile1
	gnuplot $PrgFile1
    else
	echo "Error: *** Datei(en)" $Datei $Datei1 $Datei2 " existiert nicht ***";
    fi
else
    while [ ${Wert:=0} -eq 0 ] 
    do
	echo "Welche Datei soll ausgewertet werden?";
	echo "(1) DGL-Alg";
	echo "(2) Shrink-Wrapping 1";
	echo "(3) Shrink-Wrapping 2";
	echo -n "Ihre Wahl: ";
	read Wert;
    done

    #
    # Setzen des Dateinamens.
    #
    case ${Wert:?"Eingabe fehlgeschlagen"}
    in
	1) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_00.txt";
	;;
	2) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_01.txt";
	;;
	3) Datei=${TitelDGL:?"Eingabe fehlgeschlagen"}"_02.txt";
	;;
        *) echo "Falsche Eingabe"
    esac

    #
    # Testen ob die Datei existiert.
    #
    if [ -e $Datei ]
    then
	awk -f awk/statistics.awk -v FT=$Wert $Datei
    else
	echo "Error: *** Datei" $Datei " existiert nicht ***";
    fi
fi
