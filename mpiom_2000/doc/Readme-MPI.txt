Generelles zum Starten der MPI-Version
=====================================

- Im File OCECTL ist in einer Namelist "NPROCS" die Aufteilung
  des Berechnungsgebiets in x- und y-Richtung anzugeben, d.h.
  am Anfang von OCECTL ist einzufügen:

  &NPROCS
  nprocx=...
  nprocy=...
  &end

- Starten mit

  mpirun -np <numprocs> H70_GROB_run.x

  wobei die Anzahl der Prozessoren natürlich nprocx*nprocy sein muss
  (Sie kann auch nprocx*nprocy+1 sein, in diesem Fall wird im
  "Debug"-Modus gearbeitet, der sollte aber im Normalfall nicht
  verwendet werden. Siehe unten!)

- Der Output erscheint unter

  oceout          für Prozessor 0
  oceout_001      für Prozessor 1
  oceout_002      für Prozessor 2
    ...
  bgcout
  bgcout_001
  bgcout_002
    ...

  siehe hierzu auch Anmerkungen weiter unten!


Vorgehensweise bei der Parallelisierung
=======================================

Verwendete Module
-----------------

Die MPI-Kommunikation ist im Modul mo_mpi.F90 gekapselt,
das ich von Luis Kornblüh übernommen habe (geringfügig
modifiziert), die für die Parallelisierung von MPIOM
wichtigen Variablen und Funktionen befinden sich in
mo_parallel.F90

Generell muss am ganz Anfang des Programms p_start (aus mo_mpi,
zum Starten von MPI) aufgerufen werden, dann müssen die Variablen
nprocx und nprocy (aus mo_parallel) vom Hauptprogramm besetzt werden
(zumindest auf dem I/O pe) und dann muss p_deco aus mo_parallel
gerufen werden, um die Domain-Decomposition zu machen.

- mo_mpi

mo_mpi stellt die folgenden, für den Anwender relevanten Variablen
zu Verfügung:

p_pe     Nummer des eigenen PE
p_io     Nummer des PEs, der I/O macht (i.d.R. 0)
p_nprocs Totale Anzahl von Prozessoren

ausserdem werden die folgenden Funtionen aus mo_mpi direkt
im Programm verwendet:

p_start  Initialisierung
p_abort  Abbruch des Programms (wird von einem PE gerufen)
p_stop   Kontrolliertes Ende (muss von allen PEs gerufen werden)

p_bcast  Broadcast
p_barrier Barrier

mo_mpi stellt darüber hinaus eine Reihe weiterer nützlicher
Funktionen zur Verfügung, die i.d.R. auch in MPIOM benutzt
werden können, lediglich globale Summen (p_sum) sollten nicht
direkt benutzt werden, da das zu Fehlern im Debug-Modus führt.


- mo_parallel

mo_parallel stellt die folgenden Variablen zu Verfügung:

nprocx     Anzahl der Prozessoren in x-Richtung
nprocy     Anzahl der Prozessoren in y-Richtung
nprocxy    Anzahl der Prozessoren, die parallel arbeiten
           nprocxy = nprocx*nprocy
           Im Debug-Modus ist p_nprocs == nprocxy+1, sonst
           ist p_nprocs == nprocxy
p_ioff     Offset zum globalen Gebiet in i-Richtung
p_joff     Offset zum globalen Gebiet in j-Richtung
have_g_is  Flag ob der Prozessor den globalen linken Rand besitzt
have_g_ie  Flag ob der Prozessor den globalen rechten Rand besitzt
have_g_js  Flag ob der Prozessor den globalen unteren Rand besitzt
have_g_je  Flag ob der Prozessor den globalen oberen Rand besitzt

sowie die folgenden Funktionen

p_deco         Initialisierung und Gebietszerlegung
gather_arr     Sammeln eines Globalen Arrays
scatter_arr    Verteilen eines globalen Arrays
bounds_exch    Randaustausch
read_slice     Lesen eines 2d-Feldes uauf dem IO-PE und verteilen an alle
write_slice    Sammeln eines 2d-Feldes und Schreiben durch den IO-PE
global_sum     Globale Summe
global_sum     Globales Maximum
global_min     Globales Minimum
stop_all       Ausgabe einer Fehlermeldung und Aufruf von p_abort


Gebietsgrenzen
--------------

IE und JE sind jetzt Variablen und geben die Gebietsgrenzen auf dem
jeweiligen Prozessor an, der Offset eines Prozessors zum Gebietsanfang
ist in den Variablen p_ioff und p_joff (in mo_parallel) gespeichert.
Die globalen Grenzen des Berechnungsgebeits sind jetzt in IE_G und JE_G,
was nach wie vor Parameter sind (wobei hier nichts dagegen spräche,
auch diese zu Variablen zu machen).

Da IE und JE jetzt Variable sind, werden die meisten Variablen,
die früher in COMMON-Blöcken waren, jetzt dynamisch allokiert.


Globale Felder
--------------

Einige Felder werden auf allen PEs als globale Felder gehalten
(und sind dann doppelt vorhanden), globale Felder haben generell
das Suffix _G um sie von den entsprechenden lokalen Feldern
unterscheiden zu können.


Unterschied globales/lokales i,j
---------------------------------

Beim Rechnen mit globalen Indices ist jeweils die Differenz
p_ioff bzw. p_joff zu den lokalen Indices zu beachten, so schaut
z.B. der Input der Flüsse in OCTHER.F90 jetzt so aus:

      DO N=1,NUMRIV
!
        I=IRIVI(N)-p_ioff
        J=IRIVJ(N)-p_joff
        IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
          ...
          SAO(I,J,1)= ...
          ZO(I,J)=...
          ...
        ENDIF
      ENDDO

oder die Summation der Zeitserien in der BGC

      DO k = 1,4
      DO l=1,nts
         i = its1(l)-p_ioff
         j = jts1(l)-p_joff
         IF(i<=1 .OR. i>=kpie .OR. j<=1 .OR. j>=kpje) CYCLE
         ......
         ts1(1,l+1,lts1) = ts1(1,l+1,lts1) + ...
         ......
      ENDDO
      ENDDO

wobei dann unmittelbar vor der Ausgabe die Globale Summe über alle
Prozessoren gebildet wird.

Zu beachten ist der subtile Unterschied in der Behandlung der Ränder:

Während im ersten Fall die Ränder (i==1, i==ie, j==1, j==je)
mitgenommen werden (sonst müsste man danach einen Randaustausch
machen, falls ein Fluss an der Grenze zweier Berechnungsgebiete liegt)
werden sie im zweiten Fall übersprungen (sonst würden Punkte doppelt
gezählt, falls sie auf der Grenze liegen).


Globale Ränder
--------------

Wenn Randwerte an den J-Rändern gesetzt werden müssen, darf
das nicht mehr mit

      DO I=1,IE
        X(I,1) = ...
        X(I,JE) = ...
      ENDDO

erfolgen, sondern diese Werte dürfen nur von den Prozessoren
gesetzt werden, die auch den globalen Rand besitzen.
Hierfür gibt es in mo_parallel die logischen Variablen
have_g_js und have_g_je (sowie have_g_is und have_g_ie für
die i-Ränder, die allerdings nicht direkt gebraucht werden)
und der Code zum Besetzen der Ränder schaut dann z.B. so aus

       DO I=1,IE
        IF(have_g_js) THEN
          TAUWATU(I,1)=0.
          TAUWATV(I,1)=0.
        ENDIF
        IF(have_g_je) THEN
          TAUWATU(I,JE)=0.
          TAUWATV(I,JE)=0.
        ENDIF
       ENDDO


Randaustausch
-------------

Ein Randaustausch hat überall dort stattzufinden, wo auch früher
die PERIO-Aufrufe waren.

Generell habe ich Wert daruaf gelegt, dass immer das gesamte Gebiet
einschliesslich der Ränder mit gültigen Werten besetzt ist, d.h. es wird
immer dann wenn nur das innere Gebiet berechnet wird, sofort
ein Randaustausch durchgeführt.

Der Randaustausch geschieht jetzt einheitlich mit einem Aufruf
von bounds_exch, sowohl für 2D als auch für 3D Felder, wobei
man bounds_exch bis zu 10 2D-Felder gleichzeitig mitgeben kann
die dann auf einmal ausgetauscht werden (um die Performance zu steigern).

Intern existieren 2 verschiedene Routinen für 2D und 3D-Felder,
die über dasselbe Interface angesprochen werden.


Globale Summen, Maxima, Minima
-----------------------------

Globale Summen werden durch einen Aufruf von global_sum,
global_max, global_min gebildet, die folgende Argumente
akzeptieren:

global_sum:   Bis zu 10 skalare REAL oder INTEGER Variablen
              oder ein 1D-REAL array oder ein 2D-REAL
global_max:   Bis zu 10 skalare REAL Variablen
global_min:   Bis zu 10 skalare REAL Variablen

Zu Beachten ist, dass das Ergebnis von global_sum auf verschiedenen
Prozessorzahlen jeweils leicht anders ausfällt, selbst wenn das
Programm intern mit exakt identischen Werten rechnet.


I/O
---

Der I/O auf stdout wird für jeden Prozessor auf gesonderte Files
ausgegeben  (oceout, oceout_001, ... bgcout, bgcout_001, ....)

Im Prinzip interessieren nur die Files von PE 0 (oceout und bgcout)
es ist aber zu Beachten, dass gewisser Debug-Output jetzt in
verschiedenen Files landen kann, z.B. wenn in BODEN.F90
ausgegeben wird, ob die Topographie verändert wird:

!
!     FLACHE TEILE DER TOPOGRAPHIE KORRIGIEREN
!
      DO J=1,JE
       DO I=1,IE
        IF(DEPTO(I,J).GT.1.AND.DEPTO(I,J).LT.43.5)THEN
         DEPALT=DEPTO(I,J)
         DEPTO(I,J)=MIN(DEPTO(I,J),TOPHILF(I,J))
         DEPTO(I,J)=MAX(DEPTO(I,J),DZW(1)+DZW(2))
         IF(DEPALT.NE.DEPTO(I,J))THEN
          WRITE(IO_STDOUT,*)'TOPOGRAPHIE VERAENDERT: '                  &
     &        ,I+p_ioff,J+p_joff,DEPALT,DEPTO(I,J),TOPHILF(I,J)
         ENDIF
        ENDIF
       ENDDO
      ENDDO

Diese Meldungen landen nun verteilt auf die Files oceout, oceout_001, ...

Ich bin davon ausgegangen, dass diese Meldungen eher zu Debug-Zwecken
dienen und wollte deshalb keinen großen Aufwand hineinstecken,
sie in ein File zu bekommen.

Falls das Programm mit Fehler aussteigt, sollte man allerdings
ALLE Ausgabefiles auf Fehler untersuchen.

Falls es sehr stört, dass solche Meldungen in verschiedenen Files
landen, müssste man sich hier noch eine Lösung überlegen.


Aller anderer I/O wird nur vom I/O-PE gemacht (der in allen mir bekannten
Architekturen der PE 0 sein dürfte), hierzu müssen die Daten
vor dem Schreiben zusammengesammelt bzw. nach dem Lesen verteilt werden.

Das hat zur Folge, dass sich insbesondere die Routinen zum Lesen
und Schreiben der Restart-Files wesentlich verändert haben!!!

Eine Hilfe hierfür sind die Routinen read_slice und write_slice
(in mo_parallel, für normalen unformattierten I/O) bzw.
read_netcdf_var und write_netcdf_var für NETCDF.


Solver
------

Direkte Solver sind nur schwer mit MPI zu parallelisieren,
bei der Gebietsgröße wäre es auch fraglich, ob man hier
durch Parallelisierung schneller wird.
Ich lasse deshalb den direkten Solver komplett auf Prozessor 0
ablaufen und verteile nachher die Ergebnisse.

Iterative Solver sind trivial zu parallelisieren, allerdings
war der iterative Solver so programmiert, dass eine Parallelisierung
durch die Verwendung von verscheidenen Feldern sehr schwierig
war. Ich habe deshalb den iterativen Solver umgeschrieben,
dass er nur noch auf einem Feld arbeitet, was auch der Übersichtlichkeit
sehr zugute kommt. Algorithmus und Ergebnisse sind identisch zur
vorherigen Version.

Falls der direkte Solver aufgrund der Tatsache, dass er nicht parallel
läuft, jemals zur Performance-Bremse werden sollte (was ich aber nicht
glaube), sollte man den iterativen versuchen, wobei allerdings auch hier
die häufigen Randausgleiche (im Vergleich zur Arbeit, die gemacht wird)
wahrscheinlich den Speedup bremsen.



Debug Modus
-----------

Wenn die gesamte Anzahl der Prozessoren auf nprocx*nprocy+1
gesetzt wird, läuft der letzte Prozessor auf dem Gesamtgebiet
und überprüft bei jedem Randaustausch ob das entsprechende
Feld auf den parallel rechnenden Prozessoren und auf dem
Einzelprozessor EXAKT übereinstimmt.
Eine Überprüfung auf Gleicheit kann ausserdem durch einen
Aufruf der Routine para_check an beliebeiger Stelle
gemacht werden.

Das ist sehr hilfreich, um z.B. vergessenen Randausgleichen
auf die Spur zu kommen.

Für den Debug-Modus sollte allerdings die Optimierung ganz
ausgeschaltet werden, da ibs. der NEC-Compiler offensichtlich
Optimierungen trifft, bei denen abhängig von der Länge
eines DO-Loops Ergebnisse mit kleinen numerischen Unterschieden
herauskommen, was dann zu einer Fehlermeldung im Debug-Modus führt.
