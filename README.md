# TUW-statistics-calculator
A statistical calculater. Contains functions for every calculation described in the script of course 325.062, Vienna University of Technology.

Explanation (german):
1) Auswählen der benötigten Funktion, entweder über die Menüleiste oder über das DropDown Menü. Die folgenden Einträge erscheinen:
   - Eine kurze Erklärung der Funktion (hier besteht Ausbaubedarf)
   - Eine Option, die Variable, der berechnet werden soll, auszuwählen. Bei vielen Konfidenzintervallen lässt sich z.B. nach "n", Standarabweichung, alpha und dem Radius der Intervalls auflösen, bei den Tests kann man nur "Ausgang" wählen.
2) Auswählen der Variable, nach der aufgelöst werden soll. Es erscheinen die benötigten Eingangswerte sowie Eingabefelder um die Werte einzugeben. Wenn man mit der Maus über den Namen der Eingansgwerte fährt, erscheint eine kurze Beschreibung. (siehe Hinweis 1)
3) Werte einsetzen und "Calculate" klicken.
4) Das Resultat wird im Fenster angezeigt. (siehe Hinweis 2)

Hinweis 1: Bei bestimmten Funktionen geschieht die Eingabe mancher Eingangswerte (z.B. bei der Regressionsanalyse die Stichprobenwerte) nicht über Eingabefelder sondern über Tabellen. Die Tabellen haben zu Beginn nur eine Zeile und Spalte, durch Klicken auf die unterste (rechteste) Spalte (Zeile) werden neue Spalten (Zeilen) hinzugefügt.

Hinweis 2: Bei den aufwändigeren Funktionen wie ANOVA und Regressionsanalyse werden neben dem Resultat auch viele Werte von Zwischenberechnungen in dem Kommandozeilenfenster, das sich mit der Anwendung öffnet, angezeigt (Ich werde schauen, dass diese Werte auch irgendwann im Fenster angezeigt werden). In der Kommandozeile sollten auch ggf. Fehlermeldungen angezeigt werden. Wenn man mich auf diese im Reiter "Issues" auf Github aufmerksam macht, werde ich schauen, dass ich das Problem beseitige.

Am besten man probiert das Programm einfach einmal aus, in dem meisten Fällen sollte aus selbsterklärend sein.

Vorhandenen Funktionen: Alle Konfidenzintervalle aus Kapitel 3, Hypothesentests aus Kapitel 4, einfaktorielle Anova sowie Fisher LSD und Tukey-Test (5) und Regressionsanalyse inklusive Konfidenz- Prädiktionsintervallen und Test für einen Paramter, eine Gruppe von Parametern und globale Signifikanz (6)

Funktionen die noch fehlen: Anpassungstests (Kap 4), mehrfaktorielle Anova, der spezielle Tukey-Test für die Regressionsanalyse aus Kapitel 6.3.2 und die "Modelldiagnose durch Residuenanalyse" (6.4)

Natürlich kann ich nicht für die Richtigkeit der Resultate garantieren, es ist aber alles mit Beispielen aus dem Skriptum und der Aufgabensammlung getestet worden.
