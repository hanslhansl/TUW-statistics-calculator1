from scipy.stats import norm, t, chi2, binom, f
import math
import numpy as np
import tkintertable as tkt
import pandas as pd
from statsmodels.sandbox.stats.multicomp import MultiComparison, get_tukeyQcrit2

#Base class
class function:
    name = "Default name"
    description = "Default description"
    def __init__(self):
        self.variableNames = ["x"]
        self.extendedExplanations = {"" : ""}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if variableDict[key] == "":
                raise Exception()
            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception()

        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        pass

    def round(self, value):
        return round(value, 5)

#Konfidenzintervall für Mittelwert, Standardabweichung (σ) der Grundgesamtheit bekannt
class Konfidenzintervall_my_sigma_bekannt(function):
    name = "Konfidenzintervall für Mittelwert, Standardabweichung bekannt"
    description = "Konfidenzintervall für den Mittelwert bei bekannter Standardabweichung σ der Grundgesamtheit."
    def __init__(self):
        self.variableNames = ["n", "σ", "α", "radius"]
        self.necessaryValues = {"n" : ["σ", "α", "radius"], "σ" : ["n", "α", "radius"], "α" : ["n", "σ", "radius"], "radius" : ["n", "σ", "α"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "σ" : "Standardabweichung der Grundgesamtheit", "α" : "Signifikanzniveau des Konfidenzintervalls", "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "n":
            print("looking for n")
            return {looking_for : float((variableDict["σ"] * norm.ppf(variableDict["α"]/2) / variableDict["radius"])**2)}

        elif looking_for == "σ":
            print("looking for σ")
            return {looking_for : float(abs(variableDict["radius"] * variableDict["n"]**0.5 / norm.ppf(variableDict["α"]/2)))}

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(norm.cdf(-variableDict["radius"] * variableDict["n"]**0.5 / variableDict["σ"]) * 2)}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs(variableDict["σ"] * norm.ppf(variableDict["α"]/2) / variableDict["n"]**0.5))}


#Konfidenzintervall für Mittelwert, Standardabweichung (σ) der Grundgesamtheit nicht bekannt
class Konfidenzintervall_my_sigma_unbekannt(function):
    name = "Konfidenzintervall für Mittelwert, Standardabweichung unbekannt"
    description = "Konfidenzintervall für den Mittelwert bei unbekannter Standardabweichung σ der Grundgesamtheit."
    def __init__(self):
        self.variableNames = ["n", "s", "α", "radius"]
        self.necessaryValues = {"s" : ["n", "α", "radius"], "α" : ["n", "s", "radius"], "radius" : ["n", "s", "α"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "s" : "Standardabweichung der Stichprobe", "α" : "Signifikanzniveau des Konfidenzintervalls", "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        if variableDict["n"] != None:
            variableDict["v"] = variableDict["n"] - 1

        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "s":
            print("looking for s")
            return {looking_for : float(abs(variableDict["radius"] * math.sqrt(variableDict["n"]) / t.ppf(variableDict["α"]/2, variableDict["v"])))}

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(2 * t.cdf(-variableDict["radius"] * variableDict["n"]**0.5 / variableDict["s"], variableDict["v"]))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs(variableDict["s"] * t.ppf(variableDict["α"]/2, variableDict["v"]) / math.sqrt(variableDict["n"]) ))}


#Konfidenzintervall für Varianz (σ^2)
class Konfidenzintervall_varianz(function):
    name = "Konfidenzintervall für Varianz"
    description = "Konfidenzintervall für die Varianz einer Grundgesamtheit"
    def __init__(self):
        self.variableNames = ["n", "s²", "α", "radius"]
        self.necessaryValues = {"s²" : ["n", "α", "radius"], "radius" : ["n", "s²", "α"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "s²" : "Varianz der Stichprobe", "α" : "Signifikanzniveau des Konfidenzintervalls", "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        if variableDict["n"] != None:
            variableDict["v"] = variableDict["n"] - 1

        if looking_for == "radius":
            print("looking for radius")
            return {"x1" : variableDict["v"] * variableDict["s²"] / chi2.ppf(1-variableDict["α"]/2, variableDict["v"]), "x2" : variableDict["v"] * variableDict["s²"] / chi2.ppf(variableDict["α"]/2, variableDict["v"])}

        elif looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "s²":
            print("looking for s²")
            return {looking_for : float(2 * variableDict["radius"] / variableDict["v"] / (1/(chi2.ppf(variableDict["α"]/2, variableDict["v"])) - 1/(chi2.ppf(1-variableDict["α"]/2, variableDict["v"])) ) )}

        elif looking_for == "α":
            raise Exception("There is no explicit solution for the variable 'α'!")
            print("looking for α")
            return


#Konfidenzintervall für Anteile (Binomial)
class Konfidenzintervall_Anteile(function):
    name = "Konfidenzintervall für Anteile (Binomial)"
    description = """Konfidenzintervall für die relative Häufigkeit eines Merkmals in einer Grundgesamtheit z.B. für eine Binomialverteilung.
    Gilt für große Stichproben. Faustregel: n*p≥9 und n*q≥9"""
    def __init__(self):
        self.variableNames = ["n", "p", "α", "radius"]
        self.necessaryValues = {"n" : ["p", "α", "radius"], "p" : ["n", "α", "radius"], "α" : ["n", "p", "radius"], "radius" : ["n", "p", "α"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "p" : "relative Häufigkeit des Merkmals in der Stichprobe", "α" : "Signifikanzniveau des Konfidenzintervalls", "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        if variableDict["p"] != None:
            variableDict["q"] = 1 - variableDict["p"]

        if looking_for == "n":
            print("looking for n")
            return {looking_for : float(variableDict["p"] * variableDict["q"] * norm.ppf(variableDict["α"]/2)**2 / variableDict["radius"]**2)}

        elif looking_for == "p":
            print("looking for σ")
            return {"p or q" : float(0.5 - math.sqrt(0.25 - variableDict["radius"]**2 * variableDict["n"] / norm.ppf(variableDict["α"]/2)**2))}

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(2 * norm.cdf( -variableDict["radius"] / math.sqrt(variableDict["p"] * variableDict["q"] / variableDict["n"]) ))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs(norm.ppf(variableDict["α"]/2) * math.sqrt(variableDict["p"] * variableDict["q"] / variableDict["n"])))}


#Prädiktionsintervall für nächste Messung
class Praediktionsintervall(function):
    name = "Prädiktionsintervall für nächste Messung"
    description = "Die nächste Messung fällt mit einer Wahrscheinlichkeit von 1-α in das Prädiktionsintervall."
    def __init__(self):
        self.variableNames = ["n", "s", "α", "radius"]
        self.necessaryValues = {"s" : ["n", "α", "radius"], "α" : ["n", "s", "radius"], "radius" : ["n", "s", "α"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "s" : "Standardabweichung der Stichprobe", "α" : "Signifikanzniveau des Konfidenzintervalls", "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        variableDict["v"] = variableDict["n"] - 1

        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "s":
            print("looking for s")
            return {looking_for : float(abs(variableDict["radius"] / t.ppf(variableDict["α"]/2, variableDict["v"]) / math.sqrt(1 + 1 / variableDict["n"])))}

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(2 * t.cdf(-variableDict["radius"] / variableDict["s"] / math.sqrt(1 + 1 / variableDict["n"]), variableDict["v"]))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs( variableDict["s"] * t.ppf(variableDict["α"]/2, variableDict["v"]) * math.sqrt(1 + 1 / variableDict["n"]) ))}


#Toleranzintervall für zukünftige Messungen
class Toleranzintervall(function):
    name = "Toleranzintervall für zukünftige Messungen"
    description = "k% aller Messungen fallen mit einer Wahrscheinlichkeit von 1-α in das Toleranzintervall.\nDiese Werte sind denen aus dem Skript (S.53) sehr ähnlich (0.1%) aber nicht ident."
    def __init__(self):
        self.variableNames = ["n", "s", "k", "α", "radius"]
        self.necessaryValues = {"s" : ["n", "k", "α", "radius"], "k" : ["n", "s", "α", "radius"], "α" : ["n", "s", "k", "radius"], "radius" : ["n", "s", "k", "α"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "s" : "Standardabweichung der Stichprobe", "k" : "Anteil innerhalb des Intervalls", "α" : "Signifikanzniveau des Konfidenzintervalls", "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        variableDict["v"] = variableDict["n"] - 1

        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "s":
            print("looking for s")
            return {looking_for : float(abs(variableDict["radius"] / norm.ppf((1 + variableDict["k"]) / 2) / math.sqrt( (variableDict["v"] * (1 + 1 / variableDict["n"])) / (chi2.ppf(variableDict["α"], variableDict["v"])))))}

        elif looking_for == "k":
            print("looking for k")
            return {looking_for : float(abs((2 * norm.cdf(variableDict["radius"] / variableDict["s"] / math.sqrt((variableDict["v"] * (1 + 1 / variableDict["n"])) / (chi2.ppf(variableDict["α"], variableDict["v"]))))) - 1))}

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(abs(chi2.cdf(variableDict["v"] * (1 + 1 / variableDict["n"]) * variableDict["s"]**2 * norm.ppf((1 + variableDict["k"]) / 2)**2 / variableDict["radius"]**2, variableDict["v"])))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs(variableDict["s"] *  norm.ppf((1 + variableDict["k"]) / 2) * math.sqrt((variableDict["v"] * (1 + 1 / variableDict["n"])) / (chi2.ppf(variableDict["α"], variableDict["v"])))))}


#Hypothesentest für Anteile (Binomial)
class Hypothesentest_Anteil_binomial_grosser_Umfang(function):
    name = "Hypothesentest für Anteile (Binomial) mit großem Stichprobenumfang"
    description = """Hypothesentest für Anteile (Binomial) mit großem Stichprobenumfang
    H₀: p = p₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α: H₀ wird nicht abgelehnt"""
    def __init__(self):
        self.variableNames = ["n", "p₀", "p\u0302", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n", "p₀", "p\u0302", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "p₀" : "Aussage der Hypothese", "p\u0302" : "relative Häufigkeit des Merkmals in der Stichprobe",
        "α" : "Signifikanzniveau des Tests", "Seite" : "'links', 'rechts' oder 'beide'", "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                pass
            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")

        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")
            zStern = (variableDict["p\u0302"] - variableDict["p₀"]) / (math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]) / variableDict["n"]))

            if variableDict["Seite"] in ("links", "l"):
                PWert = norm.cdf(zStern)
                krit = norm.ppf(variableDict["α"])
                result = zStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = norm.cdf(-zStern)
                krit = norm.ppf(1 - variableDict["α"])
                result = zStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(norm.cdf(zStern), 1-norm.cdf(zStern))
                krit1 = norm.ppf(variableDict["α"]/2) # links
                krit2 = norm.ppf(1 - variableDict["α"]/2) # rechts
                result = (zStern <= krit1) or (zStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Anteile (Binomial): Fehler zweiter Art (β)
class Hypothesentest_Anteil_binomial_grosser_Umfang_Fehler_zweiter_Art(function):
    name = "Hypothesentest für Anteile (Binomial) mit großem Stichprobenumfang: Fehler zweiter Art (β)"
    description = "Wahrscheinlichkeit für einen Fehler zweiter Art (β) bei einem binomialverteilten Hypothesentest mit großem Stichprobenumfang"
    def __init__(self):
        self.variableNames = ["n", "p₀", "p'", "α", "Seite", "β"]
        self.necessaryValues = {"n" : ["p₀", "p'", "α", "Seite", "β"], "β" : ["n", "p₀", "p'", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe", "p₀" : "Aussage der Hypothese", "p'" : "tatsächliche Erfolgswahrscheinlichkeit",
        "α" : "Signifikanzniveau des Tests", "Seite" : "'links', 'rechts' oder 'beide'", "β" : "Wahrscheinlichkeit β"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "β":
            print("looking for β")

            if variableDict["Seite"] in ("links", "l"):
                oben = variableDict["p₀"] - variableDict["p'"] - norm.ppf(variableDict["α"]) * math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]) / variableDict["n"])
                unten = math.sqrt(variableDict["p'"] * (1-variableDict["p'"]) / variableDict["n"])
                return {looking_for : (1 - norm.cdf(oben / unten))}

            elif variableDict["Seite"] in ("rechts", "r"):
                oben = variableDict["p₀"] - variableDict["p'"] + norm.ppf(variableDict["α"]) * math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]) / variableDict["n"])
                unten = math.sqrt(variableDict["p'"] * (1-variableDict["p'"]) / variableDict["n"])
                return {looking_for : norm.cdf(oben / unten)}

            elif variableDict["Seite"] in ("beide", "b"):
                oben1 = variableDict["p₀"] - variableDict["p'"] + norm.ppf(variableDict["α"]/2) * math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]) / variableDict["n"])
                oben2 = variableDict["p₀"] - variableDict["p'"] - norm.ppf(variableDict["α"]/2) * math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]) / variableDict["n"])
                unten = math.sqrt(variableDict["p'"] * (1-variableDict["p'"]) / variableDict["n"])
                return {looking_for : (norm.cdf(oben1 / unten) -  norm.cdf(oben2 / unten))}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")

        if looking_for == "n":
            print("looking for n")

            if variableDict["Seite"] in ("links", "l", "rechts", "r"):
                links = norm.ppf(variableDict["α"]) * math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]))
                rechts = norm.ppf(variableDict["β"]) * math.sqrt(variableDict["p'"] * (1-variableDict["p'"]))
                unten = (variableDict["p'"] - variableDict["p₀"])
                return {looking_for : ((links + rechts)/unten)**2}

            elif variableDict["Seite"] in ("beide", "b"):
                links = norm.ppf(variableDict["α"]/2) * math.sqrt(variableDict["p₀"] * (1-variableDict["p₀"]))
                rechts = norm.ppf(variableDict["β"]) * math.sqrt(variableDict["p'"] * (1-variableDict["p'"]))
                unten = (variableDict["p'"] - variableDict["p₀"])
                return {looking_for : ((links + rechts)/unten)**2}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Anteile (Binomial)
class Hypothesentest_Anteil_binomial_kleiner_Umfang(function):
    name = "Hypothesentest für Anteile (Binomial) mit kleinem Stichprobenumfang (unfertig)"
    description = """Hypothesentest für Anteile (Binomial) mit kleinem Stichprobenumfang
    H₀: p = p₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α: H₀ wird nicht abgelehnt
    Kritischer Bereich ≥Obergrenze und ≤Untergrenze"""
    def __init__(self):
        self.variableNames = []
        self.necessaryValues = {}
        self.extendedExplanations = {}

    def __call__(self, variableDict):
        raise Exception(f"This function is unfinished!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        raise Exception(f"This function is unfinished!")
        if True:
            print("looking for Ausgang")


#Hypothesentest für Anteile (Binomial): Fehler zweiter Art (β)
class Hypothesentest_Anteil_binomial_kleiner_Umfang_Fehler_zweiter_Art(function):
    name = "Hypothesentest für Anteile (Binomial) mit kleinem Stichprobenumfang: Fehler zweiter Art (β)"
    description = """Wahrscheinlichkeit für einen Fehler zweiter Art (β) bei einem Hypothesentest für Anteile (Binomial) mit kleinem Stichprobenumfang
    Kritischer Bereich ≥Obergrenze und ≤Untergrenze
    Power = 1 - β"""
    def __init__(self):
        self.variableNames = ["n", "p'", "Untergrenze", "Obergrenze", "Seite", "β"]
        self.necessaryValues = {"β" : ["n", "p'", "Untergrenze", "Obergrenze", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe",
        "p'" : "Alternativwert für p",
        "Untergrenze" : "untere Grenze",
        "Obergrenze" : "obere Grenze",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "β" : "Wahrscheinlichkeit β"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1 and len(temp) != 2:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one or two should be!")

        if "Ausgang" not in temp:
            raise Exception()
        looking_for = "Ausgang"

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] in ("links", "l"):
                    if variableDict["Obergrenze"] not in ("", None):
                        raise Exception(f"""'Obergrenze' muss für einen linksseitigen Test leer anstatt '{variableDict["Obergrenze"]}' sein!""")
                    else:
                        variableDict["Obergrenze"] = None

                elif variableDict[key] in ("rechts", "r"):
                    if variableDict["Untergrenze"] not in ("", None):
                        raise Exception(f"""'Untergrenze' muss für einen rechtsseitigen Test leer anstatt '{variableDict["Untergrenze"]}' sein!""")
                    else:
                        variableDict["Untergrenze"] = None

                elif variableDict[key] in ("beide", "b"):
                    if variableDict["Obergrenze"] in ("", None):
                        raise Exception("'Obergrenze' darf für einen beidseitigen Test nicht leer sein!")
                    if variableDict["Untergrenze"] in ("", None):
                        raise Exception("'Untergrenze' darf für einen beidseitigen Test nicht leer sein!")
                else:
                    raise Exception(f"The value {variableDict[key]} for 'Seite' is unexpected!")
            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "β":
            print("looking for β")

            if variableDict["Seite"] in ("links", "l"):
                result = 1 - binom.cdf(variableDict["Untergrenze"], variableDict["n"], variableDict["p'"])
                return {looking_for : result, "Power" : 1 - result}

            elif variableDict["Seite"] in ("rechts", "r"):
                result = binom.cdf(variableDict["Obergrenze"]-1, variableDict["n"], variableDict["p'"])
                return {looking_for : result, "Power" : 1 - result}

            elif variableDict["Seite"] in ("beide", "b"):
                result = 1 + binom.cdf(variableDict["Obergrenze"]-1, variableDict["n"], variableDict["p'"]) - binom.cdf(variableDict["Untergrenze"], variableDict["n"], variableDict["p'"])
                return {looking_for : result, "Power" : 1 - result}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")
        return {"Problem:" : "There is no solution to that combination of parameters integrated yet!"}


#Hypothesentest für Mittelwert, Standardabweichung bekannt, eine Stichprobe
class Hypothesentest_my_normal_sigma_bekannt_eine_Probe(function):
    name = "Hypothesentest für Mittelwert, Standardabweichung bekannt, eine Stichprobe"
    description = """Hypothesentest für den Mittelwert einer Grundgesamtheit bei bekannter Standardabweichung σ und einer Stichprobe.
    H₀: p = p₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt"""
    def __init__(self):
        self.variableNames = ["n", "σ", "μ₀", "x\u0304", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n", "σ", "μ₀", "x\u0304", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe",
        "σ" : "Standardabweichung der Grundgesamtheit",
        "μ₀" : "Aussage der Hypothese",
        "x\u0304" : "Mittelwert der Stichprobe",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")
            zStern = (variableDict["x\u0304"] - variableDict["μ₀"]) * math.sqrt(variableDict["n"]) / variableDict["σ"]
            print(zStern)

            if variableDict["Seite"] in ("links", "l"):
                PWert = norm.cdf(zStern)
                krit = norm.ppf(variableDict["α"])
                result = zStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - norm.cdf(zStern)
                krit = norm.ppf(1 - variableDict["α"])
                result = zStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(norm.cdf(zStern), 1-norm.cdf(zStern))
                print(f"""alternativ: {norm.ppf(variableDict["α"]/2)} - {-norm.ppf(variableDict["α"]/2)}""")
                krit1 = norm.ppf(variableDict["α"]/2)
                krit2 = norm.ppf(1 - variableDict["α"]/2)
                result = (zStern <= krit1) or (zStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Mittelwert, Standardabweichung bekannt, eine Stichprobe: Fehler zweiter Art (β)
class Hypothesentest_my_normal_sigma_bekannt_eine_Probe_Fehler_zweiter_Art(function):
    name = "Hypothesentest für Mittelwert, Standardabweichung bekannt, eine Stichprobe: Fehler zweiter Art (β)"
    description = """Wahrscheinlichkeit für einen Fehler zweiter Art (β) bei einem Hypothesentest für den Mittelwert mit bekannter Standardabweichung und einer Stichprobe.
    links: tatsächliche Mittelwert ist kleiner als Aussage der Hypothese
    rechts: tatsächliche Mittelwert ist größer als Aussage der Hypothese
    beide: tatsächliche Mittelwert ist kleiner oder größer als Aussage der Hypothese"""
    def __init__(self):
        self.variableNames = ["n", "σ", "μ₀", "μ'", "α", "Seite", "β"]
        self.necessaryValues = {"n" : ["σ", "μ₀", "μ'", "α", "Seite", "β"], "β" : ["n", "σ", "μ₀", "μ'", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe",
        "σ" : "Standardabweichung der Grundgesamtheit",
        "μ₀" : "Aussage der Hypothese",
        "μ'" : "tatsächlicher Mittelwert",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "β" : "Wahrscheinlichkeit β"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "β":
            print("looking for β")

            if variableDict["Seite"] in ("links", "l"):
                result = 1 - norm.cdf(norm.ppf(variableDict["α"]) + (variableDict["μ₀"] - variableDict["μ'"]) * math.sqrt(variableDict["n"]) / variableDict["σ"])
                return {looking_for : result}

            elif variableDict["Seite"] in ("rechts", "r"):
                result = norm.cdf(-norm.ppf(variableDict["α"]) + (variableDict["μ₀"] - variableDict["μ'"]) * math.sqrt(variableDict["n"]) / variableDict["σ"])
                return {looking_for : result}

            elif variableDict["Seite"] in ("beide", "b"):
                result1 = norm.cdf(-norm.ppf(variableDict["α"]/2) + (variableDict["μ₀"] - variableDict["μ'"]) * math.sqrt(variableDict["n"]) / variableDict["σ"])
                result2 = norm.cdf(norm.ppf(variableDict["α"]/2) + (variableDict["μ₀"] - variableDict["μ'"]) * math.sqrt(variableDict["n"]) / variableDict["σ"])
                return {looking_for : (result1 - result2)}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")

        if looking_for == "n":
            print("looking for n")

            if variableDict["Seite"] in ("links", "l", "rechts", "r"):
                result = ((variableDict["σ"] * (-norm.ppf(variableDict["α"]) - norm.ppf(variableDict["β"]))) / (variableDict["μ₀"] - variableDict["μ'"]))**2
                print(f"""(({variableDict["σ"]} * ({-norm.ppf(variableDict["α"])} - {norm.ppf(variableDict["β"])})) / ({variableDict["μ₀"]} - {variableDict["μ'"]}))**2""")
                return {looking_for : result}

            elif variableDict["Seite"] in ("beide", "b"):
                result = ((variableDict["σ"] * (-norm.ppf(variableDict["α"]/2) - norm.ppf(variableDict["β"]))) / (variableDict["μ₀"] - variableDict["μ'"]))**2
                return {looking_for : result}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Mittelwert, Standardabweichung unbekannt, eine Stichprobe
class Hypothesentest_my_normal_sigma_unbekannt_eine_Probe(function):
    name = "Hypothesentest für Mittelwert, Standardabweichung unbekannt, eine Stichprobe"
    description = """Hypothesentest für den Mittelwert einer Grundgesamtheit bei unbekannter Standardabweichung σ und einer Stichprobe.
    H₀: p = p₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt"""
    def __init__(self):
        self.variableNames = ["n", "s", "x\u0304", "μ₀", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n", "s", "x\u0304", "μ₀", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe",
        "s" : "Standardabweichung der Stichprobe",
        "x\u0304" : "Mittelwert der Stichprobe",
        "μ₀" : "Aussage der Hypothese",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")
            tStern = (variableDict["x\u0304"] - variableDict["μ₀"]) * math.sqrt(variableDict["n"]) / variableDict["s"]
            variableDict["v"] = variableDict["n"] - 1

            if variableDict["Seite"] in ("links", "l"):
                PWert = t.cdf(tStern, variableDict["v"])
                krit = t.ppf(variableDict["α"], variableDict["v"])
                result = tStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - t.cdf(tStern, variableDict["v"])
                krit = t.ppf(1 - variableDict["α"], variableDict["v"])
                result = tStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(t.cdf(tStern, variableDict["v"]), 1-t.cdf(tStern, variableDict["v"]))
                PWert = 2*t.cdf(tStern, variableDict["v"])
                krit1 = t.ppf(variableDict["α"]/2, variableDict["v"]) # links
                krit2 = t.ppf(1 - variableDict["α"]/2, variableDict["v"]) # rechts
                result = (tStern <= krit1) or (tStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Varianz, eine Stichprobe
class Hypothesentest_varianz_normal_eine_Probe(function):
    name = "Hypothesentest für Varianz, eine Stichprobe"
    description = """Hypothesentest für die Varianz σ² einer normalverteilten Grundgesamtheit bei einer Stichprobe.
    H₀: σ² = σ₀²
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt"""
    def __init__(self):
        self.variableNames = ["n", "s²", "σ₀²", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n", "s²", "σ₀²", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe",
        "s²" : "Varianz der Stichprobe",
        "σ₀²" : "Aussage der Hypothese",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")

            variableDict["v"] = variableDict["n"] - 1
            chi2Stern = variableDict["v"] * variableDict["s²"] / variableDict["σ₀²"]
            print("chi2Stern: ", chi2Stern)

            if variableDict["Seite"] in ("links", "l"):
                PWert = chi2.cdf(chi2Stern, variableDict["v"])
                krit = chi2.ppf(variableDict["α"], variableDict["v"])
                result = chi2Stern <= krit
                print("chi2krit: ", chi2.ppf(variableDict["α"], variableDict["v"]))
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : chi2Stern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : chi2Stern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - chi2.cdf(chi2Stern, variableDict["v"])
                krit = chi2.ppf(1-variableDict["α"], variableDict["v"])
                result = chi2Stern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : chi2Stern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : chi2Stern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(chi2.cdf(chi2Stern, variableDict["v"]), 1-chi2.cdf(chi2Stern, variableDict["v"]))
                krit1 = chi2.ppf(variableDict["α"]/2, variableDict["v"]) # links
                krit2 = chi2.ppf(1-variableDict["α"]/2, variableDict["v"]) # rechts
                result = (chi2Stern <= krit1) or (chi2Stern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : chi2Stern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : chi2Stern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Differenz von zwei Mittelwerten, Standardabweichung bekannt, zwei Stichproben
class Hypothesentest_my_Differenz_normal_varianzen_bekannt_zwei_Proben(function):
    name = "Hypothesentest für Differenz zweier Mittelwerte, Varianzen bekannt, zwei Stichproben"
    description = """Hypothesentest für die Differenz zwischen zwei Mittelwerten zweier Grundgesamtheit bei bekannten Varianzen σ₁² und σ₂² und zwei Stichproben.
    H₀: μ₁ - μ₂ = δ₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "x\u0304₁", "x\u0304₂", "σ₁²", "σ₂²", "δ₀", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n₁", "n₂", "x\u0304₁", "x\u0304₂", "σ₁²", "σ₂²", "δ₀", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "x\u0304₁" : "Mittelwert der ersten Stichprobe",
        "x\u0304₂" : "Mittelwert der zweiten Stichprobe",
        "σ₁²" : "Varianz der ersten Grundgesamtheit",
        "σ₂²" : "Varianz der zweiten Grundgesamtheit",
        "δ₀" : "Aussage der Hypothese",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")
            zStern = (variableDict["x\u0304₁"] - variableDict["x\u0304₂"] - variableDict["δ₀"]) / math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"])
            print(zStern)

            if variableDict["Seite"] in ("links", "l"):
                PWert = norm.cdf(zStern)
                krit = norm.ppf(variableDict["α"])
                result = zStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - norm.cdf(zStern)
                krit = -norm.ppf(variableDict["α"])
                result = zStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(norm.cdf(zStern), 1-norm.cdf(zStern))
                krit1 = norm.ppf(variableDict["α"]/2)
                krit2 = -norm.ppf(variableDict["α"]/2)
                result = (zStern <= krit1) or (zStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest für Differenz von zwei Mittelwerten, Standardabweichung bekannt, zwei Stichproben: Fehler zweiter Art (β)
class Hypothesentest_my_Differenz_normal_varianzen_bekannt_zwei_Proben_Fehler_zweiter_Art(function):
    name = "Hypothesentest für Differenz zweier Mittelwerte, Varianzen bekannt, zwei Stichproben: Fehler zweiter Art (β) (muss getestet werden)"
    description = """Wahrscheinlichkeit für einen Fehler zweiter Art (β) bei einem Hypothesentest für die Differenz zwischen zwei Mittelwerten zweier Grundgesamtheiten bei bekannten Varianzen und zwei Stichproben.
    links: tatsächliche Differenz ist kleiner als Aussage der Hypothese
    rechts: tatsächliche Differenz ist größer als Aussage der Hypothese
    beide: tatsächliche Differenz ist kleiner oder größer als Aussage der Hypothese"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "σ₁²", "σ₂²", "δ₀", "δ'", "α", "Seite", "β"]
        self.necessaryValues = {"β" : ["n₁", "n₂", "σ₁²", "σ₂²", "δ₀", "δ'", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "σ₁²" : "Varianz der ersten Grundgesamtheit",
        "σ₂²" : "Varianz der zweiten Grundgesamtheit",
        "δ₀" : "Aussage der Hypothese",
        "δ'" : "tatsächliche Differenz",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "β" : "Wahrscheinlichkeit β"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "β":
            print("looking for β")
            variableDict["σ"] = math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"])

            if variableDict["Seite"] in ("links", "l"):
                result = 1 - norm.cdf(norm.ppf(variableDict["α"]) + (variableDict["δ₀"] - variableDict["δ'"]) / variableDict["σ"])
                return {looking_for : result}

            elif variableDict["Seite"] in ("rechts", "r"):
                result = norm.cdf(-norm.ppf(variableDict["α"]) + (variableDict["δ₀"] - variableDict["δ'"]) / variableDict["σ"])
                return {looking_for : result}

            elif variableDict["Seite"] in ("beide", "b"):
                result1 = norm.cdf(-norm.ppf(variableDict["α"]/2) + (variableDict["δ₀"] - variableDict["δ'"]) / variableDict["σ"])
                result2 = norm.cdf(norm.ppf(variableDict["α"]/2) + (variableDict["δ₀"] - variableDict["δ'"]) / variableDict["σ"])
                return {looking_for : result1 - result2}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Konfidenzintervall für Mittelwert, Varianzen (σ²) der Grundgesamtheiten bekannt, zwei Stichproben
class Konfidenzintervall_my_Differenz_normal_varianzen_bekannt_zwei_Proben(function):
    name = "Konfidenzintervall für Differenz zweier Mittelwerte, Varianzen bekannt, zwei Stichproben"
    description = """Konfidenzintervall für die Differenz zwischen zwei Mittelwerten zweier Grundgesamtheiten bei bekannten Varianzen und zwei Stichproben.
    Das Intervall ist gegeben durch den Radius und den Mittelpunkt x\u0304₁ - x\u0304₂"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "σ₁²", "σ₂²", "α", "radius"]
        self.necessaryValues = {"α" : ["n₁", "n₂", "σ₁²", "σ₂²", "radius"], "radius" : ["n₁", "n₂", "σ₁²", "σ₂²", "α"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "σ₁²" : "Varianz der ersten Grundgesamtheit",
        "σ₂²" : "Varianz der zweiten Grundgesamtheit",
        "α" : "Signifikanzniveau des Konfidenzintervalls",
        "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "σ":
            raise Exception("There is no explicit solution for the variable 'σ'!")
            print("looking for σ")
            return

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(2 * norm.cdf(-variableDict["radius"] / math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"])))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs(norm.ppf(variableDict["α"]/2) * math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"])))}


#Hypothesentest für Differenz von zwei Mittelwerten, Standardabweichung unbekannt, zwei Stichproben
class Hypothesentest_my_Differenz_normal_varianz_unbekannt_zwei_Proben(function):
    name = "Hypothesentest für Differenz zweier Mittelwerte, Varianzen unbekannt, zwei Stichproben (muss getestet werden)"
    description = """Hypothesentest für die Differenz zwischen zwei Mittelwerten zweier Grundgesamtheit bei bekannten Varianzen σ₁² und σ₂² und zwei Stichproben.
    H₀: μ₁ - μ₂ = δ₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt
    Geht für n₁, n₂ > 2000 in den 'Hypothesentest für Differenz zweier Mittelwerte, Varianzen bekannt, zwei Stichproben' über."""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "x\u0304₁", "x\u0304₂", "s₁²", "s₂²", "δ₀", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n₁", "n₂", "x\u0304₁", "x\u0304₂", "s₁²", "s₂²", "δ₀", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "x\u0304₁" : "Mittelwert der ersten Stichprobe",
        "x\u0304₂" : "Mittelwert der zweiten Stichprobe",
        "s₁²" : "Varianz der ersten Stichprobe",
        "s₂²" : "Varianz der zweiten Stichprobe",
        "δ₀" : "Aussage der Hypothese",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        oben = (variableDict["s₁²"] / variableDict["n₁"] + variableDict["s₂²"] / variableDict["n₂"])**2
        unten = ((variableDict["s₁²"] / variableDict["n₁"])**2 / (variableDict["n₁"]-1)) + ((variableDict["s₂²"] / variableDict["n₂"])**2 / (variableDict["n₂"]-1))
        variableDict["v"] = math.floor(oben / unten)

        if looking_for == "Ausgang":
            print("looking for Ausgang")
            tStern = (variableDict["x\u0304₁"] - variableDict["x\u0304₂"] - variableDict["δ₀"]) / math.sqrt(variableDict["s₁²"] / variableDict["n₁"] + variableDict["s₂²"] / variableDict["n₂"])

            if variableDict["Seite"] in ("links", "l"):
                PWert = t.cdf(tStern, variableDict["v"])
                krit = t.ppf(variableDict["α"], variableDict["v"])
                result = tStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - t.cdf(zStern, variableDict["v"])
                krit = -t.ppf(variableDict["α"], variableDict["v"])
                result = zStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(t.cdf(tStern, variableDict["v"]), 1-t.cdf(tStern, variableDict["v"]))
                krit1 = t.ppf(variableDict["α"]/2, variableDict["v"])
                krit2 = -t.ppf(variableDict["α"]/2, variableDict["v"])
                result = (tStern <= krit1) or (tStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Konfidenzintervall für Differenz von zwei Mittelwerten, Varianzen (σ) der Grundgesamtheiten unbekannt, zwei Stichproben
class Konfidenzintervall_my_Differenz_normal_varianzen_unbekannt_zwei_Proben(function):
    name = "Konfidenzintervall für Differenz zweier Mittelwerte, Varianzen unbekannt, zwei Stichproben"
    description = """Konfidenzintervall für die Differenz zwischen zwei Mittelwerten zweier Grundgesamtheiten bei unbekannten Varianzen und zwei Stichproben.
    Das Intervall ist gegeben durch den Radius und den Mittelpunkt x\u0304₁ - x\u0304₂"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "s₁²", "s₂²", "α", "radius"]
        self.necessaryValues = {"α" : ["n₁", "n₂", "s₁²", "s₂²", "radius"], "radius" : ["n₁", "n₂", "s₁²", "s₂²", "α"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "s₁²" : "Varianz der ersten Stichprobe",
        "s₂²" : "Varianz der zweiten Stichprobe",
        "α" : "Signifikanzniveau des Konfidenzintervalls",
        "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        oben = (variableDict["s₁²"] / variableDict["n₁"] + variableDict["s₂²"] / variableDict["n₂"])**2
        unten = ((variableDict["s₁²"] / variableDict["n₁"])**2 / (variableDict["n₁"]-1)) + ((variableDict["s₂²"] / variableDict["n₂"])**2 / (variableDict["n₂"]-1))
        variableDict["v"] = math.floor(oben / unten)

        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "σ":
            raise Exception("There is no explicit solution for the variable 'σ'!")
            print("looking for σ")
            return

        elif looking_for == "α":
            print("looking for α")
            return {looking_for : float(2 * t.cdf(-variableDict["radius"] / math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"]), variableDict["v"]))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : float(abs(t.ppf(variableDict["α"]/2, variableDict["v"]) * math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"])))}


#Hypothesentest für Differenz von zwei Mittelwerten, Varianzen (σ) der Grundgesamtheiten unbekannt aber ident (Pooled Data), zwei Stichproben
class Hypothesentest_my_Differenz_normal_pooled_data_zwei_Proben(function):
    name = "Hypothesentest für Differenz zweier Mittelwerte, Varianzen unbekannt, Pooled Data, zwei Stichproben"
    description = """Hypothesentest für die Differenz zwischen zwei Mittelwerten zweier Grundgesamtheiten bei unbekannten, jedoch identen Varianzen (=Pooled Data) und zwei Stichproben.
    Zuerst wird die geschätzte, gemeinsame Varianz 'sₚ²' errechnet.
    Dann wird basierend auf dieser Varianz ein klassischer t-Test wie auch in 'Hypothesentest für Mittelwert, Standardabweichung unbekannt, eine Stichprobe' durchgeführt."""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "s₁²", "s₂²", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n₁", "n₂", "s₁²", "s₂²", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "s₁²" : "Varianz der ersten Grundgesamtheit",
        "s₂²" : "Varianz der zweiten Grundgesamtheit",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):

        if looking_for == "Ausgang":
            print("looking for Ausgang")
            variableDict["v"] = variableDict["n₁"] + variableDict["n₂"] - 2
            variableDict["sₚ²"] = (variableDict["n₁"]-1) / variableDict["v"] * variableDict["s₁²"] + (variableDict["n₂"]-1) / variableDict["v"] * variableDict["s₂²"]

            if variableDict["Seite"] in ("links", "l"):
                PWert = t.cdf(tStern, variableDict["v"])
                krit = t.ppf(variableDict["α"], variableDict["v"])
                result = tStern <= krit
                if not result:
                    return {"sₚ²" : variableDict["sₚ²"], looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {"sₚ²" : variableDict["sₚ²"], looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - t.cdf(tStern, variableDict["v"])
                krit = t.ppf(1 - variableDict["α"], variableDict["v"])
                result = tStern >= krit
                if not result:
                    return {"sₚ²" : variableDict["sₚ²"], looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {"sₚ²" : variableDict["sₚ²"], looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(t.cdf(tStern, variableDict["v"]), 1-t.cdf(tStern, variableDict["v"]))
                krit1 = t.ppf(variableDict["α"]/2, variableDict["v"]) # links
                krit2 = t.ppf(1 - variableDict["α"]/2, variableDict["v"]) # rechts
                result = (tStern <= krit1) or (tStern >= krit2)
                if not result:
                    return {"sₚ²" : variableDict["sₚ²"], looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {"sₚ²" : variableDict["sₚ²"], looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest Differenz von zwei Anteilen (Binomial), zwei Stichproben
class Hypothesentest_Anteile_Differenz_grosser_Umfang_zwei_Proben(function):
    name = "Hypothesentest für Differenz zweier Anteile (Binomial), großer Umfang, zwei Stichproben"
    description = """Hypothesentest für die Differenz zwischen zwei Anteilen zweier Grundgesamtheiten bei großem Umfang und zwei Stichproben.
    H₀: p\u0302₁ - p\u0302₂ = δ₀"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "p\u0302₁", "p\u0302₂", "δ₀", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n₁", "n₂", "p\u0302₁", "p\u0302₂", "δ₀", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "p\u0302₁" : "relative Häufigkeit des Merkmals in der ersten Stichprobe",
        "p\u0302₂" : "relative Häufigkeit des Merkmals in der zweiten Stichprobe",
        "δ₀" : "Aussage der Hypothese",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")

            variableDict["p\u0302"] = (variableDict["n₁"] / (variableDict["n₁"] + variableDict["n₂"]) * variableDict["p\u0302₁"]
            + variableDict["n₂"] / (variableDict["n₁"] + variableDict["n₂"]) * variableDict["p\u0302₂"])
            zStern = (variableDict["p\u0302₁"] - variableDict["p\u0302₂"]) / math.sqrt(variableDict["p\u0302"] * (1-variableDict["p\u0302"]) * (1/variableDict["n₁"] + 1/variableDict["n₂"]))

            if variableDict["Seite"] in ("links", "l"):
                PWert = norm.cdf(zStern)
                krit = norm.ppf(variableDict["α"])
                result = zStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - norm.cdf(zStern)
                krit = norm.ppf(1 - variableDict["α"])
                result = zStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(norm.cdf(zStern), 1-norm.cdf(zStern))
                krit1 = norm.ppf(variableDict["α"]/2)
                krit2 = norm.ppf(1 - variableDict["α"]/2)
                result = (zStern <= krit1) or (zStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : zStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Hypothesentest Differenz von zwei Anteilen (Binomial), zwei Stichproben: Fehler zweiter Art (β)
class Hypothesentest_Anteile_Differenz_grosser_Umfang_zwei_Proben_Fehler_zweiter_Art(function):
    name = "Hypothesentest für Differenz zweier Anteile (Binomial), großer Umfang, zwei Stichproben: Fehler zweiter Art (β) (muss getestet werden)"
    description = """Wahrscheinlichkeit für einen Fehler zweiter Art (β) bei einem Hypothesentest für die Differenz zwischen zwei Anteilen zweier Grundgesamtheiten bei großem Umfang und zwei Stichproben.
    links: tatsächliche Differenz ist kleiner als Aussage der Hypothese
    rechts: tatsächliche Differenz ist größer als Aussage der Hypothese
    beide: tatsächliche Differenz ist kleiner oder größer als Aussage der Hypothese
    Auflösen für 'n' geht nur für links- bzw. rechtsseitige Tests. Werte stimmen nur, wenn n = n₁ = n₂."""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "p₁'", "p₂'", "α", "Seite", "β"]
        self.necessaryValues = {"n" : ["p₁'", "p₂'", "α", "Seite", "β"], "β" : ["n₁", "n₂", "p₁'", "p₂'", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "p₁'" : "erste tatsächliche Erfolgswahrscheinlichkeit",
        "p₂'" : "zweite tatsächliche Erfolgswahrscheinlichkeit",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "β" : "Wahrscheinlichkeit β"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1 and len(temp) != 2:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one or two should be!")

        if len(temp) == 1:
            if "β" in temp:
                pass
            else:
                raise Exception()
        elif len(temp) == 2:
            if "n₁" in temp and "n₂" in temp:
                pass
            else:
                raise Exception()

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            elif key in ("n₁", "n₂"):
                if variableDict["n₁"] == "" and variableDict["n₂"] == "":
                    variableDict["n₁"] = None
                    variableDict["n₂"] = None
                    looking_for = "n"
                else:
                    raise Exception()
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")

        return self.defaultCall(variableDict, looking_for)


    def defaultCall(self, variableDict, looking_for):
        if looking_for == "β":
            print("looking for β")
            variableDict["p"] = (variableDict["n₁"] * variableDict["p₁'"] + variableDict["n₂"] * variableDict["p₂'"]) / (variableDict["n₁"] + variableDict["n₂"])
            variableDict["q"] = (variableDict["n₁"] * (1-variableDict["p₁'"]) + variableDict["n₂"] * (1-variableDict["p₂'"])) / (variableDict["n₁"] + variableDict["n₂"])
            variableDict["σ"] = math.sqrt(variableDict["p₁'"] * (1-variableDict["p₁'"]) / variableDict["n₁"] + variableDict["p₂'"] * (1-variableDict["p₂'"]) / variableDict["n₂"])
            sqrt = math.sqrt(variableDict["p"] * variableDict["q"] * (1/variableDict["n₁"] + 1/variableDict["n₂"]))

            if variableDict["Seite"] in ("links", "l"):
                result = 1 - norm.cdf((norm.ppf(variableDict["α"]) * sqrt - variableDict["p₁'"] + variableDict["p₂'"])/variableDict["σ"])
                return {looking_for : result}

            elif variableDict["Seite"] in ("rechts", "r"):
                result = norm.cdf((-norm.ppf(variableDict["α"]) * sqrt - variableDict["p₁'"] + variableDict["p₂'"])/variableDict["σ"])
                return {looking_for : result}

            elif variableDict["Seite"] in ("beide", "b"):
                result1 = norm.cdf((-norm.ppf(variableDict["α"]/2) * sqrt - variableDict["p₁'"] + variableDict["p₂'"])/variableDict["σ"])
                result2 = norm.cdf((norm.ppf(variableDict["α"]/2) * sqrt - variableDict["p₁'"] + variableDict["p₂'"])/variableDict["σ"])
                return {looking_for : result1 - result2}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")

        if looking_for == "n":
            if variableDict["Seite"] in ("links", "l") or variableDict["Seite"] in ("rechts", "r"):
                print(f"""α: {norm.ppf(1-variableDict["α"])}, β: {norm.ppf(1-variableDict["β"])}""")
                result = ((norm.ppf(1-variableDict["α"]) * math.sqrt((variableDict["p₁'"] + variableDict["p₂'"]) * ((1-variableDict["p₁'"]) + (1-variableDict["p₂'"])) /2)
                + norm.ppf(1-variableDict["β"]) * math.sqrt(variableDict["p₁'"] * (1-variableDict["p₁'"]) + variableDict["p₂'"] * (1-variableDict["p₂'"])))**2
                / (variableDict["p₁'"] - variableDict["p₂'"])**2)

                return {looking_for : result}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Konfidenzintervall für Differenz von zwei Anteilen (Binomial), zwei Stichproben
class Konfidenzintervall_Anteile_Differenz_grosser_Umfang_zwei_Proben(function):
    name = "Konfidenzintervall für Differenz zweier Anteilen (Binomial), großer Umfang, zwei Stichproben"
    description = """Konfidenzintervall für die Differenz von zwei Anteilen zweier Grundgesamtheiten bei zwei Stichproben.
    Das Intervall ist gegeben durch den Radius und den Mittelpunkt p\u0302₁ - p\u0302₂"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "p\u0302₁", "p\u0302₂", "α", "radius"]
        self.necessaryValues = {"radius" : ["n₁", "n₂", "p\u0302₁", "p\u0302₂", "α"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "p\u0302₁" : "Anteil der ersten Stichprobe",
        "p\u0302₂" : "Anteil der zweiten Stichprobe",
        "α" : "Signifikanzniveau des Konfidenzintervalls",
        "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "σ":
            raise Exception("There is no explicit solution for the variable 'σ'!")
            print("looking for σ")
            return

        elif looking_for == "α":
            raise Exception("There is no explicit solution for the variable 'σ'!")
            print("looking for α")
            return {looking_for : float(2 * t.cdf(-variableDict["radius"] / math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"]), variableDict["v"]))}

        elif looking_for == "radius":
            print("looking for radius")
            return {looking_for : (float(abs(norm.ppf(variableDict["α"]/2) * math.sqrt(variableDict["p\u0302₁"] * (1-variableDict["p\u0302₁"]) / variableDict["n₁"]
            + variableDict["p\u0302₂"] * (1-variableDict["p\u0302₂"]) / variableDict["n₂"]))))}


#Hypothesentest für Gleichheit von zwei Varianzen (σ²) zweier Grundgesamtheiten, zwei Stichproben
class Hypothesentest_varianzen_gleich_zwei_Proben(function):
    name = "Hypothesentest für Gleichheit von zwei Varianzen zweier Grundgesamtheiten, zwei Stichproben"
    description = """Hypothesentest für die Gleichheit von zwei Varianzen zweier Grundgesamtheiten bei zwei Stichproben.
    H₀: σ₁ = σ₂
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt
    Geht für n₁, n₂ > 2000 in den 'Hypothesentest für Differenz zweier Mittelwerte, Varianzen bekannt, zwei Stichproben' über."""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "s₁²", "s₂²", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n₁", "n₂", "s₁²", "s₂²", "α", "Seite"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "s₁²" : "Varianz der ersten Stichprobe",
        "s₂²" : "Varianz der zweiten Stichprobe",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):

        if looking_for == "Ausgang":
            print("looking for Ausgang")
            fStern = variableDict["s₁²"] / variableDict["s₂²"]
            variableDict["v₁"] = variableDict["n₁"] - 1
            variableDict["v₂"] = variableDict["n₂"] - 1

            if variableDict["Seite"] in ("links", "l"):
                PWert = f.cdf(fStern, variableDict["v₁"], variableDict["v₂"])
                krit = f.ppf(variableDict["α"], variableDict["v₁"], variableDict["v₂"])
                result = fStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - f.cdf(fStern, variableDict["v₁"], variableDict["v₂"])
                krit = f.ppf(1-variableDict["α"], variableDict["v₁"], variableDict["v₂"])
                result = fStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(f.cdf(fStern, variableDict["v₁"], variableDict["v₂"]), 1-f.cdf(fStern, variableDict["v₁"], variableDict["v₂"]))
                krit1 = f.ppf(variableDict["α"]/2, variableDict["v₁"], variableDict["v₂"])
                krit2 = f.ppf(1-variableDict["α"]/2, variableDict["v₁"], variableDict["v₂"])
                result = (fStern <= krit1) or (fStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : fStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : fStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#Konfidenzintervall für Quotienten zweier Varianzen, zwei Stichproben
class Konfidenzintervall_my_Quotient_normal_varianzen_unbekannt_zwei_Proben(function):
    name = "Konfidenzintervall für Quotienten zweier Varianzen, zwei Stichproben"
    description = """Konfidenzintervall für den Quotienten zwischen zwei Varianzen zweier Grundgesamtheiten bei zwei Stichproben.
    Das Intervall ist gegeben durch den Radius und den Mittelpunkt s₁² / s₂²"""
    def __init__(self):
        self.variableNames = ["n₁", "n₂", "s₁²", "s₂²", "α", "radius"]
        self.necessaryValues = {"radius" : ["n₁", "n₂", "s₁²", "s₂²", "α"]}
        self.extendedExplanations = {"n₁" : "Umfang der ersten Stichprobe",
        "n₂" : "Umfang der zweiten Stichprobe",
        "s₁²" : "Varianz der ersten Stichprobe",
        "s₂²" : "Varianz der zweiten Stichprobe",
        "α" : "Signifikanzniveau des Konfidenzintervalls",
        "radius" : "Radius des Konfidenzintervalls"}

    def defaultCall(self, variableDict, looking_for):
        variableDict["v₁"] = variableDict["n₁"] - 1
        variableDict["v₂"] = variableDict["n₂"] - 1

        if looking_for == "n":
            raise Exception("There is no explicit solution for the variable 'n'!")
            print("looking for n")
            return

        elif looking_for == "σ":
            raise Exception("There is no explicit solution for the variable 'σ'!")
            print("looking for σ")
            return

        elif looking_for == "α":
            raise Exception("There is no explicit solution for the variable 'σ'!")
            print("looking for α")
            return {looking_for : float(2 * t.cdf(-variableDict["radius"] / math.sqrt(variableDict["σ₁²"] / variableDict["n₁"] + variableDict["σ₂²"] / variableDict["n₂"]), variableDict["v"]))}

        elif looking_for == "radius":
            print("looking for radius")
            x1 = variableDict["s₁²"] / variableDict["s₂²"] * f.ppf(variableDict["α"]/2, variableDict["v"], variableDict["v₁"], variableDict["v₂"])
            x2 = variableDict["s₁²"] / variableDict["s₂²"] * f.ppf(1-variableDict["α"]/2, variableDict["v₁"], variableDict["v₂"])
            return {"x1" : x1, "x2" : x2}


#
class Hypothesentest_gepaarter_t_test(function):
    name = "Hypothesentest für gepaarte Stichproben Mittelwert, zwei Stichproben"
    description = """Hypothesentest für gepaarte Stichproben Mittelwert, zwei Stichproben.
    H₀: μ₁ - μ₂ = μ = δ₀
    P-Wert: Niedrigstes Signifikanzniveau (α), bei dem H₀ gerade noch abgelehnt wird.
    P > α → H₀ wird nicht abgelehnt"""
    def __init__(self):
        self.variableNames = ["n", "sd", "d\u0304", "δ₀", "α", "Seite", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["n", "sd", "d\u0304", "δ₀", "α", "Seite"]}
        self.extendedExplanations = {"n" : "Umfang der Stichprobe",
        "sd" : "Standardabweichung der Stichprobe",
        "d\u0304" : "Mittelwert der Stichprobe",
        "δ₀" : "Aussage der Hypothese",
        "α" : "Signifikanzniveau des Tests",
        "Seite" : "'links', 'rechts' oder 'beide'",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")
            tStern = (variableDict["d\u0304"] - variableDict["δ₀"]) * math.sqrt(variableDict["n"]) / variableDict["sd"]
            variableDict["v"] = variableDict["n"] - 1

            if variableDict["Seite"] in ("links", "l"):
                PWert = t.cdf(tStern, variableDict["v"])
                krit = t.ppf(variableDict["α"], variableDict["v"])
                result = tStern <= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("rechts", "r"):
                PWert = 1 - t.cdf(tStern, variableDict["v"])
                krit = t.ppf(1 - variableDict["α"], variableDict["v"])
                result = tStern >= krit
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritischer Wert" : krit, "P-Wert" : PWert}

            elif variableDict["Seite"] in ("beide", "b"):
                PWert = 2*min(t.cdf(tStern, variableDict["v"]), 1-t.cdf(tStern, variableDict["v"]))
                krit1 = t.ppf(variableDict["α"]/2, variableDict["v"]) # links
                krit2 = t.ppf(1 - variableDict["α"]/2, variableDict["v"]) # rechts
                result = (tStern <= krit1) or (tStern >= krit2)
                if not result:
                    return {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}
                elif result:
                    return {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : tStern, "Kritische Werte" : str(krit1) + "; " + str(krit2), "P-Wert" : PWert}

            else:
                raise Exception(f"The value {variableDict['Seite']} for 'Seite' is unexpected!")


#
class ANOVA_einfache_varianzanalyse(function):
    name = "ANOVA - Einfache Varianzanalyse"
    description = """B.S. 92"""
    def __init__(self):
        self.variableNames = []
        self.necessaryValues = {"Ausgang über Quadratsummen" : ["MSTr", "MSE", "N", "M", "α"], "Ausgang über Werte" : ["M", "α", "Tabelle"]}
        self.extendedExplanations = {"MSTr" : "Mittlere quadratische Abweichung zw. den Proben",
        "MSE" : "Mittlerer quadratischer Fehler",
        "N" : "Anzahl der Proben",
        "M" : "Anzahl der Werte pro Probe",
        "α" : "Signifikanzniveau des Tests",
        "Ausgang über Quadratsummen" : "'true' für nicht abgelehnt, 'false' für abgelehnt",
        "Tabelle" : ["Mittelwerte X\u0304ᵢ", "Varianzen Sᵢ²"],
        "Ausgang über Werte" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        if set(variableDict) == set(self.necessaryValues["Ausgang über Quadratsummen"]):
            looking_for = "Ausgang über Quadratsummen"
            foo = self.necessaryValues["Ausgang über Quadratsummen"]

        elif set(variableDict) == set(self.necessaryValues["Ausgang über Werte"]):
            looking_for = "Ausgang über Werte"
            foo = self.necessaryValues["Ausgang über Werte"]
        else:
            raise Exception()

        temp = set(foo) - set(variableDict)
        if len(temp) != 0:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables ({temp}) are empty whilst exactly zero should be!")

        for key in variableDict:
            if variableDict[key] == "":
                raise Exception(f"{key} must not be empty but is '{variableDict[key]}'")

            if key == "Tabelle":
                tableModel = variableDict[key]
                if type(tableModel) != tkt.TableModel:
                    raise Exception(f"{key} must not be of type 'tkt.TableModel' but is {type(tableModel)}")

                while all([tableModel.getCellRecord(tableModel.getRowCount()-1, k) == None for k in range(tableModel.getColumnCount())]):
                    tableModel.deleteRow(tableModel.getRowCount()-1)
                if tableModel.getRowCount() in (0, 1):
                    raise Exception(f"{key} is empty (={tableModel})")

                self.tableList = []
                for i in range(tableModel.getColumnCount()):
                    column = tableModel.getColCells(i)
                    if column == [""]:
                        raise Exception(f"Column {i} is empty. It needs to have at least one value!")

                    self.tableList.append([])
                    for j in range(len(column)):
                        cell = column[j]

                        try:
                            cell = float(cell.replace(",", "."))
                        except ValueError:
                            raise Exception(f"The value {cell} for {i}x{j} is unexpected!")
                        self.tableList[i].append(cell)
                variableDict[key] = self.tableList

            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")

        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang über Quadratsummen":
            print("looking for Ausgang über Quadratsummen")
            additionalResult = {}

        elif looking_for == "Ausgang über Werte":
            print("looking for Ausgang über Werte")
            variableDict["Xᵢ"] = variableDict["Tabelle"][0]
            variableDict["Sᵢ²"] = variableDict["Tabelle"][1]

            variableDict["N"] = len(variableDict["Xᵢ"])
            variableDict["X\u0304"] = sum(variableDict["Xᵢ"]) / variableDict["N"]
            foo = [(xi - variableDict["X\u0304"])**2 for xi in variableDict["Xᵢ"]]
            variableDict["MSTr"] = sum(foo) * variableDict["M"] / (variableDict["N"] - 1)
            variableDict["MSE"] = sum(variableDict["Sᵢ²"]) / variableDict["N"]

            additionalResult = {"X\u0304" : variableDict["X\u0304"]}

        if looking_for in ("Ausgang über Quadratsummen", "Ausgang über Werte"):
            fStern = variableDict["MSTr"] / variableDict["MSE"]
            SSTr = variableDict["MSTr"] * (variableDict["N"] - 1)
            SSE = variableDict["MSE"] * variableDict["N"] * (variableDict["M"] - 1)
            variableDict["v₁"] = variableDict["N"] - 1
            variableDict["v₂"] = variableDict["N"] * (variableDict["M"] - 1)

            krit = f.ppf(1-variableDict["α"], variableDict["v₁"], variableDict["v₂"])
            result = fStern >= krit

            data = {
            "1" : {"Steuung" : "Treatments", "Freiheitsgrade" : f"N-1={variableDict['v₁']}", "Quadratsumme" : f"SSTr={SSTr}", "Mittlere Quadratsumme" : f"MSTr={variableDict['MSTr']}", "f-Wert" : f"MSTr/MSE={fStern}"},
            "2" : {"Steuung" : "Fehler", "Freiheitsgrade" : f"N(M-1)={variableDict['v₂']}", "Quadratsumme" : f"SSE={SSE}", "Mittlere Quadratsumme" : f"MSE={variableDict['MSE']}", "f-Wert" : f""},
            "3" : {"Steuung" : "Total", "Freiheitsgrade" : f"NM-1={variableDict['N']*variableDict['M']-1}", "Quadratsumme" : f"SST={SSTr+SSE}", "Mittlere Quadratsumme" : f"", "f-Wert" : f""}
            }

            if not result:
                foo =  {looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit,
                "MSTr" : variableDict["MSTr"], "MSE" : variableDict["MSE"], "SSTr" : SSTr, "SSE" : SSE, "data" : data}
            elif result:
                foo = {looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit,
                "MSTr" : variableDict["MSTr"], "MSE" : variableDict["MSE"], "SSTr" : SSTr, "SSE" : SSE, "data" : data}

            return {**foo, **additionalResult}


#
class ANOVA_Fisher_LSD(function):
    name = "ANOVA - Fisher LSD Test"
    description = """Fisher least significant difference Test zwischen zwei Werten
    L ≤ N(N-1)/2 bzw. L = N(N-1)/2
    Für α ist immer der globale Wert einzusetzen, niemals der lokale, der sich durch 1-(1-N)^(1/L) ergibt.
    B.S. 96"""
    def __init__(self):
        self.variableNames = []
        self.necessaryValues = {"Ausgang ohne Bonferroni-Korrektur" : ["X\u0304₁", "X\u0304₂", "N", "M", "MSE", "α"], "Ausgang mit Bonferroni-Korrektur" : ["X\u0304₁", "X\u0304₂", "N", "M", "MSE", "α", "L"]}
        self.extendedExplanations = {"X\u0304₁" : "Durchschnitt Probe 1",
        "X\u0304₂" : "Durchschnitt Probe 2",
        "N" : "Anzahl der Proben",
        "M" : "Anzahl der Werte pro Probe",
        "MSE" : "Mittlerer quadratischer Fehler",
        "α" : "Signifikanzniveau des Tests",
        "Ausgang ohne Bonferroni-Korrektur" : "'true' für nicht abgelehnt, 'false' für abgelehnt",
        "L" : "Anzahl der Vergleiche",
        "Ausgang mit Bonferroni-Korrektur" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        if set(variableDict) == set(self.necessaryValues["Ausgang ohne Bonferroni-Korrektur"]):
            looking_for = "Ausgang ohne Bonferroni-Korrektur"
            foo = self.necessaryValues["Ausgang ohne Bonferroni-Korrektur"]

        elif set(variableDict) == set(self.necessaryValues["Ausgang mit Bonferroni-Korrektur"]):
            looking_for = "Ausgang mit Bonferroni-Korrektur"
            foo = self.necessaryValues["Ausgang mit Bonferroni-Korrektur"]
        else:
            raise Exception()

        temp = set(foo) - set(variableDict)
        if len(temp) != 0:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables ({temp}) are empty whilst exactly zero should be!")

        for key in variableDict:
            if variableDict[key] == "":
                raise Exception(f"{key} must not be empty but is '{variableDict[key]}'")

            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")

        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang ohne Bonferroni-Korrektur":
            pass

        elif looking_for == "Ausgang mit Bonferroni-Korrektur":
            variableDict["α"] = 1 - (1 - variableDict["α"]) ** (1/variableDict["L"])

        if looking_for == "Ausgang ohne Bonferroni-Korrektur" or looking_for == "Ausgang mit Bonferroni-Korrektur":
            variableDict["v"] = variableDict["N"] * (variableDict["M"] - 1)
            tStern = t.ppf(1-variableDict["α"]/2, variableDict["v"])
            q = tStern * math.sqrt(2 * variableDict["MSE"] / variableDict["M"])

            result = abs(variableDict["X\u0304₁"] - variableDict["X\u0304₂"]) > q

            if not result:
                return {looking_for : "True, X\u0304₁ und X\u0304₂ unterscheiden nicht sich signifikant", "Teststatistik" : tStern, "Radius q" : str(q)}
            elif result:
                return {looking_for : "False, X\u0304₁ und X\u0304₂ unterscheiden sich signifikant", "Teststatistik" : tStern, "Radius q" : str(q)}

        else:
            raise Exception()


#
class ANOVA_Tukeys_t_test(function):
    name = "ANOVA - Tukey's T-Test"
    description = """Tukey's T-Test zwischen zwei Werten
    Bei gleichen Stichprobenumfängen im Experiment gilt n = N*M.

    B.S. 96"""
    def __init__(self):
        self.variableNames = ["X\u0304₁", "X\u0304₂", "M₁", "M₂", "MSE", "N", "n", "α", "Ausgang"]
        self.necessaryValues = {"Ausgang" : ["X\u0304₁", "X\u0304₂", "M₁", "M₂", "N", "n", "MSE", "α"]}
        self.extendedExplanations = {"X\u0304₁" : "Durchschnitt Probe 1",
        "X\u0304₂" : "Durchschnitt Probe 2",
        "M₁" : "Anzahl der Werte Probe 1",
        "M₂" : "Anzahl der Werte Probe 2",
        "N" : "Anzahl der Proben",
        "n" : "Gesamtzahl aller Messungen",
        "MSE" : "Mittlerer quadratischer Fehler",
        "α" : "Signifikanzniveau des Tests",
        "Ausgang" : "'true' für nicht abgelehnt, 'false' für abgelehnt"}

    def __call__(self, variableDict):
        temp = set(self.variableNames) - set(variableDict)
        if len(temp) != 1:
            print("variableDict: ", variableDict)
            raise Exception(f"{len(temp)} variables are empty whilst exactly one should be!")
        looking_for = list(temp)[0]

        for key in variableDict:
            if key == "Seite":
                if variableDict[key] == "":
                    raise Exception()
                else:
                    pass
            else:
                if variableDict[key] == "":
                    raise Exception()
                else:
                    try:
                        variableDict[key] = float(variableDict[key].replace(",", "."))
                    except ValueError:
                        raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")
        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            print("looking for Ausgang")

            variableDict["v₁"] = variableDict["N"]
            variableDict["v₂"] = variableDict["n"] - variableDict["N"]

            Q_crit = get_tukeyQcrit2(alpha=variableDict["α"], k=variableDict["v₁"], df=variableDict["v₂"])

            q = Q_crit * math.sqrt(variableDict["MSE"] / 2 * (1 / variableDict["M₁"] + 1 / variableDict["M₂"]))

            result = abs(variableDict["X\u0304₁"] - variableDict["X\u0304₂"]) < q

            if result:
                return {looking_for : "True, X\u0304₁ und X\u0304₂ unterscheiden nicht sich signifikant", "Q-Wert" : Q_crit, "Radius q" : str(q)}
            elif not result:
                return {looking_for : "False, X\u0304₁ und X\u0304₂ unterscheiden sich signifikant", "Q-Wert" : Q_crit, "Radius q" : str(q)}

        else:
            raise Exception()


#
class ANOVA_einfache_varianzanalyse_vollstaendige_Daten(function):
    name = "ANOVA - Einfache Varianzanalyse mit vollständiger Dateneingabe"
    description = """Einfache Varianzanalyse mit Eingabe aller gemessener Werte
    B.S. 100"""
    def __init__(self):
        self.variableNames = []
        #self.necessaryValues = {"Ausgang über zusammengefasste Daten" : ["MSTr", "MSE", "N", "α", "kleine Tabelle"], "Ausgang über alle Daten" : ["α", "große Tabelle"]}
        self.necessaryValues = {"Ausgang über alle Daten" : ["α", "Tukey's Test", "große Tabelle"]}
        self.extendedExplanations = {"Ausgang über alle Daten" : "'true' für nicht abgelehnt, 'false' für abgelehnt",
        "α" : "Signifikanzniveau des Tests",
        "Tukey's Test" : "'true' oder 'false'",
        "große Tabelle" : [["Stichprobe", "Wert"]]}

    def __call__(self, variableDict):
        if set(variableDict) == set(self.necessaryValues["Ausgang über alle Daten"]):
            looking_for = "Ausgang über alle Daten"
            foo = self.necessaryValues["Ausgang über alle Daten"]
        else:
            print(variableDict)
            raise Exception()

        for key in variableDict:
            if key == "Tukey's Test":
                if variableDict[key].lower() in ("true", "t", "j", "ja"):
                    variableDict[key] = True
                elif variableDict[key].lower() in ("false", "f", "n", "nein", ""):
                    variableDict[key] = False
                else:
                    raise Exception(f"Unexpected value '{variableDict[key]}' for {key}")

            elif key == "große Tabelle":
                #self.tableList = [[13.1,15.0,14.0,14.4,14.0,11.6], [16.3,15.7,17.2,14.9,14.4,17.2], [13.7,13.9,12.4,13.8,14.9,13.3], [15.7,13.7,14.4,16.0,13.9,14.7], [13.5,13.4,13.2,12.7,13.4,12.3]]
                #self.tableList = [[26.8, 27.9, 23.7, 25.0, 26.3, 24.8, 25.7, 24.5], [26.4, 24.2, 28.0, 26.9, 29.1], [25.7, 27.2, 29.9, 28.5, 29.4, 28.3]]
                #variableDict[key] = self.tableList
                #continue
                tableModel = variableDict[key]
                if type(tableModel) != tkt.TableModel:
                    raise Exception(f"{key} must not be of type 'tkt.TableModel' but is {type(tableModel)}")

                while all([tableModel.getCellRecord(tableModel.getRowCount()-1, k) == None for k in range(tableModel.getColumnCount())]):
                    tableModel.deleteRow(tableModel.getRowCount()-1)
                while all([tableModel.getCellRecord(k, tableModel.getColumnCount()-1) == None for k in range(tableModel.getRowCount())]):
                    tableModel.deleteColumn(tableModel.getColumnCount()-1)

                self.tableList = []
                for i in range(tableModel.getColumnCount()):
                    column = tableModel.getColCells(i)
                    self.tableList.append([])
                    for j in range(len(column)):
                        cell = column[j]
                        if cell == "":
                            cell = None
                        else:
                            try:
                                cell = float(cell.replace(",", "."))
                            except ValueError:
                                raise Exception(f"The value {cell} for {i}x{j} is unexpected!")
                            self.tableList[i].append(cell)
                variableDict[key] = self.tableList

            elif variableDict[key] == "":
                raise Exception(f"{key} must not be empty but is '{variableDict[key]}'")

            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")

        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang über alle Daten":
            print("Ausgang über alle Daten")

            Tabelle = variableDict["große Tabelle"]
            additionalResult = {}


            variableDict["N"] = len(variableDict["große Tabelle"])

            variableDict["Mᵢ"] = []
            variableDict["X\u0304ᵢ"] = []
            variableDict["Sᵢ²"] = []
            for column in Tabelle:
                amount = len(column)
                avg = sum([num for num in column]) / amount
                var = sum([(num-avg)**2 for num in column]) / (amount-1)

                variableDict["Mᵢ"].append(amount)
                variableDict["X\u0304ᵢ"].append(avg)
                variableDict["Sᵢ²"].append(var)

            X = sum(variableDict["X\u0304ᵢ"]) / variableDict["N"]
            n = sum(variableDict["Mᵢ"])
            SSTr = sum([variableDict["Mᵢ"][i] * (variableDict["X\u0304ᵢ"][i] - X)**2 for i in range(variableDict["N"])])
            MSTr = SSTr / (variableDict["N"] - 1)
            SSE = sum([(variableDict["Mᵢ"][i] - 1) * variableDict["Sᵢ²"][i] for i in range(variableDict["N"])])
            SST = SSTr + SSE
            MSE = SSE / (n - variableDict["N"])
            variableDict["v₁"] = variableDict["N"] - 1
            variableDict["v₂"] = n - variableDict["N"]
            fStern = MSTr / MSE
            PWert = 1 - f.cdf(fStern, variableDict["v₁"], variableDict["v₂"])

            krit = f.ppf(1-variableDict["α"], variableDict["v₁"], variableDict["v₂"])
            result = fStern >= krit

            """print("Mᵢ: ", variableDict["Mᵢ"])
            print("X\u0304ᵢ: ", variableDict["X\u0304ᵢ"])
            print("Sᵢ²: ", variableDict["Sᵢ²"])
            print("X: ", X)
            print("n: ", n)
            print()
            print("SSTr: ", SSTr)
            print("SSE: ", SSE)
            print("SST: ", SST)
            print()
            print("MSTr: ", MSTr)
            print("MSE: ", MSE)"""

            data1 = {
            "1" : {"Steuung" : "Treatments", "Freiheitsgrade" : f"N-1={variableDict['v₁']}", "Quadratsumme" : f"SSTr={SSTr}", "Mittlere Quadratsumme" : f"MSTr={MSTr}", "f-Wert" : f"MSTr/MSE={fStern}"},
            "2" : {"Steuung" : "Fehler", "Freiheitsgrade" : f"n-N={variableDict['v₂']}", "Quadratsumme" : f"SSE={SSE}", "Mittlere Quadratsumme" : f"MSE={MSE}", "f-Wert" : f""},
            "3" : {"Steuung" : "Total", "Freiheitsgrade" : f"n-1={n-1}", "Quadratsumme" : f"SST={SST}", "Mittlere Quadratsumme" : f"", "f-Wert" : f""}
            }

            data2 = {}
            for i in range(len(variableDict["Mᵢ"])):
                data2[f"{i+1}"] = {"Mᵢ" : variableDict["Mᵢ"][i], "X\u0304ᵢ" : variableDict["X\u0304ᵢ"][i], "Sᵢ²" : variableDict["Sᵢ²"][i]}

            if variableDict["Tukey's Test"] == True:
                score = []
                group = []
                for i in range(len(Tabelle)):
                    score = score + Tabelle[i]
                    group = group + [f"Probe {i+1}" for j in range(len(Tabelle[i]))]
                dataFrame = pd.DataFrame({"score" : score, "group"  : group})
                tukey = MultiComparison(data=dataFrame['score'], groups=dataFrame['group']).tukeyhsd(alpha=variableDict["α"])

                i = len(variableDict["X\u0304ᵢ"]) - 1
                j = 1
                couples = []
                var1 = list(tukey.reject)
                for avg in variableDict["X\u0304ᵢ"]:
                    var2 = var1[:i]
                    var1 = var1[i:]
                    k = j+1
                    for boo in var2:
                        if not boo:
                            couples.append({j, k})
                        k += 1
                    j += 1
                    i -= 1

                sortet_indizes = np.argsort(variableDict["X\u0304ᵢ"])
                subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                text = ""
                text2 = [" " * 5 * len(sortet_indizes) for i in range(len(sortet_indizes)-1)]
                for i in sortet_indizes:
                    text = text + " < " + (f"X\u0304{i+1}").translate(subscript)

                repl = "‾‾‾‾‾"
                leng = len(repl)
                i = 0
                sortet_indizes2 = sortet_indizes.copy()
                for ind1 in sortet_indizes:
                    sortet_indizes2 = sortet_indizes2[1:]
                    j = i
                    for ind2 in sortet_indizes2:
                        if j == i:
                            text2[i] = text2[i][:leng*j+3] + "‾‾" + text2[i][leng*(j+1):]
                            j += 1
                        try:
                            couples.remove({ind1+1, ind2+1})
                            text2[i] = text2[i][:leng*j] + repl + text2[i][leng*(j+1):]
                        except ValueError:
                            pass
                        j += 1
                    i += 1

                text3 = text[3:]
                for tex in text2:
                    text3 = text3 + "\n" + tex[3:]

                print(text3)
                additionalResult["\nTukey's Test"] = (str(tukey), )



            if not result:
                foo = {"data1" : data1,
                "data2" : data2,
                looking_for : "True, H₀ wird nicht abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit, "P-Wert" : PWert, "MSTr" : MSTr, "MSE" : MSE, "SSTr" : SSTr, "SSE" : SSE}
            elif result:
                foo = {"data1" : data1,
                "data2" : data2,
                looking_for : "False, H₀ wird abgelehnt", "Teststatistik" : fStern, "Kritischer Wert" : krit, "P-Wert" : PWert, "MSTr" : MSTr, "MSE" : MSE, "SSTr" : SSTr, "SSE" : SSE}

            #print({**foo, **additionalResult})
            return {**foo, **additionalResult}



#
class Regressionsanalyse(function):
    name = "Regressionsanalyse"
    description = """Regressionsanalyse mit Eingabe aller gemessener Werte
    Wichtig: '1. Regressorv.' entspricht y, alle weiteren '2. Regressorv.', '3. Regressorv.',... entsprechen x₁, x₂,...
    Model: y = β₀ + β₁x₁ + β₂x₂ + β₃x₃ + β₄x₄ + ... + βₖxₖ
    Die Funktionen der 2. Tabelle können auch auf Eingangswerte (Regressorvektoren) angewendet werden, die nicht Teil der Stichproben sind.
    Diese zusätzlichen Regressorvektoren können unterhalb der modellrelevanten RV in die Tabelle eingetragen werden. DAS FELD FÜR y MUSS DAFÜR LEER BLEIBEN.
    Nullhypothese der Parameterhypothesentests:    βᵢ = βH₀  ↔  -t<sub>krit<\sub> ≤ t* ≤ t<sub>krit<\sub>    (nicht getestet für βH₀ != 0)

    Erklärung des Resultats:

    1. Tabelle, 1. Zeile: Werte der Parameter β₀,...
    2. Zeile: Radien der Konfidenzintervalle der Parameter mit Konfidenznivevau 1-α (Signifikanzniveau α).
    Gegeben sind Radien r, das Intervall ergibt sich zu [β-r ; β+r]. β liegt mit einer Wahrscheinlichkeit von 1-α in diesem Intervall.
    Alternativ kann das Intervall so verstanden werden:
    Wenn der zu dem β gehöhrende Wert eines Eingangsvektors (Regressorvektors) x<sub>0<\sub> um 1 zunimmt, ändert sich das Resultat y mit einer Wahrscheinlichkeit von 1-α um einen Wert in dem Intervall.
    3. Zeile: Hypothesentests der einzelnen Parameter. Angen: H₀ kann nicht abgelehnt werden, abgel: H₀ wird abgelehnt.

    2. Tabelle, 1. Zeile: Werte der Modelausgänge y\u0302<sub>1...n<\sub> für die Regressorvektoren (Zeilen in der Tabelle) x<sub>1...n<\sub>.
    2./3. Zeile: Radien der Konfidenzintervalle/Prädiktionsintervalle der Modelausgänge mit Konfidenznivevau 1-α (Signifikanzniveau α).
    Gegeben sind Radien r, das Intervall ergibt sich zu [y\u0302-r ; y\u0302+r].
    y\u0302 liegt für die Eingangswerte, die im entsprechenden Regressorvektor x (Zeilen in der Tabelle) gegeben sind, mit einer Wahrscheinlichkeit von 1-α im Konfidenzintervall.
    Die nächste Messung von y\u0302 liegt für die Eingangswerte, die im entsprechenden Regressorvektor x gegeben sind, mit einer Wahrscheinlichkeit von 1-α im Prädiktionsintervall.


    B.S. 107"""
    def __init__(self):
        self.variableNames = []
        self.necessaryValues = {"Ausgang" : ["α", "βH₀", "Tukey's Test", "Parametergruppe", "große Tabelle"]}
        self.extendedExplanations = {"Ausgang" : "",
        "α" : "Signifikanzniveau der diversen Tests",
        "βH₀" : "Nullhypothese aller Hypothesentests der Parameter",
        "Tukey's Test" : "'true' oder 'false'",
        "Parametergruppe" : "Parameter, die zusammen getestet werden. Gegeben durch Indexe von 0 bis k, getrennt durch Beistriche.",
        "große Tabelle" : [["Regressorv.", "Wert"]]}

    def __call__(self, variableDict):
        if set(variableDict) == set(self.necessaryValues["Ausgang"]):
            looking_for = "Ausgang"
            foo = self.necessaryValues["Ausgang"]
        else:
            print(variableDict)
            raise Exception()

        for key in variableDict:
            if key == "Parametergruppe":
                if variableDict[key] in ("", None):
                    variableDict[key] = None
                else:
                    foo = variableDict[key].split(",")
                    foo = [x.strip() for x in foo]
                    try:
                        for i in range(len(foo)):
                            x = foo[i]
                            foo[i] = int(x)
                    except ValueError as e:
                        raise Exception(f"Unexpected value {x} for {key}")
                    variableDict[key] = foo

            elif key == "Tukey's Test":
                if variableDict[key].lower() in ("true", "t", "j", "ja"):
                    variableDict[key] = True
                    raise Exception("Tukeys Test ist für die Regressionsanalyse noch nicht implementiert.")
                elif variableDict[key].lower() in ("false", "f", "n", "nein", ""):
                    variableDict[key] = False
                else:
                    raise Exception(f"Unexpected value '{variableDict[key]}' for {key}")

            elif key == "große Tabelle":
                #self.tableList = [[1.0, 5.0, 9.0, 11.0, 3.0], [2.0, 6.0, 11.0, 25.0, 2.0], [3.0, 7.0, 12.0, 13.0, 9.0]]
                #self.tableList = [[0.12, 0.28, 0.55, 0.68, 0.85, 1.02, 1.15, 1.34, None, 1.29, None], [4.0, 8.7, 12.7, 19.1, 21.4, 24.6, 28.9, 29.8, 100, 30.5, 20]]
                #variableDict[key] = self.tableList
                #continue
                tableModel = variableDict[key]
                if type(tableModel) != tkt.TableModel:
                    raise Exception(f"{key} must not be of type 'tkt.TableModel' but is {type(tableModel)}")

                #Leere Zeilen und Spalten löschen
                while all([tableModel.getCellRecord(tableModel.getRowCount()-1, k) == None for k in range(tableModel.getColumnCount())]):
                    tableModel.deleteRow(tableModel.getRowCount()-1)
                while all([tableModel.getCellRecord(k, tableModel.getColumnCount()-1) == None for k in range(tableModel.getRowCount())]):
                    tableModel.deleteColumn(tableModel.getColumnCount()-1)

                self.tableList = []
                for i in range(tableModel.getColumnCount()):
                    column = tableModel.getColCells(i)
                    self.tableList.append([])
                    for j in range(len(column)):
                        cell = column[j]
                        if cell == "":
                            cell = None
                            #raise Exception(f"Cell must not be empty")
                        else:
                            try:
                                cell = float(cell.replace(",", "."))
                            except ValueError:
                                raise Exception(f"The value {cell} for {i}x{j} is unexpected!")
                        self.tableList[i].append(cell)
                variableDict[key] = self.tableList

            elif variableDict[key] == "":
                raise Exception(f"{key} must not be empty but is '{variableDict[key]}'")

            else:
                try:
                    variableDict[key] = float(variableDict[key].replace(",", "."))
                except ValueError:
                    raise Exception(f"The value {variableDict[key]} for {key} is unexpected!")

        return self.defaultCall(variableDict, looking_for)

    def defaultCall(self, variableDict, looking_for):
        if looking_for == "Ausgang":
            additionalResult = {}
            subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

            Tabelle = variableDict["große Tabelle"]
            #Tabelle = [[1.0, 5.0, 9.0, 11.0, 3.0], [2.0, 6.0, 11.0, 25.0, 2.0], [3.0, 7.0, 12.0, 13.0, 9.0]]#, [4.0, 8.0, 1.0, 4.0, 7.0]]
            #Tabelle = [[0.12, 0.28, 0.55, 0.68, 0.85, 1.02, 1.15, 1.34, 1.29], [4.0, 8.7, 12.7, 19.1, 21.4, 24.6, 28.9, 29.8, 30.5]]
            print("Tabelle: ", Tabelle)
            print("================================")
            Tabelle.insert(1, [1.0 for x in Tabelle[1]])

            all_y_list = Tabelle[0]
            all_x_list_list = Tabelle[1:]

            y_Vektor = np.array(all_y_list)
            X_Matrix = np.column_stack(all_x_list_list)
            NoneRows = np.where(y_Vektor == None)[0]
            NotNoneRows = np.where(y_Vektor != None)[0]
            y_Vektor = np.delete(arr=y_Vektor, obj=NoneRows, axis=0)
            additionalRVs = np.delete(arr=X_Matrix, obj=NotNoneRows, axis=0)
            X_Matrix = np.delete(arr=X_Matrix, obj=NoneRows, axis=0)
            X_Matrix_trans = X_Matrix.transpose()
            print("y_Vektor: ", y_Vektor)
            print("================================")
            print("X_Matrix: \n", X_Matrix)
            print("================================")
            print("X_Matrix transponiert: \n", X_Matrix_trans)
            print("================================")
            print("additionalRVs: \n", additionalRVs)
            print("================================")

            n = y_Vektor.shape[0]
            if n != X_Matrix.shape[0]:
                raise Exception()
            k = X_Matrix.shape[1]-1
            p = k + 1
            v = n - p
            print("n: ", n)
            print("================================")
            print("k: ", k)
            print("================================")
            print("p: ", p)
            print("================================")

            P_Matrix_vor_invertierung = np.dot(X_Matrix_trans, X_Matrix)
            print("P_Matrix vor Invertierung: \n", P_Matrix_vor_invertierung)
            print("================================")
            P_Matrix = np.linalg.inv(P_Matrix_vor_invertierung)
            print("P_Matrix: \n", P_Matrix)
            print("================================")

            foo = np.dot(P_Matrix, X_Matrix_trans)
            beta_dach = np.dot(foo, y_Vektor)
            print("beta_dach: ", beta_dach)
            print("================================")

            y_dach_Vektor = np.dot(X_Matrix, beta_dach)
            print("y_dach_Vektor: ", y_dach_Vektor)
            print("================================")
            y_strich_Skalar = sum([y_Vektor[i] for i in range(n)]) / n # Durchschnitt
            print("y_strich_Skalar: ", y_strich_Skalar)
            print("================================")

            e_Vektor = y_Vektor - y_dach_Vektor
            print("e_Vektor: ", e_Vektor)
            print("================================")

            SSE = np.dot(e_Vektor, e_Vektor)
            print("SSE: ", SSE)
            print("================================")
            SST = sum([(y_Vektor[i] - y_strich_Skalar)**2 for i in range(n)])
            print("SST: ", SST)
            print("================================")
            SSR = sum([(y_dach_Vektor[i] - y_strich_Skalar)**2 for i in range(n)])
            print("SSR: ", SSR)
            print("================================")
            MSR = SSR / (p - 1)
            print("MSR: ", MSR)
            print("================================")
            MSE = SSE / (n - p)
            print("MSE: ", MSE)
            print("================================")
            R_quadrat_Skalar = 1 - SSE / SST # 1 für perfektes Modell, 0 für nutzloses Modell, wird ungenau bei wachsender Modellordnung k
            print("Bestimmtheitsmaß R_quadrat_Skalar: ", R_quadrat_Skalar)
            print("================================")
            R_quadrat_angepasst_Skalar = 1 - (n - 1) / (n - p) * (1 - R_quadrat_Skalar) # 1 für perfektes Modell, 0 für nutzloses Modell
            print("R_quadrat_angepasst_Skalar: ", R_quadrat_angepasst_Skalar)
            print("================================")

            Varianz_e = SSE / (n - p)
            print("Varianz_e: ", Varianz_e)
            print("================================")
            Standardabweichung_e = math.sqrt(Varianz_e)
            print("Standardabweichung_e: ", Standardabweichung_e)
            print("================================")


            #Konfidenzintervall der einzelenen betas
            #Hypothesentests der einzelenen betas gegen Hβ₀
            data1 = {"1" : {}, "2" : {}, "3" : {}}
            for i in range(len(beta_dach)):
                r = t.ppf(variableDict["α"]/2, v) * Standardabweichung_e * math.sqrt(P_Matrix[i, i])
                tStern = (beta_dach[i] - variableDict["βH₀"]) / Standardabweichung_e / math.sqrt(P_Matrix[i, i])
                data1["1"][f"β{i}".translate(subscript)] = f"{beta_dach[i]}"
                data1["2"][f"β{i}".translate(subscript)] = f"± {abs(r)}"

                krit_links = t.ppf(variableDict["α"]/2, v) # links
                krit_rechts =  t.ppf(1 - variableDict["α"]/2, v) # rechts
                if krit_links <= tStern <= krit_rechts:
                    text = f"βH₀ angen.: {krit_links} ≤ {tStern} ≤ {krit_rechts}"
                elif krit_links > tStern:
                    text = f"βH₀ abgel.: {krit_links} \u2270 {tStern} ≤ {krit_rechts}"
                elif tStern > krit_rechts:
                    text = f"βH₀ abgel.: {krit_links} ≤ {tStern} \u2270 {krit_rechts}"
                else:
                    raise Exception()
                data1["3"][f"β{i}".translate(subscript)] = text
            print("data1: ", data1)
            print("================================")

            #Konfidenzintervall der Modelausgänge
            #Prädiktionsintervall der Modelausgänge
            data2 = {"1" : {}, "2" : {}, "3" : {}}
            for i in range(len(X_Matrix)):
                x = X_Matrix[i]
                print("x1: ", x)
                r1 = t.ppf(variableDict["α"]/2, v) * Standardabweichung_e * math.sqrt(np.dot(x, np.dot(P_Matrix, x.transpose())))
                r2 = t.ppf(variableDict["α"]/2, v) * Standardabweichung_e * math.sqrt(1 + np.dot(x, np.dot(P_Matrix, x.transpose())))
                data2["1"][f"y\u0302{i+1}".translate(subscript)] = f"{y_dach_Vektor[i]}"
                data2["2"][f"y\u0302{i+1}".translate(subscript)] = f"± {abs(r1)}"
                data2["3"][f"y\u0302{i+1}".translate(subscript)] = f"± {abs(r2)}"
            data2["1"][f""] = f""
            data2["2"][f""] = f""
            data2["3"][f""] = f""
            #Konfidenzintervall für die zusätzlichen Regressionsvektoren
            #Prädiktionsintervall für die zusätzlichen Regressionsvektoren
            for j in range(len(additionalRVs)):
                x = additionalRVs[j]
                print("x2: ", x)
                r1 = t.ppf(variableDict["α"]/2, v) * Standardabweichung_e * math.sqrt(np.dot(x, np.dot(P_Matrix, x.transpose())))
                r2 = t.ppf(variableDict["α"]/2, v) * Standardabweichung_e * math.sqrt(1 + np.dot(x, np.dot(P_Matrix, x.transpose())))

                data2["1"][f"Zusätzl. Wert {j+1}"] = f"{np.dot(beta_dach, x.transpose())}"
                data2["2"][f"Zusätzl. Wert {j+1}"] = f"± {abs(r1)}"
                data2["3"][f"Zusätzl. Wert {j+1}"] = f"± {abs(r2)}"
                #data2["1"][f"y\u0302{i+j+1}".translate(subscript)] = f"{y_dach_Vektor[i]}"
                #data2["2"][f"y\u0302{i+j+1}".translate(subscript)] = f"± {abs(r1)}"
                #data2["3"][f"y\u0302{i+j+1}".translate(subscript)] = f"± {abs(r2)}"
            print("data2: ", data2)
            print("================================")

            #Test auf globale Signifikanz des Models (ANOVA)
            fStern = MSR / MSE
            krit = f.ppf(1-variableDict["α"], p - 1, n - p)
            result = fStern >= krit
            data3 = {
            "1" : {"Steuung" : "Treatments", "Freiheitsgrade" : f"p-1={p - 1}", "Quadratsumme" : f"SSR={SSR}", "Mittlere Quadratsumme" : f"MSR={MSR}", "f-Wert" : f"MSR/MSE={fStern}"},
            "2" : {"Steuung" : "Fehler", "Freiheitsgrade" : f"n-p={n - p}", "Quadratsumme" : f"SSE={SSE}", "Mittlere Quadratsumme" : f"MSE={MSE}", "f-Wert" : f""},
            "3" : {"Steuung" : "Total", "Freiheitsgrade" : f"n-1={n-1}", "Quadratsumme" : f"SST={SST}", "Mittlere Quadratsumme" : f"", "f-Wert" : f""}
            }
            print("fStern: ", fStern)
            print("================================")
            print("krit: ", krit)
            print("================================")
            print("data3: ", data3)
            print("================================")
            if not result:
                additionalResult[looking_for] = "True, H₀ wird nicht abgelehnt, alle Parameter β unterscheiden sich nicht signifikant von 0"
                additionalResult["Teststatistik"] = fStern
                additionalResult["Kritischer Wert"] = krit
            elif result:
                additionalResult[looking_for] = "False, H₀ wird abgelehnt, mindestens ein Parameter β unterscheidet sich signifikant von 0"
                additionalResult["Teststatistik"] = fStern
                additionalResult["Kritischer Wert"] = krit


            #Tukey's Test
            if variableDict["Tukey's Test"] == True:
                score = []
                group = []
                for i in range(len(Tabelle)):
                    score = score + Tabelle[i]
                    group = group + [f"Probe {i+1}" for j in range(len(Tabelle[i]))]
                dataFrame = pd.DataFrame({"score" : score, "group"  : group})
                tukey = MultiComparison(data=dataFrame['score'], groups=dataFrame['group']).tukeyhsd(alpha=variableDict["α"])

                i = len(variableDict["X\u0304ᵢ"]) - 1
                j = 1
                couples = []
                var1 = list(tukey.reject)
                for avg in variableDict["X\u0304ᵢ"]:
                    var2 = var1[:i]
                    var1 = var1[i:]
                    k = j+1
                    for boo in var2:
                        if not boo:
                            couples.append({j, k})
                        k += 1
                    j += 1
                    i -= 1

                sortet_indizes = np.argsort(variableDict["X\u0304ᵢ"])
                subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
                text = ""
                text2 = [" " * 5 * len(sortet_indizes) for i in range(len(sortet_indizes)-1)]
                for i in sortet_indizes:
                    text = text + " < " + (f"X\u0304{i+1}").translate(subscript)

                repl = "‾‾‾‾‾"
                leng = len(repl)
                i = 0
                sortet_indizes2 = sortet_indizes.copy()
                for ind1 in sortet_indizes:
                    sortet_indizes2 = sortet_indizes2[1:]
                    j = i
                    for ind2 in sortet_indizes2:
                        if j == i:
                            text2[i] = text2[i][:leng*j+3] + "‾‾" + text2[i][leng*(j+1):]
                            j += 1
                        try:
                            couples.remove({ind1+1, ind2+1})
                            text2[i] = text2[i][:leng*j] + repl + text2[i][leng*(j+1):]
                        except ValueError:
                            pass
                        j += 1
                    i += 1

                text3 = text[3:]
                for tex in text2:
                    text3 = text3 + "\n" + tex[3:]

                print(text3)
                additionalResult["\nTukey's Test"] = (str(tukey), )

            if variableDict["Parametergruppe"] != None:
                #Test für Paramtergruppe
                print("\n\n================================")
                print("================================")
                print("Test für Paramtergruppe")
                print("================================")
                print("================================")
                Params = variableDict["Parametergruppe"]
                NotParams = [i for i in range(X_Matrix.shape[1]) if i not in Params]
                X_Matrix_r = np.delete(arr=X_Matrix, obj=NotParams, axis=1)
                print("X_Matrix_r: \n", X_Matrix_r)
                print("================================")
                X_Matrix_r_trans = X_Matrix_r.transpose()
                print("X_Matrix_r_trans: ", X_Matrix_r_trans)
                print("================================")
                P_Matrix_r = np.linalg.inv((np.dot(X_Matrix_r_trans, X_Matrix_r)))
                print("P_Matrix_r: ", P_Matrix_r)
                print("================================")
                beta_dach_r = np.dot(P_Matrix_r, np.dot(X_Matrix_r_trans, y_Vektor))
                print("beta_dach_r: ", beta_dach_r)
                print("================================")
                r = len(beta_dach_r)
                print("r: ", r)
                print("================================")
                SSE_r = np.dot(y_Vektor.transpose(), y_Vektor) - np.dot(beta_dach_r.transpose(), np.dot(X_Matrix_r.transpose(), y_Vektor))
                print("SSE_r: ", SSE_r)
                print("================================")
                fStern = ((SSE_r - SSE) / (p - r)) / (SSE / (n - p))
                print("fStern: ", fStern)
                print("================================")
                krit = f.ppf(1-variableDict["α"], p - r, n - p)
                result = fStern >= krit
                print("krit: ", krit)
                print("================================")
                print("data3: ", data3)
                print("================================")
                if not result:
                    additionalResult["Ausgang Paramtergruppe"] = f"True, H₀ wird nicht abgelehnt, alle Parameter β aus der Paramtergruppe unterscheiden sich nicht signifikant von 0, im Modell kann auf β ∈ {Params} verzichtet werden"
                    additionalResult["Teststatistik Paramtergruppe"] = fStern
                    additionalResult["Kritischer Wert Paramtergruppe"] = krit
                elif result:
                    additionalResult["Ausgang Paramtergruppe"] = f"False, H₀ wird abgelehnt, mindestens ein Parameter β aus der Paramtergruppe unterscheidet sich signifikant von 0, β ∈ {Params} beinflussen das Resultat signifikant"
                    additionalResult["Teststatistik Paramtergruppe"] = fStern
                    additionalResult["Kritischer Wert Paramtergruppe"] = krit










            print("additionalResult: ", additionalResult)
            print("================================")

            foo = {"data1" : data1,
            "data2" : data2,
            "data3" : data3}

            return {**foo, **additionalResult}










if __name__ == "__main__":
    import formelsammlung
