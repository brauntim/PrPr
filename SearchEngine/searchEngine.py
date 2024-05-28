import json
import logging
from main import SubstanceExtractor


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%m- %Y %H:%M:%S")

file_logger = logging.FileHandler('searchEngine.log')
file_logger.setLevel(logging.DEBUG)
file_logger.setFormatter(formatter)

logger.addHandler(file_logger)



def load_substances():

    print("\n(1) Neu Laden \n(2) Neu hizugefügte Substanzen laden \n"
          "(3) Neu hizugefügte substanzen laden, Änderungen anpassen \n")

    input_load_option = input("Eingabe: ")

    if(input_load_option == "1"):
        url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
        extractor = SubstanceExtractor(url)
        extractor.run()
        print("Datenextraktion abgeschlossen. Die Ergebnisse sind in substancesV4.json gespeichert.")

    elif(input_load_option == "2"):
        print("2")
        logger.debug("inputLoading")
    elif(input_load_option == "3"):
        print("3")

def filter_smiles():
    input_smiles = input("Smiles eingeben: ")
    contains = False

    for substance in data:
        if substance["smiles"] == input_smiles:
            print("\nEintrag mit dieser Smiles gefunden:")
            print(substance)
            contains = True
    print("\n")

    if not contains:
        print("Smiles \"" + input_smiles + "\" nicht in Datei enthalten.")


def filter_formular():
    input_formular = input("Summenformel eingeben: ")
    contains = False

    for substance in data:
        if substance["formular"] == input_formular:

            print("\nEintrag mit dieser Summenformel gefunden")
            print(substance)
            contains = True
    print("\n")

    if not contains:
        print("Formel \"" + input_formular + "\" nicht in Datei enthalten.")


def filter_mass():
    input_mass_min = input("minimale Masse eingeben: ")
    input_mass_max = input("maximale Masse eingeben: ")

    logger.debug("nach Masse filtern")

    for substance in data:
        if input_mass_max >= substance["molecular_mass"] >= input_mass_min:
            print("Eintrag mit dieser Masse:")
            print(substance)
            print("\n")


def search_start():
    while True:
        print("(1) Inkrementelles Laden \n(2) Nach Smiles filtern \n(3) Nach Summenformel filtern \n(4) Nach Masse filtern \n"
              "(5) Beenden\n")

        operation = input("Eingabe: ")
        operation = int(operation)

        if operation == 1:
            load_substances()

            return operation

        elif operation == 2:
            filter_smiles()

            return operation

        elif operation == 3:
            filter_formular()

            return operation

        elif operation == 4:
            filter_mass()

            return operation

        elif operation == 5:
            return operation
        else:
            operation = input("Bitte eingabe wiederholen: ")
            return operation



with open("substancesV4.json") as file:
    data = json.load(file)

print("Willkommen zu dieser Suchmaschine für Designerdrogen")
print("______________________________________________________")

end = 0

while end != 5:
    end = search_start()

print("Suchmaschine beendet")

