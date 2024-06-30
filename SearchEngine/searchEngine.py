import json
import logging

import main

# from main import SubstanceExtractor

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

    input_load_option = int(input("Eingabe: "))

    if input_load_option == 1:
        main.start_scraping("substances.json")
        print("Substanzen neu geladen.")

    elif input_load_option == 2:

        current_substances = json.load(open("substances.json"))
        print("Current: ", current_substances)

        new_substances = load_new_substances()
        print("New: ", new_substances)

        added, _ = compare_data(current_substances, new_substances)

        current_substances.extend(added)
        save_status(current_substances)
        #print(f"Neue Substanzen: {formatted_print(added)}")
        print("Neue Substanzen: ")
        formatted_print(added)
        logger.debug("inputLoading")

    elif input_load_option == 3:

        current_substances = json.load(open("substances.json"))
        #print("Current: ", current_substances)

        new_substances = load_new_substances()
        #print("New: ", new_substances)

        added, modified = compare_data(current_substances, new_substances)



        for item in added:
            current_substances.append(item)
        for item in modified:
            for i, old_item in enumerate(current_substances):
                if old_item['smiles'] == item['smiles']:
                    current_substances[i] = item

        save_status(current_substances)

        print("Neue Substanzen: ")
        formatted_print(added)
        print(f"Geänderte Substanzen: ")
        formatted_print(modified)



def compare_data(current, new):
    added = []
    modified = []

    current_smiles = {details['smiles']: details for details in current}
    new_smiles = {details['smiles']: details for details in new}

    for smiles, new_details in new_smiles.items():
        if smiles not in current_smiles:
            added.append(new_details)

        else:
            old_details = current_smiles[smiles]

            if old_details != new_details:
                modified.append(new_details)

    return added, modified


def load_new_substances():
    main.start_scraping("new_substances.json")
    with open("new_substances.json") as new_file:
        new_substances = json.load(new_file)

    if new_substances:
        return new_substances
    else:
        print(f"Fehler beim Abrufen der Webseite: Status Code {new_substances.status_code}")
        return None

def save_status(substances):

    with open("substances.json", 'w') as f:
        json.dump(substances, f, indent=4)

def formatted_print(substance_data):
    formatted = json.dumps(substance_data, indent=4)
    print(formatted)

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
    input_mass_min = float(input("minimale Masse eingeben: "))
    input_mass_max = float(input("maximale Masse eingeben: "))

    logger.debug("nach Masse filtern")

    for substance in data:
        if input_mass_max >= substance["molecular_mass"] >= input_mass_min:
            print("Eintrag mit dieser Masse:")
            print(substance)
            print("\n")


def search_start():
    while True:
        print(
            "(1) Inkrementelles Laden \n(2) Nach Smiles filtern \n(3) Nach Summenformel filtern \n(4) Nach Masse filtern \n"
            "(5) Beenden\n")

        operation = int(input("Eingabe: "))

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

            with open("test.json", "w") as testfile:
                json.dump(data, testfile, indent=4)

            return operation
        else:
            operation = input("Bitte eingabe wiederholen: ")
            return operation


if __name__ == "__main__":
    try:
        with open("substances.json") as file:
            data = json.load(file)
    except FileNotFoundError:
        main.start_scraping("substances.json")
        with open("substances.json") as file:
            data = json.load(file)

    print("Willkommen zu dieser Suchmaschine für Designerdrogen")
    print("______________________________________________________")

    end = 0

    while end != 5:
        end = search_start()

    print("Suchmaschine beendet")
