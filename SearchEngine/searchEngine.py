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

    input_load_option = input("Eingabe: ")

    if input_load_option == "1":
        main.start_scraping()

    elif input_load_option == "2":

        save_substances(data)

        main.start_scraping()
        current_substances = json.load(open("current_substances.json"))

        with open("substancesV4.json") as new_file:
            new_substances = json.load(new_file)

        added, _ = compare_data(current_substances, new_substances)

        current_substances.update(added)
        save_substances(current_substances)
        print(current_substances)

        logger.debug("inputLoading")

    elif input_load_option == "3":
        print("3")


def save_substances(data):
    with open("current_substances.json", "w") as save_file:
        json.dump(data, save_file, indent=4)


def compare_data(current, new):
    added = {}
    modified = {}

    current_smiles = {details['smiles']: details for details in current.values()}
    new_smiles = {details['smiles']: details for details in new.values()}

    for smiles, new_details in new_smiles.items():
        if smiles not in current_smiles:
            added[smiles] = new_details

        elif current_smiles[smiles] != new_details:
            modified[smiles] = new_details

        else:
            old_details = current[id]
        if (
                old_details["smiles1"] != new_details["smiles"] or
                old_details["names"] != new_details["names"] or
                old_details["iupac_names"] != new_details["iupac_names"] or
                old_details["formular"] != new_details["formular"] or
                old_details["molecular_mass"] != new_details["molecular_mass"] or
                old_details["Inchi"] != new_details["Inchi"] or
                old_details["InchiKey"] != new_details["InchiKey"] or
                old_details["cas_num"] != new_details["cas_num"] or
                old_details["category"] != new_details["category"] or
                old_details["source_name"] != new_details["source_name"] or
                old_details["source_url"] != new_details["source_url"] or
                old_details["valid"] != new_details["valid"] or
                old_details["deleted"] != new_details["deleted"] or
                old_details["last_changed_at"] != new_details["last_changed_at"] or
                old_details["version"] != new_details["version"]
        ):
            modified[id] = new_details

    return added, modified


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
        print(
            "(1) Inkrementelles Laden \n(2) Nach Smiles filtern \n(3) Nach Summenformel filtern \n(4) Nach Masse filtern \n"
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

            with open("test.json", "w") as testfile:
                json.dump(data, testfile, indent=4)

            return operation
        else:
            operation = input("Bitte eingabe wiederholen: ")
            return operation


if __name__ == "__main__":
    with open("current_substances.json") as file:
        data = json.load(file)

    print("Willkommen zu dieser Suchmaschine für Designerdrogen")
    print("______________________________________________________")

    end = 0

    while end != 5:
        end = search_start()

    print("Suchmaschine beendet")
