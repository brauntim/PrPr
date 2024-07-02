import json
import logging
import os
import main

# Erstelle Logs und JSONS Ordner, falls nicht vorhanden
os.makedirs('logs', exist_ok=True)
os.makedirs('jsons', exist_ok=True)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%m-%Y %H:%M:%S")

file_logger = logging.FileHandler('logs/searchEngine.log')
file_logger.setLevel(logging.DEBUG)
file_logger.setFormatter(formatter)

logger.addHandler(file_logger)


def list_json_files():
    json_files = [f for f in os.listdir('jsons') if f.endswith('.json')]
    if not json_files:
        logger.warning("Keine JSON-Dateien im JSONS-Ordner gefunden.")
        return []
    return json_files


def get_filename_from_user():
    while True:
        json_files = list_json_files()
        if not json_files:
            return None
        print("Verfügbare JSON-Dateien:")
        for idx, file in enumerate(json_files, start=1):
            print(f"{idx}: {file}")

        try:
            file_index = int(input("Bitte geben Sie die Nummer der zu durchsuchenden JSON-Datei an: ")) - 1
            if 0 <= file_index < len(json_files):
                return json_files[file_index]
            else:
                print("Ungültige Nummer. Bitte erneut versuchen.")
        except ValueError:
            print("Ungültige Eingabe. Bitte eine Nummer eingeben.")


def formatted_print(substance_data):
    formatted = json.dumps(substance_data, indent=4)
    print(formatted)


def filter_smiles(data):
    logger.info("Nach Smiles filtern")
    input_smiles = input("Smiles eingeben: ")
    if not isinstance(input_smiles, str):
        input_smiles = input("Smiles eingeben: ")
    contains = False
    count = 1

    for substance in data:
        if substance["smiles"] == input_smiles:
            print(f"\n{count}. Eintrag mit diesem Smiles: ")
            formatted_print(substance)
            contains = True
            count += 1
    print("\n")

    if not contains:
        print(f"Smiles \"{input_smiles}\" nicht in Datei enthalten.")


def filter_formular(data):
    logger.info("Nach Formel filtern")
    input_formular = input("Summenformel eingeben: ")
    if not isinstance(input_formular, str):
        input_formular = input("Summenformel eingeben: ")
    contains = False
    count = 1

    for substance in data:
        if substance["formula"] == input_formular:
            print(f"\n{count}. Eintrag mit dieser Summenformel: ")
            formatted_print(substance)
            contains = True
            count += 1
    print("\n")

    if not contains:
        print(f"Formel \"{input_formular}\" nicht in Datei enthalten.")


def filter_mass(data):
    logger.info("Nach Masse filtern")
    contains = False
    count = 1

    while True:
        try:
            input_mass_min = float(input("Minimale Masse eingeben: "))
            break
        except ValueError:
            print("Ungültige Eingabe. Bitte eine Zahl eingeben.")

    while True:
        try:
            input_mass_max = float(input("Maximale Masse eingeben: "))
            break
        except ValueError:
            print("Ungültige Eingabe. Bitte eine Zahl eingeben.")

    for substance in data:
        if input_mass_max >= float(substance["molecular_mass"]) >= input_mass_min:
            contains = True
            print(f"\n{count}. Eintrag mit dieser Masse: ")
            formatted_print(substance)
            count += 1
            print("\n")

    if not contains:
        print(f"\nKeine Substanz zwischen {input_mass_min} und {input_mass_max}\n")


def search_start(data):
    while True:
        print(
            "(1) Nach Smiles filtern \n(2) Nach Summenformel filtern \n(3) Nach Masse filtern \n"
            "(4) Datei auswählen \n(5) Beenden\n")

        while True:
            try:
                operation = int(input("Eingabe: "))
                break
            except ValueError:
                print("Ungültige Eingabe. Bitte eine ganze Zahl eingeben.")

        if operation == 1:
            filter_smiles(data)

        elif operation == 2:
            filter_formular(data)

        elif operation == 3:
            filter_mass(data)

        elif operation == 4:
            global filename
            filename = get_filename_from_user()
            with open(f'jsons/{filename}') as file:
                info = json.load(file)
            search_start(info)

        elif operation == 5:
            return operation

        else:
            print("Ungültige Eingabe. Bitte erneut versuchen.")


if __name__ == "__main__":
    filename = get_filename_from_user()
    if filename:
        try:
            with open(f'jsons/{filename}') as file:
                data = json.load(file)
        except FileNotFoundError:
            print("Datei nicht gefunden, bitte inkrementelles Laden wählen.")
            data = []

        print("Willkommen zu dieser Suchmaschine für Designerdrogen")
        print("______________________________________________________\n")

        end = 0

        while end != 5:
            end = search_start(data)

        print("Suchmaschine beendet")
    else:
        print("Keine JSON-Dateien gefunden. Beenden.")