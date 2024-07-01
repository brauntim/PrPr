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


def load_substances(filename):
    print("\n(1) Neu Laden \n(2) Neu hizugefügte Substanzen laden \n"
          "(3) Neu hizugefügte substanzen laden, Änderungen anpassen \n")

    input_load_option = int(input("Eingabe: "))

    if input_load_option == 1:
        logger.info("Alle Substanzen neu laden")
        # alle Substanzen werden neu von der Website geholt
        main.start_scraping(filename)
        print("Substanzen neu geladen.")

    elif input_load_option == 2:
        logger.info("Neu hinzugefuegte Substanzen laden")
        current_substances = json.load(open(f'jsons/{filename}'))

        # Neue Daten zum Vergleichen werden von der Website gescraped
        new_substances = load_new_substances()

        # Alte Substanzen werden mit den neuen Substanzen verglichen
        added, _ = compare_data(current_substances, new_substances)

        # Nur neue Substanzen werden an die Datei angehängt
        current_substances.extend(added)
        save_status(current_substances, filename)

        print("Neue Substanzen: ")
        formatted_print(added)
        print("\n")

    elif input_load_option == 3:

        logger.info("Neu hinzugefuegte und geaenderte Substanzen laden")
        current_substances = json.load(open(f'jsons/{filename}'))

        new_substances = load_new_substances()

        added, modified = compare_data(current_substances, new_substances)

        # hängt neue Substanzen an die alte Json-Datei
        for item in added:
            current_substances.append(item)

        # für alle Modifizierten Elemente wird, wenn die Smiles übereinstimmt, die alte Substanz überschrieben
        for item in modified:
            for i, old_item in enumerate(current_substances):
                # wenn die Smiles der Modifizierten Substanz übereinstimmt wird diese überschreiben
                if old_item['smiles'] == item['smiles']:
                    current_substances[i] = item

        save_status(current_substances, filename)

        print("Neue Substanzen: ")
        formatted_print(added)
        print("\n")

        print(f"Geänderte Substanzen: ")
        formatted_print(modified)
        print("\n")

    # Lösche new_substances.json nach dem Vergleich
    if os.path.exists("jsons/new_substances.json"):
        os.remove("jsons/new_substances.json")
        logger.info("new_substances.json gelöscht")

        
def compare_data(current, new):
    added = []
    modified = []

    # Smiles wird in in der neuen und alten Json als Primary-Key genutzt
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
    # zum Vergleichen werden alle Substanzen von der Website geholt
    main.start_scraping("new_substances.json")
    with open(f"jsons/new_substances.json") as new_file:
        new_substances = json.load(new_file)

    if new_substances:
        return new_substances
    else:
        print("Fehler beim Abrufen der Webseite")
        return None


def save_status(substances, filename):
    # aktuellen Stand der Substanzen speichern
    logger.info("Datei gespeichert")
    with open(f'jsons/{filename}', 'w') as f:
        json.dump(substances, f, indent=4)


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


def search_start(data, filename):
    while True:
        print(
            "(1) Inkrementelles Laden \n(2) Nach Smiles filtern \n(3) Nach Summenformel filtern \n(4) Nach Masse filtern \n"
            "(5) Beenden\n")

        while True:
            try:
                operation = int(input("Eingabe: "))
                break
            except ValueError:
                print("Ungültige Eingabe. Bitte eine ganze Zahl eingeben.")

        if operation == 1:
            load_substances(filename)
            with open(f'jsons/{filename}') as file:
                data = json.load(file)  # Refresh data after loading substances

        elif operation == 2:
            filter_smiles(data)

        elif operation == 3:
            filter_formular(data)

        elif operation == 4:
            filter_mass(data)

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
        print("______________________________________________________")

        end = 0

        while end != 5:
            end = search_start(data, filename)

        print("Suchmaschine beendet")
    else:
        print("Keine JSON-Dateien gefunden. Beenden.")
