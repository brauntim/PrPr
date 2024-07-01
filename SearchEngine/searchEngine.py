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
        logger.info("Alle Substanzen neu laden")
        # alle Substanzen werden neu von der Website geholt
        main.start_scraping(filename)
        print("Substanzen neu geladen.")

    elif input_load_option == 2:
        logger.info("Neu hinzugefuegte Substanzen laden")
        current_substances = json.load(open(filename))

        # neue Daten zum Vergleichen werden von der Website gescraped
        new_substances = load_new_substances()

        # alte Substanzen werden mit den neuen Substanzen verglichen
        added, _ = compare_data(current_substances, new_substances)

        # nur neue Substanzen werden an die Datei angehängt
        current_substances.extend(added)
        save_status(current_substances)

        print("Neue Substanzen: ")
        formatted_print(added)
        print("\n")

    elif input_load_option == 3:

        logger.info("Neu hinzugefuegte und geaenderte Substanzen laden")
        current_substances = json.load(open(filename))

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

        save_status(current_substances)

        print("Neue Substanzen: ")
        formatted_print(added)
        print("\n")

        print(f"Geänderte Substanzen: ")
        formatted_print(modified)
        print("\n")


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
    with open("new_substances.json") as new_file:
        new_substances = json.load(new_file)

    if new_substances:
        return new_substances
    else:
        print("Fehler beim Abrufen der Webseite")
        return None


def save_status(substances):
    # aktuellen Stand der Substanzen speichern
    logger.info("Datei gespeichert")
    with open(filename, 'w') as f:
        json.dump(substances, f, indent=4)


def formatted_print(substance_data):
    formatted = json.dumps(substance_data, indent=4)
    print(formatted)


def filter_smiles():
    logger.info("Nach Smiles filtern")
    input_smiles = input("Smiles eingeben: ")
    if not isinstance(input_smiles, str):
        input_smiles = input("Smiles eingeben: ")
    contains = False

    for substance in data:
        # geht durch alle Substanzen und checkt, wo die Smiles übereinstimmt
        if substance["smiles"] == input_smiles:
            print("\nEintrag mit dieser Smiles gefunden:")
            formatted_print(substance)
            contains = True
    print("\n")

    if not contains:
        print("Smiles \"" + input_smiles + "\" nicht in Datei enthalten.")


def filter_formular():
    logger.info("Nach Formel filtern")
    input_formular = input("Summenformel eingeben: ")
    if not isinstance(input_formular, str):
        input_formular = input("Summenformel eingeben: ")
    contains = False

    for substance in data:
        # geht durch alle Substanzen und checkt, wo die Formel übereinstimmt
        if substance["formula"] == input_formular:
            print("\nEintrag mit dieser Summenformel gefunden")
            formatted_print(substance)
            contains = True
    print("\n")

    if not contains:
        print("Formel \"" + input_formular + "\" nicht in Datei enthalten.")


def filter_mass():
    logger.info("nach Masse filtern")
    contains = False
    while True:
        try:
            input_mass_min = float(input("minimale Masse eingeben: "))
            break  # Wenn die Konvertierung erfolgreich ist, die Schleife beenden
        except ValueError:
            print("Ungültige Eingabe. Bitte eine Zahl eingeben.")

    while True:
        try:
            input_mass_max = float(input("maximale Masse eingeben: "))
            break  # Wenn die Konvertierung erfolgreich ist, die Schleife beenden
        except ValueError:
            print("Ungültige Eingabe. Bitte eine Zahl eingeben.")

    for substance in data:
        # checkt für jede Substanz ob sie in dem jeweiligen Massebereich liegt
        if input_mass_max >= float(substance["molecular_mass"]) >= input_mass_min:
            contains = True
            print("Eintrag mit dieser Masse:")
            formatted_print(substance)
            print("\n")

    if not contains:
        print(f"\nKeine Substanz zwischen {input_mass_min} und {input_mass_max}\n")


def search_start():
    while True:
        print(
            "(1) Inkrementelles Laden \n(2) Nach Smiles filtern \n(3) Nach Summenformel filtern \n(4) Nach Masse filtern \n"
            "(5) Beenden\n")

        while True:
            try:
                operation = int(input("Eingabe: "))
                break  # Wenn die Konvertierung erfolgreich ist, die Schleife beenden
            except ValueError:
                print("Ungültige Eingabe. Bitte eine ganze Zahl eingeben.")

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


if __name__ == "__main__":
    filename = "Tim_Jonas_Policija.json"
    try:
        with open(filename) as file:
            data = json.load(file)
    except FileNotFoundError:
        main.start_scraping(filename)
        with open(filename) as file:
            data = json.load(file)
        logger.info("Neu gescraped")

    print("Willkommen zu dieser Suchmaschine für Designerdrogen")
    print("______________________________________________________")

    end = 0

    while end != 5:
        end = search_start()

    print("Suchmaschine beendet")
