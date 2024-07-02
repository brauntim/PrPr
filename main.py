#
#   webscraper for https://www.policija.si/apps/nfl_response_web/seznam.php
#
import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging
import os

# Erstelle Logs Ordner, falls nicht vorhanden
os.makedirs('logs', exist_ok=True)
os.makedirs('jsons', exist_ok=True)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Log-Format für Dateilogs
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%m-%Y %H:%M:%S")

# Dateilogs konfigurieren
file_logger = logging.FileHandler('logs/scraper.log')
file_logger.setLevel(logging.DEBUG)
file_logger.setFormatter(formatter)

# Konsolenlogs konfigurieren
console_logger = logging.StreamHandler()
console_logger.setLevel(logging.INFO)
console_logger.setFormatter(formatter)

# Füge beide Logger hinzu
logger.addHandler(file_logger)
logger.addHandler(console_logger)


class SubstanceExtractor:
    def __init__(self, url):
        self.url = url
        self.substances = []

    # HTML-Daten von der URL abrufen
    def fetch_data(self):
        logger.info("Scraping started")
        response = requests.get(self.url)
        response.raise_for_status()  # Fehler auslösen, wenn die Anfrage fehlschlägt
        return response.text

    # Konvertiere SMILES in kanonisches SMILES-Format
    def convert_to_canonical_smiles(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles

    # Daten von OPSIN-Datenbank abrufen
    def fetch_opsin_data(self, iupac_name):
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{iupac_name}.json"
        response = requests.get(opsin_url)
        if response.status_code == 200:
            data = response.json()
            smiles_str = data.get('smiles', '')
            inchi_str = data.get('stdinchi', '')
            smiles_str = self.convert_to_canonical_smiles(smiles_str)
            return smiles_str, inchi_str
        else:
            return "", ""

    # Eine einzelne Zeile der HTML-Tabelle parsen
    def parse_row(self, row):
        columns = row.find_all('td')
        if len(columns) < 17:  # Überprüfe, ob die Zeile genügend Spalten hat
            return None

        # Extrahiere Daten aus den Tabellenzellen
        category = columns[0].text.strip()
        name = columns[1].text.strip()
        pdf_link = columns[1].find('a')['href'] if columns[1].find('a') else self.url
        pdf_link = requests.compat.urljoin(self.url, pdf_link)
        iupac_name = columns[3].text.strip()
        formula = columns[5].text.strip()
        molecular_mass = columns[6].text.strip()
        std_inchi_key = columns[12].text.strip()
        date_of_entry = columns[15].text.strip()
        report_updated = columns[16].text.strip()

        last_changed_at = report_updated if report_updated else date_of_entry
        deleted = '<strike>' in name  # Überprüfe, ob die Substanz als gelöscht markiert ist

        try:
            # Versuche, SMILES und InChI von OPSIN zu holen
            smiles, std_inchi = self.fetch_opsin_data(iupac_name)
        except Exception as e:
            logger.error(f"Error fetching data for IUPAC name {iupac_name}: {e}")
            smiles, std_inchi = "", ""

        # Erstelle ein Dictionary mit den extrahierten und abgerufenen Daten
        substance_data = {
            "version": "1.0",
            "smiles": smiles,
            "names": [name],
            "iupac_names": [iupac_name],
            "formula": formula,
            "inchi": std_inchi,
            "inchi_key": std_inchi_key,
            "molecular_mass": float(molecular_mass) if molecular_mass else None,
            "cas_num": None,
            "categories": [category],
            "source": {
                "name": "Policija",
                "url": pdf_link
            },
            "validated": None,
            "deleted": deleted,
            "last_modified": last_changed_at,
            "details": {}
        }
        return substance_data

    # Den gesamten HTML-Inhalt parsen
    def parse_html(self, html):
        soup = BeautifulSoup(html, 'html.parser')
        for index, row in enumerate(soup.select('table tr'), start=1):
            logger.info("Element " + str(index) + " scraped.")
            try:
                substance_data = self.parse_row(row)
                if substance_data:
                    self.substances.append(substance_data)
            except Exception as e:
                print(f"Error parsing row: {index}")
                logger.error(f"Fehler beim  Scrapen: {index}")

    # Unicode-Zeichen ersetzen
    def replace_unicode_characters(self, data):
        replacements = {
            "\u2010": "-",  # Hyphen
            "\u200b": "",  # Zero Width Space
        }
        if isinstance(data, str):
            for key, value in replacements.items():
                data = data.replace(key, value)
            return data
        elif isinstance(data, list):
            return [self.replace_unicode_characters(item) for item in data]
        elif isinstance(data, dict):
            if 'iupac_names' in data:
                data['iupac_names'] = self.replace_unicode_characters(data['iupac_names'])
            return {key: self.replace_unicode_characters(value) for key, value in data.items()}
        else:
            return data

    # Daten validieren
    def validate_data(self):
        for substance in self.substances:
            # Hole den SMILES-String der Substanz, falls vorhanden, ansonsten leere Zeichenkette
            smiles = substance.get("smiles", "")
            # Wenn ein SMILES-String vorhanden ist
            if smiles:
                # Erstelle ein Molekülobjekt aus dem SMILES-String
                molecule = Chem.MolFromSmiles(smiles)
                if molecule:
                    # Berechne die molekulardes Moleküls
                    calculated_molecular_mass = Descriptors.MolWt(molecule)
                    # Berechne die chemische Formel des Moleküls
                    calculated_formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)

                    try:
                        molecular_mass_valid = abs(
                            float(substance["molecular_mass"]) - calculated_molecular_mass) < 0.99
                    except ValueError:
                        molecular_mass_valid = False
                    # Überprüfe, ob die angegebene Formel mit der berechneten übereinstimmt
                    formula_valid = substance["formula"] == calculated_formula
                    substance["validated"] = molecular_mass_valid and formula_valid
                else:
                    # Wenn das Molekül nicht erstellt werden konnte, setze die Validität auf None
                    substance["validated"] = None
            else:
                # Wenn kein SMILES-String vorhanden ist, setze die Validität auf None
                substance["validated"] = None


    # Daten in eine JSON-Datei speichern
    def save_to_json(self, filename):
        updated_substances = self.replace_unicode_characters(self.substances)
        with open(f"jsons/{filename}", 'w', encoding='utf-8') as json_file:
            json.dump(updated_substances, json_file, indent=4)
        logger.info(f"Saving substances to {filename}")
        logger.info("Scraping finished")

    # Scraper ausführen
    def run(self, filename):
        html = self.fetch_data()
        self.parse_html(html)
        self.validate_data()
        self.save_to_json(filename)


def load_options(filename):
    print("\n(1) Neu Laden \n(2) Neu hizugefügte Substanzen laden \n"
          "(3) Neu hizugefügte substanzen laden, Änderungen anpassen \n")

    input_load_option = int(input("Eingabe: "))

    if input_load_option == 1:
        logger.info("Alle Substanzen neu laden")
        # alle Substanzen werden neu von der Website geholt
        start_scraping(filename)
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
    start_scraping("new_substances.json")
    with open(f"jsons/new_substances.json") as new_file:
        new_substances = json.load(new_file)

    if new_substances:
        return new_substances
    else:
        print("Fehler beim Abrufen der Webseite")
        return None


def save_status(substances, filename):
    # aktuellen Stand der Substanzen speichern
    logger.info(f"Datei {filename} gespeichert")
    with open(f'jsons/{filename}', 'w') as f:
        json.dump(substances, f, indent=4)


def formatted_print(substance_data):
    formatted = json.dumps(substance_data, indent=4)
    print(formatted)


# Scraper starten
def start_scraping(filename):
    url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
    extractor = SubstanceExtractor(url)
    extractor.run(filename)
    print(f"Datenextraktion abgeschlossen. Die Ergebnisse sind in jsons/{filename} gespeichert.")


if __name__ == "__main__":
    filename = 'ProjektPolicija.json'
    try:
        with open("jsons/"+filename) as file:
            data = json.load(file)
            logger.info("Ladeoptionen aufrufen")
        load_options(filename)

    except FileNotFoundError:
        logger.info("Website wird neu gescraped")
        start_scraping(filename)
