#
#   webscraper for https://www.policija.si/apps/nfl_response_web/seznam.php
#

import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import Descriptors

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%m- %Y %H:%M:%S")

file_logger = logging.FileHandler('scraper.log')
file_logger.setLevel(logging.DEBUG)
file_logger.setFormatter(formatter)

conole_logger = logging.StreamHandler()
conole_logger.setLevel(logging.INFO)
conole_logger.setFormatter(formatter)

logger.addHandler(file_logger)
logger.addHandler(conole_logger)


# Diese Klasse handhabt das Abrufen, Parsen und Speichern von Substanzdaten von einer gegebenen URL.
class SubstanceExtractor:

    # Initialisiert eine neue Instanz der SubstanceExtractor-Klasse.
    def __init__(self, url):
        self.url = url
        self.substances = []

    # GET-Anfrage an die Website und gibt HTML-Inhalt zurück
    def fetch_data(self):
        logger.info("Scraping started")
        response = requests.get(self.url)
        response.raise_for_status()
        return response.text

    def convert_to_canonical_smiles(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles

    # Funktion zum Abrufen von StdInChI und SMILES von OPSIN-API
    def fetch_opsin_data(self, iupac_name):
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{iupac_name}.json"
        response = requests.get(opsin_url)
        if response.status_code == 200:
            data = response.json()
            smiles_str = data.get('smiles', '')
            inchi_str = data.get('stdinchi', '')
            # Konvertieren des SMILES-Strings in einen kanonischen SMILES-String
            smiles_str = self.convert_to_canonical_smiles(smiles_str)
            return smiles_str, inchi_str
        else:
            return "", ""

    # Funktion parst eine einzelne Zeile und gibt die Substanzdaten zurück.
    def parse_row(self, row):
        columns = row.find_all('td')
        if len(columns) < 17:  # Sicherstellen, dass die Zeile genügend Spalten hat
            return None

        # Extrahieren der Daten aus den entsprechenden Spalten
        category = columns[0].text.strip()
        name = columns[1].text.strip()
        iupac_name = columns[3].text.strip()
        formular = columns[5].text.strip()
        molecular_mass = columns[6].text.strip()
        std_inchi_key = columns[12].text.strip()
        date_of_entry = columns[15].text.strip()
        report_updated = columns[16].text.strip()

        # Bestimmen des Datums der letzten Änderung
        if report_updated:
            last_changed_at = report_updated
        else:
            last_changed_at = date_of_entry

        # Überprüfen, ob die Substanz gelöscht wurde
        deleted = '<strike>' in name

        # Abrufen von SMILES und InChI über OPSIN-API mittels IUPAC
        try:
            smiles, std_inchi = self.fetch_opsin_data(iupac_name)
        except Exception as e:
            logger.error(f"Error fetching data for IUPAC name {iupac_name}: {e}")
            print(f"Error fetching data for IUPAC name {iupac_name}: {e}")
            smiles, std_inchi = "", ""

        # Erstellen eines Dictionaries mit den extrahierten Daten
        substance_data = {
            "smiles": smiles,
            "names": name,
            "iupac_name": iupac_name,
            "formular": formular,
            "molecular_mass": molecular_mass,
            "Inchi": std_inchi,
            "InchiKey": std_inchi_key,
            "cas_num": "",  # Nicht vorhanden
            "category": category,
            "source_name": "Policija",
            "source_url": self.url,
            "valid": None,
            "deleted": deleted,
            "last_changed_at": last_changed_at,
            "version": "1.0",
            "details": {
                # Details unseres Teams
            }
        }
        return substance_data

    # Parst den HTML-Inhalt und extrahiert die Substanzdaten
    def parse_html(self, html):
        soup = BeautifulSoup(html, 'html.parser')
        for index, row in enumerate(soup.select('table tr'), start=1):
            logger.info("Element " + str(index) + " scraped.")
            if index >= 20:
                break
            try:
                substance_data = self.parse_row(row)
                if substance_data:
                    self.substances.append(substance_data)
            except Exception as e:
                print(f"Error parsing row: {index}")
                logger.error(f"Fehler beim  Scrapen: {index}")

    # Funktion zum Ersetzen von "\u2010" durch "-" in JSON Datei
    def replace_unicode_characters(self, data):
        replacements = {
            "\u2010": "-",
            "\u200b": "",
        }
        if isinstance(data, str):
            for key, value in replacements.items():
                data = data.replace(key, value)
            return data
        elif isinstance(data, list):
            return [self.replace_unicode_characters(item) for item in data]
        elif isinstance(data, dict):
            if 'iupac_name' in data:
                data['iupac_name'] = self.replace_unicode_characters(data['iupac_name'])
            return {key: self.replace_unicode_characters(value) for key, value in data.items()}
        else:
            return data

    # Funktion zur Berechnung von Formelgewicht und Summenformel aus SMILES
    def validate_data(self, substance):
        for substance in self.substances:
            smiles = substance.get("smiles", "")
            if smiles:
                molecule = Chem.MolFromSmiles(smiles)
                if molecule:
                    calculated_molecular_mass = Descriptors.MolWt(molecule)
                    calculated_formular = Chem.rdMolDescriptors.CalcMolFormula(molecule)
                    try:
                        molecular_mass_valid = abs(float(substance[
                                                             "molecular_mass"]) - calculated_molecular_mass) < 0.05  # Absolute Differenz, Wert notfalls anpassen für die Toleranz
                    except ValueError:
                        molecular_mass_valid = False

                    formular_valid = substance["formular"] == calculated_formular
                    substance["valid"] = molecular_mass_valid and formular_valid
                else:
                    substance["valid"] = None
            else:
                substance["valid"] = None

                # Speichert die extrahierten Substanzdaten in einer JSON-Datei.

    def save_to_json(self, filename):
        updated_substances = self.replace_unicode_characters(self.substances)
        with open(filename, 'w', encoding='utf-8') as json_file:
            json.dump(updated_substances, json_file, indent=4)
        logger.info(f"Saving substances to {filename}")
        logger.info("Scraping finished")

    # Führt den gesamten Extraktionsprozess aus.
    def run(self):
        html = self.fetch_data()
        self.parse_html(html)
        self.validate_data(self.substances)
        self.save_to_json('substances.json')



def start_scraping():
    url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
    extractor = SubstanceExtractor(url)
    extractor.run()
    print("Datenextraktion abgeschlossen. Die Ergebnisse sind in substancesV4.json gespeichert.")


if __name__ == "__main__":
    start_scraping()
