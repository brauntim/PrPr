#
#   webscraper for https://www.policija.si/apps/nfl_response_web/seznam.php
#

import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import Descriptors
# import time

# Diese Klasse handhabt das Abrufen, Parsen und Speichern von Substanzdaten von einer gegebenen URL.
class SubstanceExtractor:

    # Initialisiert eine neue Instanz der SubstanceExtractor-Klasse.
    def __init__(self, url):
        self.url = url
        self.substances = []

    # GET-Anfrage an die Website und gibt HTML-Inhalt zurück
    def fetch_data(self):
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
            # print(smiles)
            # print(std_inchi)
            # time.sleep(10)
        except Exception as e:
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
            "cas_num": "",  # Wird aus anderen Informationen abgeleitet
            "category": category,
            "source_name": "Policija",
            "source_url": self.url,
            "valid": None,  # Nicht genügend Informationen, um dies zu bestimmen
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
        for index,row in enumerate(soup.select('table tr'), start=1):
            try:
                substance_data = self.parse_row(row)
                if substance_data:
                    self.substances.append(substance_data)
            except Exception as e:
                print(f"Error parsing row: {index}" )

    # Funktion zum Ersetzen von "\u2010" durch "-" in JSON Datei
    def replace_unicode_hyphen(self, data):
        if isinstance(data, str):
            return data.replace("\\u2010", "-")
        elif isinstance(data, list):
            return [self.replace_unicode_hyphen(item) for item in data]
        elif isinstance(data, dict):
            return {key: self.replace_unicode_hyphen(value) for key, value in data.items()}
        else:
            return data

    # Speichert die extrahierten Substanzdaten in einer JSON-Datei.
    def save_to_json(self, filename):
        updated_substances = self.replace_unicode_hyphen(self.substances)
        with open(filename, 'w', encoding='utf-8') as json_file:
            json.dump(updated_substances, json_file, indent=4)

    # Führt den gesamten Extraktionsprozess aus.
    def run(self):
        html = self.fetch_data()
        self.parse_html(html)
        self.save_to_json('substancesV4.json')

if __name__ == "__main__":
    url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
    extractor = SubstanceExtractor(url)
    extractor.run()
    print("Datenextraktion abgeschlossen. Die Ergebnisse sind in substancesV4.json gespeichert.")
