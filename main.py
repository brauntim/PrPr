#
#   webscraper for https://www.policija.si/apps/nfl_response_web/seznam.php
#

import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import inchi
import time

#    Diese Klasse handhabt das Abrufen, Parsen und Speichern von Substanzdaten von einer gegebenen URL.
class SubstanceExtractor:

#   Initialisiert eine neue Instanz der SubstanceExtractor-Klasse.
    def __init__(self, url):
        self.url = url
        self.substances = []

#   GET-Anfrage an die Website und gibt HTML-Inhalt zurück
    def fetch_data(self):
        response = requests.get(self.url)
        response.raise_for_status()
        return response.text
    
    def get_inchi_from_inchikey(self, inchi_key):
        inchi_url = f'https://cactus.nci.nih.gov/chemical/structure/{inchi_key}/inchi'
        response = requests.get(inchi_url)
        if response.status_code == 200:
            return response.text.strip()
        return None

    def convert_inchi_to_smiles(self, inchi):
        mol = Chem.MolFromInchi(inchi)
        if mol:
            return Chem.MolToSmiles(mol)
        return None


#   Parst den HTML-Inhalt und extrahiert die Substanzdaten
    def parse_html(self, html):
        soup = BeautifulSoup(html, 'html.parser')
        for row in soup.select('table tr'):
            columns = row.find_all('td')
            if len(columns) < 17:  # Sicherstellen, dass die Zeile genügend Spalten hat
                continue

            # Extrahieren der Daten aus den entsprechenden Spalten
            category = columns[0].text.strip()
            name = columns[1].text.strip()
            # smiles = columns[2].text.strip()
            iupac_name = columns[3].text.strip()
            formular = columns[5].text.strip()
            molecular_mass = columns[6].text.strip()
            inchi_key = columns[12].text.strip()
            date_of_entry = columns[15].text.strip()
            report_updated = columns[16].text.strip()
            # details =

           # inchi = self.get_inchi_from_inchikey(inchi_key)
           # smiles = self.convert_inchi_to_smiles(inchi) 

            # Bestimmen des Datums der letzten Änderung
            if report_updated:
                last_changed_at = report_updated
            else:
                last_changed_at = date_of_entry
                
            # Überprüfen, ob die Substanz gelöscht wurde
            deleted = '<strike>' in name

            # Erstellen eines Dictionaries mit den extrahierten Daten
            substance_data = {
                "smiles": "",
                "names": name,
                "iupac_name": iupac_name,
                "formular": formular,
                "molecular_mass": molecular_mass,
                "Inchi": "",  
                "InchiKey": inchi_key,
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

            # Hinzufügen der Substanzdaten zur Liste
            self.substances.append(substance_data)


#   Speichert die extrahierten Substanzdaten in einer JSON-Datei.
    def save_to_json(self, filename):
        with open(filename, 'w') as json_file:
            json.dump(self.substances, json_file, indent=4)

#   Führt den gesamten Extraktionsprozess aus.
    def run(self):
        html = self.fetch_data()
        self.parse_html(html)
        self.save_to_json('substances.json')

    
if __name__ == "__main__":
    url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
    extractor = SubstanceExtractor(url)
    extractor.run()
    print("Datenextraktion abgeschlossen. Die Ergebnisse sind in substances.json gespeichert.")
