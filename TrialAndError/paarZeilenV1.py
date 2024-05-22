import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import inchi
import time

class SubstanceExtractor:

    def __init__(self, url):
        self.url = url
        self.substances = []

    def fetch_data(self):
        response = requests.get(self.url)
        response.raise_for_status()
        return response.text
    
    def parse_html(self, html):
        soup = BeautifulSoup(html, 'html.parser')
        rows = soup.select('table tr')
        for index, row in enumerate(rows):
            if index >= 10:  # Nur die ersten 10 Zeilen betrachten
                break

            columns = row.find_all('td')
            if len(columns) < 17:
                continue

            category = columns[0].text.strip()
            name = columns[1].text.strip()
            iupac_name = columns[3].text.strip()
            formular = columns[5].text.strip()
            molecular_mass = columns[6].text.strip()
            inchi_key = columns[12].text.strip()
            date_of_entry = columns[15].text.strip()
            report_updated = columns[16].text.strip()

            
            inchi_url = f'https://cactus.nci.nih.gov/chemical/structure/{inchi_key}/inchi'
            response = requests.get(inchi_url)
            inchi = response.text.strip()


            mol = Chem.MolFromInchi(inchi)
            smiles = Chem.MolToSmiles(mol)
        


            if report_updated:
                last_changed_at = report_updated
            else:
                last_changed_at = date_of_entry
                
            deleted = '<strike>' in name



            substance_data = {
                "smiles": smiles,
                "names": name,
                "iupac_name": iupac_name,
                "formular": formular,
                "molecular_mass": molecular_mass,
                "Inchi": inchi,  
                "InchiKey": inchi_key,
                "cas_num": "",
                "category": category,
                "source_name": "Policija",
                "source_url": self.url,
                "valid": None,
                "deleted": deleted,
                "last_changed_at": last_changed_at,
                "version": "1.0",
                "details": {}
            }

            self.substances.append(substance_data)



    def save_to_json(self, filename):
        with open(filename, 'w') as json_file:
            json.dump(self.substances, json_file, indent=4)

    def run(self):
        html = self.fetch_data()
        self.parse_html(html)
        self.save_to_json('paarZeilen2.json')

if __name__ == "__main__":
    url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
    extractor = SubstanceExtractor(url)
    extractor.run()
    print("Datenextraktion abgeschlossen. Die Ergebnisse sind in substances.json gespeichert.")
