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

formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%d-%m-%Y %H:%M:%S")

file_logger = logging.FileHandler('logs/scraper.log')
file_logger.setLevel(logging.DEBUG)
file_logger.setFormatter(formatter)

console_logger = logging.StreamHandler()
console_logger.setLevel(logging.INFO)
console_logger.setFormatter(formatter)

logger.addHandler(file_logger)
logger.addHandler(console_logger)


class SubstanceExtractor:
    def __init__(self, url):
        self.url = url
        self.substances = []

    def fetch_data(self):
        logger.info("Scraping started")
        response = requests.get(self.url)
        response.raise_for_status()
        return response.text

    def convert_to_canonical_smiles(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles

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

    def parse_row(self, row):
        columns = row.find_all('td')
        if len(columns) < 17:
            return None

        category = columns[0].text.strip()
        name = columns[1].text.strip()
        iupac_name = columns[3].text.strip()
        formula = columns[5].text.strip()
        molecular_mass = columns[6].text.strip()
        std_inchi_key = columns[12].text.strip()
        date_of_entry = columns[15].text.strip()
        report_updated = columns[16].text.strip()

        last_changed_at = report_updated if report_updated else date_of_entry
        deleted = '<strike>' in name

        try:
            smiles, std_inchi = self.fetch_opsin_data(iupac_name)
        except Exception as e:
            logger.error(f"Error fetching data for IUPAC name {iupac_name}: {e}")
            smiles, std_inchi = "", ""

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
                "url": self.url
            },
            "validated": None,
            "deleted": deleted,
            "last_modified": last_changed_at,
            "details": {}
        }
        return substance_data

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
            if 'iupac_names' in data:
                data['iupac_names'] = self.replace_unicode_characters(data['iupac_names'])
            return {key: self.replace_unicode_characters(value) for key, value in data.items()}
        else:
            return data

    def validate_data(self):
        for substance in self.substances:
            smiles = substance.get("smiles", "")
            if smiles:
                molecule = Chem.MolFromSmiles(smiles)
                if molecule:
                    calculated_molecular_mass = Descriptors.MolWt(molecule)
                    calculated_formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)
                    try:
                        molecular_mass_valid = abs(
                            float(substance["molecular_mass"]) - calculated_molecular_mass) < 0.99
                    except ValueError:
                        molecular_mass_valid = False

                    formula_valid = substance["formula"] == calculated_formula
                    substance["validated"] = molecular_mass_valid and formula_valid
                else:
                    substance["validated"] = None
            else:
                substance["validated"] = None

    def save_to_json(self, filename):
        updated_substances = self.replace_unicode_characters(self.substances)
        with open(f"jsons/{filename}", 'w', encoding='utf-8') as json_file:
            json.dump(updated_substances, json_file, indent=4)
        logger.info(f"Saving substances to {filename}")
        logger.info("Scraping finished")

    def run(self, filename):
        html = self.fetch_data()
        self.parse_html(html)
        self.validate_data()
        self.save_to_json(filename)


def start_scraping(filename):
    url = "https://www.policija.si/apps/nfl_response_web/seznam.php"
    extractor = SubstanceExtractor(url)
    extractor.run(filename)
    print(f"Datenextraktion abgeschlossen. Die Ergebnisse sind in jsons/{filename} gespeichert.")


if __name__ == "__main__":
    start_scraping("Tim_Jonas_Policija.json")