import requests
from bs4 import BeautifulSoup
import json
from rdkit import Chem
from rdkit.Chem import MolToSmiles
from datetime import datetime
import pytz

BASE_URL = 'https://www.policija.si/apps/nfl_response_web/seznam.php'
SUBSTANCE_BASE_URL = 'https://www.policija.si/apps/nfl_response_web/'

def fetch_substances():
    response = requests.get(BASE_URL)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')
    table = soup.find('table', {'id': 'seznam'})
    if not table:
        raise ValueError("Table with id 'seznam' not found")
    rows = table.find_all('tr')[1:]  # Skip the header row
    substances = []
    for row in rows:
        cols = row.find_all('td')
        if len(cols) >= 2:
            substance_name = cols[0].text.strip()
            categories = [cat.strip() for cat in cols[1].text.split(',')]
            substance_link = cols[0].find('a')
            if substance_link and 'href' in substance_link.attrs:
                substance_url = SUBSTANCE_BASE_URL + substance_link['href']
                substances.append({
                    'name': substance_name,
                    'categories': categories,
                    'url': substance_url
                })
            else:
                substances.append({
                    'name': substance_name,
                    'categories': categories,
                    'url': None
                })
        else:
            print(f"Skipping row due to insufficient columns: {row}")
    return substances

def fetch_substance_details(substance):
    if not substance['url']:
        return {}
    response = requests.get(substance['url'])
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')
    details = {}
    for row in soup.find_all('tr'):
        cells = row.find_all('td')
        if len(cells) == 2:
            key = cells[0].text.strip()
            value = cells[1].text.strip()
            details[key] = value
    return details

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return MolToSmiles(mol)
    return None

def fetch_inchi_from_inchikey(inchikey):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{inchikey}/stdinchi"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

def transform_to_json(substances):
    data = []
    for substance in substances:
        details = fetch_substance_details(substance)
        canonical_smiles = canonicalize_smiles(details.get('SMILES', ''))
        inchi = fetch_inchi_from_inchikey(details.get('StdInChIKey', ''))
        
        berlin_tz = pytz.timezone('Europe/Berlin')
        last_changed_at = details.get('report updated', None)
        if not last_changed_at:
            last_changed_at = details.get('entry of report', None)
        if last_changed_at:
            last_changed_at = datetime.strptime(last_changed_at, '%d.%m.%Y').astimezone(berlin_tz).isoformat()
        else:
            last_changed_at = datetime.now(berlin_tz).isoformat()

        entry = {
            "smiles": canonical_smiles,
            "names": [substance['name']] + details.get('Synonyme', '').split(','),
            "iupac_name": details.get('Formeller Name', ''),
            "formular": details.get('Molekulare Summenformel', ''),
            "inchi": inchi,
            "inchi_key": details.get('StdInChIKey', ''),
            "molecular_mass": float(details.get('Mw (g/mol) per base form NPS1', 0.0)),
            "cas_num": details.get('CAS-Nummer', ''),
            "category": substance['categories'],
            "source_name": "Policija",
            "source_url": BASE_URL,
            "valid": None,
            "deleted": False,
            "last_changed_at": last_changed_at,
            "version": 1.0,
            "details": details
        }
        data.append(entry)
    return data

def main():
    substances = fetch_substances()
    data = transform_to_json(substances)
    with open('substances.json', 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()
