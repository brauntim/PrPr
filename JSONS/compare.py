import json
import os

# Einlesen einer JSON-Datei
def read_json(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        data = json.load(file)
    return data

# Schreiben einer JSON-Datei
def write_json(file_path, data):
    with open(file_path, 'w', encoding='utf-8') as file:
        json.dump(data, file, indent=2)

# Vergleich der InChI-Keys und Logging
def compare_and_log(inchi_key_a, inchi_key_b, log_file):
    if inchi_key_a == inchi_key_b:
        log_message = f"Der InChI-Key '{inchi_key_a}' von Datei A ist ebenfalls in Datei B vorhanden."
        print(log_message)
        with open(log_file, 'a', encoding='utf-8') as log:
            log.write(log_message + '\n')

# Liste der JSON-Dateien im Ordner 'jsons'
def list_json_files():
    json_files = [f for f in os.listdir('jsons') if f.endswith('.json')]
    if not json_files:
        print("Keine JSON-Dateien im JSONS-Ordner gefunden.")
        return []
    return json_files

# Benutzer zur Auswahl einer JSON-Datei auffordern
def get_filename_from_user(prompt):
    json_files = list_json_files()
    if not json_files:
        return None
    print(prompt)
    for idx, file in enumerate(json_files, start=1):
        print(f"{idx}: {file}")
    while True:
        try:
            file_index = int(input("Bitte geben Sie die Nummer der Datei ein: ")) - 1
            if 0 <= file_index < len(json_files):
                return json_files[file_index]
            else:
                print("Ungültige Nummer. Bitte erneut versuchen.")
        except ValueError:
            print("Ungültige Eingabe. Bitte eine Nummer eingeben.")

# Auswahl der JSON-Dateien für Datei A und Datei B
file_a = get_filename_from_user("Bitte wählen Sie die JSON-Datei für Datei A aus:")
if not file_a:
    exit("Keine Datei für Datei A ausgewählt. Beenden.")
file_b = get_filename_from_user("Bitte wählen Sie die JSON-Datei für Datei B aus:")
if not file_b:
    exit("Keine Datei für Datei B ausgewählt. Beenden.")
log_file = 'log.txt'

# Daten aus den Dateien lesen
try:
    data_a = read_json(f'jsons/{file_a}')
    data_b = read_json(f'jsons/{file_b}')
except FileNotFoundError:
    print("Eine der angegebenen Dateien wurde nicht gefunden. Stellen Sie sicher, dass die Dateinamen korrekt sind.")
    exit()

# Log-Datei vorbereiten
with open(log_file, 'w', encoding='utf-8') as log:
    log.write("Vergleich der InChI-Keys zwischen Datei A und Datei B:\n\n")

# Durch die Daten in Datei A iterieren und InChI-Keys vergleichen
for entry_a in data_a:
    inchi_key_a = entry_a.get('inchi_key')

    # Durch die Daten in Datei B iterieren und vergleichen
    for entry_b in data_b:
        inchi_key_b = entry_b.get('inchi_key')
        compare_and_log(inchi_key_a, inchi_key_b, log_file)

# Durch die Daten in Datei B iterieren und Daten an Datei A anhängen, wenn nicht vorhanden
for entry_b in data_b:
    inchi_key_b = entry_b.get('inchi_key')
    found_in_a = False

    # Überprüfen, ob InChI-Key bereits in Datei A vorhanden ist
    for entry_a in data_a:
        inchi_key_a = entry_a.get('inchi_key')
        if inchi_key_a == inchi_key_b:
            found_in_a = True
            break

    # Wenn InChI-Key nicht in Datei A gefunden wurde, an Datei A anhängen
    if not found_in_a:
        data_a.append(entry_b)
        log_message = f"Der InChI-Key '{inchi_key_b}' von Datei B wurde zu Datei A hinzugefügt."
        print(log_message)
        with open(log_file, 'a', encoding='utf-8') as log:
            log.write(log_message + '\n')

# Aktualisierte Daten in Datei A schreiben
write_json(f'jsons/{file_a}', data_a)

# Ausgabe für den Benutzer
print(f"Vergleich abgeschlossen. Aktualisierte Daten wurden in 'jsons/{file_a}' geschrieben.")
print(f"Das Log der Vergleiche und Ergänzungen wurde in '{log_file}' gespeichert.")
