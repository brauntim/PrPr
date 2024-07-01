import json

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

# Eingabe der Dateinamen von Datei A und Datei B durch den Benutzer
file_a = input("Bitte geben Sie den Dateinamen für Datei A (JSON) ein: ")
file_b = input("Bitte geben Sie den Dateinamen für Datei B (JSON) ein: ")
log_file = 'log.txt'

# Daten aus den Dateien lesen
try:
    data_a = read_json(file_a)
    data_b = read_json(file_b)
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
write_json(file_a, data_a)

# Ausgabe für den Benutzer
print(f"Vergleich abgeschlossen. Aktualisierte Daten wurden in '{file_a}' geschrieben.")
