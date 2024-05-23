import json

def filter_smiles():
    input_smiles = input("Smiles eingeben: ")
    contains = False

    for substance in data:
        if substance["smiles"] == input_smiles:
            print("Eintrag mit dieser Smiles gefunden")
            print(substance)
            print("\n")
            contains = True

    if contains == False:
        print("Smiles \"" + input_smiles + "\" nicht in Datei enthalten.")

def filter_formular():
    input_formular = input("Summenformel eingeben: ")
    contains = False

    for substance in data:
        if substance["formular"] == input_formular:
            print("Eintrag mit dieser Summenformel gefunden")
            print(substance)
            print("\n")
            contains = True

    if contains == False:
        print("Formel \"" + input_formular + "\" nicht in Datei enthalten.")

def filter_mass():
    input_mass_min = input("minimale Masse eingeben: ")
    input_mass_max = input("maximale Masse eingeben: ")

    for substance in data:
        if input_mass_max >= substance["molecular_mass"] >= input_mass_min:
            print("Eintrag mit dieser Masse:")
            print(substance)
            print("\n")

def searchStart():
    while True:
        print("Alle Daten ausgeben: 1 \nNach Smiles filtern: 2 \nNach Summenformel filtern: 3 \nNach Masse filtern: 4\n"
              "Beenden: 5\n")

        operation = input("Eingabe: ")
        operation = int(operation)

        if operation == 1:
            for info in data:
                print(info)

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


with open("substances.json") as file:
    data = json.load(file)

    print("Willkommen zu dieser Suchmaschine für Designerdrogen")
    print("______________________________________________________")

end = 0

while end != 5:
    end = searchStart()

print("Suchmaschine beendet")






# for name in data:
#  print(name['names'])
#   print("")
#    print(name['category'])
