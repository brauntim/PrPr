
# PrPr

This repository contains a Python script designed to scrape data from the [Policija website](https://www.policija.si/) and convert chemical names to SMILES and InChI using the [OPSIN API](https://opsin.ch.cam.ac.uk/).

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [File Usage Order](#file-usage-order)
- [Contributing](#contributing)

## Installation

### Windows

1. Open your terminal (cmd).
2. Install the required packages:
    ```sh
    pip install requests
    pip install bs4
    pip install rdkit
    ```
3. Alternatively, for RDKit:
    ```sh
    conda install -c conda-forge rdkit
    ```

### Linux and MacOS

1. Open your terminal.
2. Install the required packages:
    ```sh
    pip install requests bs4 rdkit
    ```
3. Alternatively, for RDKit:
    ```sh
    conda install -c conda-forge rdkit
    ```

## Usage

To run the scraper, execute the following command in your terminal:
```sh
python main.py
```
This will start the data extraction process and save the results to `Tim_Jonas_Policija.json`.

To merge data from JSON files without duplicates, execute the following command:
```sh
python jsons/forging.py
```
This will compare data from the generated JSON files and create a merged file without duplicates.

To use the search engine functionalities within the project, execute the following command:
```sh
python searchengine.py
```
This script provides search engine related functionalities for the project.

## Project Structure

- **main.py**: The main script to run the web scraper.
- **logs/**: Directory containing log files.
- **jsons/**: Directory containing JSON files and `forging.py`.
    - **forging.py**: Compares data from two JSON files and merges them into one file without duplicates.
- **searchengine.py**: Python script related to the search engine functionality.
- **Programmierprojekt_SS2024_Schildgen.pdf**: Project documentation in PDF format.
- **tempCodeRunnerFile.py**: Temporary code runner file.
- - **__pycache__/**: Directory containing cached bytecode files.
- **.idea/**: Directory containing project-specific configuration files for IntelliJ IDEA/PyCharm.

## File Usage Order

1. **main.py**: Start with this script to scrape data from the Policija website. It will generate the initial JSON files.
2. **forging.py**: Use this script after running `main.py` to merge data from the generated JSON files into one comprehensive file without duplicates.
3. **searchengine.py**: This script can be used as needed for search engine related functionalities within the project.

## Contributing

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Commit your changes (`git commit -am 'Add new feature'`).
4. Push to the branch (`git push origin feature-branch`).
5. Create a new Pull Request.

## Authors

Created by Tim Braun and Jonas Holzapfel.
