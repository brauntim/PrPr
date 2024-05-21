import requests
from bs4 import BeautifulSoup
import os

# URL of the monographs page
url = 'https://www.cfsre.org/nps-discovery/monographs'

# Function to get the HTML content of the page
def get_html(url):
    response = requests.get(url)
    response.raise_for_status()
    return response.text

# Function to parse the HTML and extract PDF links
def extract_pdfs(html):
    soup = BeautifulSoup(html, 'html.parser')
    pdf_links = []
    # Assuming PDFs are linked with <a> tags and contain 'pdf' in href
    for link in soup.find_all('a', href=True):
        if 'pdf' in link['href']:
            pdf_links.append(link['href'])
    return pdf_links

# Function to download PDFs
def download_pdfs(pdf_links, download_folder='pdfs'):
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)
    for link in pdf_links:
        filename = os.path.join(download_folder, os.path.basename(link))
        response = requests.get(link)
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f'Downloaded: {filename}')

# Main script
if __name__ == '__main__':
    html_content = get_html(url)
    pdf_links = extract_pdfs(html_content)
    download_pdfs(pdf_links)
