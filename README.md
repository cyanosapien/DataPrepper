# DataPrepper
Shiny app for file format conversion and taxonomy correction in OTU counts tables (16S or ITS). 

## ✨ Features
- Accepts an OTU count table for either 16S or ITS data.
- Cleans and reformats taxonomy strings, removes low-count OTUs, and excludes chloroplast/mitochondrial entries.
- Optionally corrects genus-level taxonomy using NCBI Entrez.
- Generates a metadata table with inferred treatments and replicate numbers.
- Lets users edit metadata, rename samples, and apply bulk treatment/replicate changes.
- Provides downloads for the processed OTU table, metadata, and taxonomy-correction log.
- Includes a second workflow that converts a “% of Library” OTU file into a transposed format suitable for PAST.
- Displays previews and processing status directly in the app.


## 📦 Installation

Access on the Web
- https://somera-lab.shinyapps.io/DataPrepper/

Run Locally 
- Download the app-DataPrepper_v2.1.R file
- Open it in R studio and click "Run App"

## 🚀 Usage

- Click on the “Browse…” button under “Upload OTU.txt File” and select one of the OTU counts files to process. You can also find the file in your file navigator and then drag and drop it into the web app over the “Browse” button.
  - For bacterial analysis, open the appropriate zip file from the table above, then navigate through the folders: analysisfiles>prokaryote>bacteria> and select the file ending in “OTU.txt”.
  - For fungal analysis, open the appropriate zip file from the table above, then navigate through the folders: analysisfiles>eukaryote>fungi> and select the file ending in “OTU.txt”.
- Select the appropriate Data Type button: 16S (bacterial) or ITS (fungal).
- Select with or without taxonomy correction (via the NCBI Entrez database).
- Click the button “Process OTU File Data”. You can preview the data on the right side of the webpage. Ensure that all your treatment names are included and correct, and that the taxa information is correctly formatted under the header “O_T_U”. Here is an example: 
bacteria/firmicutes/bacilli/bacillales/bacillaceae/bacillus/bacillus circulans
- Use the download buttons to receive your re-formatted OTU counts file, metadata file, and taxonomy correction log (if taxonomy correction were selected). 

## 🤝 Contributing

## 📄 License
Distributed under the MIT License. See `LICENSE` for more information.

## ✉️ Contact
Dr. Martinez, ermartinez17@gmail.com
