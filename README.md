# TCGA,DrugFinder

Github: [GitHub - akcaymd/TCGA-DrugFinder](https://github.com/akcaymd/TCGA-DrugFinder)

This project aims to develop a bioinformatics analysis tool that integrates cancer clinical data with drug-level molecular information to evaluate treatment outcomes.

The program generates summary statistics for each drug MeSH term, including patient counts, treatment categories, smoking distribution, and mean survival ranked in decreasing order. The system supports side-by-side statistical comparison of selected compounds and is designed to function through both command-line and graphical user interfaces, enabling reproducible and systematic analysis of cancer drug response patterns.

## How to run?

1. Create a python virtual environment inside project folder.

```bash
python -m venv .venv
```

2. Activate created environment.

Linux / macOS

```bash
source venv/bin/activate
```

Windows (PowerShell)

```bash
venv\Scripts\Activate.ps1
```

Windows (CMD)

```bash
venv\Scripts\activate.bat
```

3. Install required packages.

```bash
pip install -r requirements.txt
```

4. To run app on cli:

```
python cgi-bin/app.py
```

4. To run app on web run this command and go to `http://localhost:8000/cgi-bin/home.cgi

```
python3 -m http.server --cgi 8000
```

## References and Tools Used

- [GDC Application Programming Interface (API) \| NCI Genomic Data Commons](http://gdc.cancer.gov/developers/gdc-application-programming-interface-api)
- [Entrez Molecular Sequence Database System](https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html)
