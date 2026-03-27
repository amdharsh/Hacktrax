# BioDetect Backend — Setup & Run Guide

## Requirements
- Python 3.10+

## Install
```bash
cd biodetect_backend
pip install -r requirements.txt
```

## Run
```bash
uvicorn main:app --reload --port 8000
```

## API Endpoints

### GET /diseases
Returns all 8 diseases in the ICD-10 database with their biomarkers.

### GET /biomarkers
Returns all unique biomarkers and which diseases they apply to.

### POST /analyse
Main endpoint. Send patient info + biomarker values, get ranked disease scores.

**Example request body:**
```json
{
  "patient_name": "Jane Doe",
  "patient_id": "PT-001",
  "age": 52,
  "sex": "Female",
  "biomarkers": {
    "HbA1c": 7.8,
    "FPG": 145,
    "CA-125": 80,
    "Troponin I": 0.02,
    "BNP": 60
  }
}
```

**Example response:**
```json
{
  "patient": { "name": "Jane Doe", ... },
  "primary_diagnosis": {
    "disease": "Type 2 Diabetes Mellitus",
    "icd10": "E11",
    "confidence": 76,
    "risk": "high",
    "specialist": "Endocrinologist",
    "recommendation": "...",
    "markers": { ... }
  },
  "all_results": [ ... ],
  "database_version": "ICD-10-2024",
  "total_diseases_screened": 8
}
```

## ICD-10 Database Coverage

| Disease | ICD-10 | Key Biomarkers |
|---|---|---|
| Ovarian Cancer | C56 | CA-125, HE4, MUC-16 |
| Type 2 Diabetes | E11 | HbA1c, FPG, C-Peptide |
| Alzheimer's Disease | G30 | Amyloid-b42, p-tau 181, NFL |
| Cardiovascular Disease | I25 | Troponin I, BNP, CRP, LDL |
| Chronic Kidney Disease | N18 | Creatinine, eGFR, Urine Albumin |
| Liver Disease | K74 | ALT, AST, Bilirubin, Albumin |
| Hypothyroidism | E03 | TSH, Free T4, Free T3 |
| Anemia (Iron Deficiency) | D50 | Hemoglobin, Serum Ferritin, TIBC |

## Connecting to BioDetect Frontend
Change the frontend's fetch call from local logic to:
```
POST http://localhost:8000/analyse
```
