"""
BioDetect FastAPI Backend
ICD-10 aligned biomarker-to-disease mapping engine
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional

app = FastAPI(title="BioDetect API", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ─────────────────────────────────────────────
# ICD-10 ALIGNED DISEASE DATABASE
# ─────────────────────────────────────────────
DISEASE_DB = {
    "Ovarian Cancer": {
        "icd10": "C56",
        "icd10_label": "Malignant neoplasm of ovary",
        "biomarkers": {
            "CA-125": {
                "unit": "U/mL", "reference": "<35", "inverted": False,
                "thresholds": [
                    {"min": 35,  "max": 70,  "level": "warn",     "score": 0.45},
                    {"min": 70,  "max": 200, "level": "high",     "score": 0.70},
                    {"min": 200, "max": None,"level": "critical", "score": 0.92},
                ],
                "normal_score": 0.05, "weight": 0.45,
            },
            "HE4": {
                "unit": "pmol/L", "reference": "<70", "inverted": False,
                "thresholds": [
                    {"min": 70,  "max": 90,  "level": "warn",     "score": 0.40},
                    {"min": 90,  "max": 150, "level": "high",     "score": 0.70},
                    {"min": 150, "max": None,"level": "critical", "score": 0.90},
                ],
                "normal_score": 0.05, "weight": 0.35,
            },
            "MUC-16": {
                "unit": "U/mL", "reference": "<35", "inverted": False,
                "thresholds": [
                    {"min": 35,  "max": 65,  "level": "warn", "score": 0.40},
                    {"min": 65,  "max": None,"level": "high", "score": 0.72},
                ],
                "normal_score": 0.04, "weight": 0.20,
            },
        },
        "category": "Oncology", "specialist": "Oncologist / Gynecologist", "urgency": "urgent",
    },

    "Type 2 Diabetes Mellitus": {
        "icd10": "E11",
        "icd10_label": "Type 2 diabetes mellitus",
        "biomarkers": {
            "HbA1c": {
                "unit": "%", "reference": "<5.7", "inverted": False,
                "thresholds": [
                    {"min": 5.7, "max": 6.5, "level": "warn",     "score": 0.50},
                    {"min": 6.5, "max": 8.0, "level": "high",     "score": 0.80},
                    {"min": 8.0, "max": None,"level": "critical", "score": 0.95},
                ],
                "normal_score": 0.05, "weight": 0.55,
            },
            "FPG": {
                "unit": "mg/dL", "reference": "<100", "inverted": False,
                "thresholds": [
                    {"min": 100, "max": 126, "level": "warn",     "score": 0.45},
                    {"min": 126, "max": 200, "level": "high",     "score": 0.75},
                    {"min": 200, "max": None,"level": "critical", "score": 0.92},
                ],
                "normal_score": 0.05, "weight": 0.35,
            },
            "C-Peptide": {
                "unit": "ng/mL", "reference": "0.8-3.1", "inverted": False,
                "thresholds": [
                    {"min": 3.1, "max": 5.0, "level": "warn", "score": 0.45},
                    {"min": 5.0, "max": None,"level": "high", "score": 0.70},
                ],
                "normal_score": 0.05, "weight": 0.10,
            },
        },
        "category": "Endocrinology", "specialist": "Endocrinologist", "urgency": "moderate",
    },

    "Alzheimer's Disease": {
        "icd10": "G30",
        "icd10_label": "Alzheimer's disease",
        "biomarkers": {
            "Amyloid-b42": {
                "unit": "pg/mL", "reference": ">800", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 400, "level": "critical", "score": 0.88},
                    {"min": 400, "max": 600, "level": "high",     "score": 0.65},
                    {"min": 600, "max": 800, "level": "warn",     "score": 0.35},
                ],
                "normal_score": 0.05, "weight": 0.45,
            },
            "p-tau 181": {
                "unit": "pg/mL", "reference": "<19", "inverted": False,
                "thresholds": [
                    {"min": 19, "max": 25, "level": "warn",     "score": 0.45},
                    {"min": 25, "max": 40, "level": "high",     "score": 0.72},
                    {"min": 40, "max": None,"level": "critical","score": 0.92},
                ],
                "normal_score": 0.05, "weight": 0.40,
            },
            "NFL": {
                "unit": "pg/mL", "reference": "<10", "inverted": False,
                "thresholds": [
                    {"min": 10, "max": 20, "level": "warn", "score": 0.40},
                    {"min": 20, "max": None,"level": "high","score": 0.70},
                ],
                "normal_score": 0.04, "weight": 0.15,
            },
        },
        "category": "Neurology", "specialist": "Neurologist", "urgency": "moderate",
    },

    "Cardiovascular Disease": {
        "icd10": "I25",
        "icd10_label": "Chronic ischaemic heart disease",
        "biomarkers": {
            "Troponin I": {
                "unit": "ng/mL", "reference": "<0.04", "inverted": False,
                "thresholds": [
                    {"min": 0.01, "max": 0.04, "level": "warn",     "score": 0.40},
                    {"min": 0.04, "max": 0.10, "level": "high",     "score": 0.72},
                    {"min": 0.10, "max": None, "level": "critical", "score": 0.95},
                ],
                "normal_score": 0.05, "weight": 0.40,
            },
            "BNP": {
                "unit": "pg/mL", "reference": "<100", "inverted": False,
                "thresholds": [
                    {"min": 100, "max": 400, "level": "warn",     "score": 0.45},
                    {"min": 400, "max": 900, "level": "high",     "score": 0.75},
                    {"min": 900, "max": None,"level": "critical", "score": 0.93},
                ],
                "normal_score": 0.05, "weight": 0.35,
            },
            "CRP": {
                "unit": "mg/L", "reference": "<3.0", "inverted": False,
                "thresholds": [
                    {"min": 3.0, "max": 10.0, "level": "warn", "score": 0.35},
                    {"min": 10.0,"max": None, "level": "high", "score": 0.65},
                ],
                "normal_score": 0.04, "weight": 0.15,
            },
            "LDL": {
                "unit": "mg/dL", "reference": "<100", "inverted": False,
                "thresholds": [
                    {"min": 100, "max": 130, "level": "warn",     "score": 0.30},
                    {"min": 130, "max": 160, "level": "high",     "score": 0.55},
                    {"min": 160, "max": None,"level": "critical", "score": 0.78},
                ],
                "normal_score": 0.04, "weight": 0.10,
            },
        },
        "category": "Cardiology", "specialist": "Cardiologist", "urgency": "urgent",
    },

    "Chronic Kidney Disease": {
        "icd10": "N18",
        "icd10_label": "Chronic kidney disease",
        "biomarkers": {
            "Creatinine": {
                "unit": "mg/dL", "reference": "0.6-1.2", "inverted": False,
                "thresholds": [
                    {"min": 1.2, "max": 2.0, "level": "warn",     "score": 0.45},
                    {"min": 2.0, "max": 4.0, "level": "high",     "score": 0.72},
                    {"min": 4.0, "max": None,"level": "critical", "score": 0.92},
                ],
                "normal_score": 0.05, "weight": 0.40,
            },
            "eGFR": {
                "unit": "mL/min/1.73m2", "reference": ">60", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 30, "level": "critical", "score": 0.90},
                    {"min": 30,  "max": 45, "level": "high",     "score": 0.70},
                    {"min": 45,  "max": 60, "level": "warn",     "score": 0.45},
                ],
                "normal_score": 0.05, "weight": 0.45,
            },
            "Urine Albumin": {
                "unit": "mg/g", "reference": "<30", "inverted": False,
                "thresholds": [
                    {"min": 30,  "max": 300, "level": "warn", "score": 0.40},
                    {"min": 300, "max": None,"level": "high", "score": 0.72},
                ],
                "normal_score": 0.03, "weight": 0.15,
            },
        },
        "category": "Nephrology", "specialist": "Nephrologist", "urgency": "moderate",
    },

    "Liver Disease": {
        "icd10": "K74",
        "icd10_label": "Fibrosis and cirrhosis of liver",
        "biomarkers": {
            "ALT": {
                "unit": "U/L", "reference": "<40", "inverted": False,
                "thresholds": [
                    {"min": 40,  "max": 80,  "level": "warn",     "score": 0.40},
                    {"min": 80,  "max": 200, "level": "high",     "score": 0.68},
                    {"min": 200, "max": None,"level": "critical", "score": 0.90},
                ],
                "normal_score": 0.04, "weight": 0.30,
            },
            "AST": {
                "unit": "U/L", "reference": "<40", "inverted": False,
                "thresholds": [
                    {"min": 40,  "max": 80,  "level": "warn",     "score": 0.40},
                    {"min": 80,  "max": 200, "level": "high",     "score": 0.68},
                    {"min": 200, "max": None,"level": "critical", "score": 0.90},
                ],
                "normal_score": 0.04, "weight": 0.30,
            },
            "Bilirubin": {
                "unit": "mg/dL", "reference": "<1.2", "inverted": False,
                "thresholds": [
                    {"min": 1.2, "max": 3.0, "level": "warn", "score": 0.42},
                    {"min": 3.0, "max": None,"level": "high", "score": 0.75},
                ],
                "normal_score": 0.04, "weight": 0.25,
            },
            "Albumin": {
                "unit": "g/dL", "reference": "3.5-5.0", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 2.5, "level": "critical", "score": 0.88},
                    {"min": 2.5, "max": 3.5, "level": "high",     "score": 0.60},
                ],
                "normal_score": 0.03, "weight": 0.15,
            },
        },
        "category": "Hepatology", "specialist": "Hepatologist / Gastroenterologist", "urgency": "moderate",
    },

    "Hypothyroidism": {
        "icd10": "E03",
        "icd10_label": "Other hypothyroidism",
        "biomarkers": {
            "TSH": {
                "unit": "mIU/L", "reference": "0.4-4.0", "inverted": False,
                "thresholds": [
                    {"min": 4.0, "max": 10.0, "level": "warn", "score": 0.50},
                    {"min": 10.0,"max": None, "level": "high", "score": 0.85},
                ],
                "normal_score": 0.05, "weight": 0.55,
            },
            "Free T4": {
                "unit": "ng/dL", "reference": "0.8-1.8", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 0.5, "level": "high", "score": 0.75},
                    {"min": 0.5, "max": 0.8, "level": "warn", "score": 0.45},
                ],
                "normal_score": 0.04, "weight": 0.35,
            },
            "Free T3": {
                "unit": "pg/mL", "reference": "2.3-4.2", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 2.0, "level": "high", "score": 0.70},
                    {"min": 2.0, "max": 2.3, "level": "warn", "score": 0.40},
                ],
                "normal_score": 0.03, "weight": 0.10,
            },
        },
        "category": "Endocrinology", "specialist": "Endocrinologist", "urgency": "routine",
    },

    "Anemia (Iron Deficiency)": {
        "icd10": "D50",
        "icd10_label": "Iron deficiency anemia",
        "biomarkers": {
            "Hemoglobin": {
                "unit": "g/dL", "reference": ">12 (F) / >13.5 (M)", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 8.0,  "level": "critical", "score": 0.90},
                    {"min": 8.0, "max": 11.0, "level": "high",     "score": 0.70},
                    {"min": 11.0,"max": 12.0, "level": "warn",     "score": 0.40},
                ],
                "normal_score": 0.04, "weight": 0.40,
            },
            "Serum Ferritin": {
                "unit": "ng/mL", "reference": "12-150", "inverted": True,
                "thresholds": [
                    {"min": None,"max": 12, "level": "high", "score": 0.78},
                    {"min": 12,  "max": 30, "level": "warn", "score": 0.45},
                ],
                "normal_score": 0.04, "weight": 0.40,
            },
            "TIBC": {
                "unit": "ug/dL", "reference": "250-370", "inverted": False,
                "thresholds": [
                    {"min": 370, "max": 450, "level": "warn", "score": 0.40},
                    {"min": 450, "max": None,"level": "high", "score": 0.70},
                ],
                "normal_score": 0.04, "weight": 0.20,
            },
        },
        "category": "Hematology", "specialist": "Hematologist", "urgency": "moderate",
    },
}


# ─────────────────────────────────────────────
# SCORING ENGINE
# ─────────────────────────────────────────────

def score_biomarker(value: float, marker_def: dict):
    inverted = marker_def.get("inverted", False)
    thresholds = marker_def["thresholds"]

    for threshold in thresholds:
        lo = threshold["min"]
        hi = threshold["max"]
        if inverted:
            in_range = (lo is None or value >= lo) and (hi is None or value < hi)
        else:
            in_range = (lo is None or value >= lo) and (hi is None or value < hi)
        if in_range:
            return threshold["score"], threshold["level"]

    return marker_def["normal_score"], "normal"


def compute_disease_score(submitted: dict, disease_def: dict):
    biomarker_defs = disease_def["biomarkers"]
    total_weight = 0
    weighted_score = 0
    marker_results = {}
    matched_count = 0

    for marker_name, marker_def in biomarker_defs.items():
        if marker_name in submitted:
            value = submitted[marker_name]
            score, level = score_biomarker(value, marker_def)
            weight = marker_def["weight"]
            weighted_score += score * weight
            total_weight += weight
            matched_count += 1
            marker_results[marker_name] = {
                "value": value,
                "unit": marker_def["unit"],
                "reference": marker_def["reference"],
                "score": round(score, 3),
                "level": level,
            }

    if total_weight == 0:
        return None

    final_score = weighted_score / total_weight
    coverage = matched_count / len(biomarker_defs)
    final_score = final_score * (0.6 + 0.4 * coverage)

    return {
        "score": min(round(final_score, 4), 0.99),
        "coverage": round(coverage, 2),
        "markers_matched": matched_count,
        "markers_total": len(biomarker_defs),
        "markers": marker_results,
    }


def get_risk_label(score: float) -> str:
    if score >= 0.80: return "critical"
    if score >= 0.65: return "high"
    if score >= 0.45: return "moderate"
    if score >= 0.25: return "low"
    return "minimal"


def get_recommendation(score: float, urgency: str) -> str:
    if score >= 0.80 or urgency == "urgent":
        return "Immediate specialist referral strongly recommended. Please seek medical attention within 24-48 hours."
    if score >= 0.60:
        return "Monitoring recommended. Schedule a follow-up with your physician within 1-2 weeks."
    if score >= 0.35:
        return "Borderline values detected. Repeat testing in 4-6 weeks and discuss with your doctor."
    return "Values within or near normal range. Continue routine annual health screening."


# ─────────────────────────────────────────────
# MODELS
# ─────────────────────────────────────────────

class AnalyseRequest(BaseModel):
    patient_name: str
    patient_id: Optional[str] = None
    age: Optional[int] = None
    sex: Optional[str] = None
    biomarkers: dict


# ─────────────────────────────────────────────
# ROUTES
# ─────────────────────────────────────────────

@app.get("/")
def root():
    return {"service": "BioDetect API", "version": "2.0", "status": "running"}


@app.get("/diseases")
def list_diseases():
    return {
        "count": len(DISEASE_DB),
        "diseases": [
            {
                "name": name,
                "icd10": d["icd10"],
                "icd10_label": d["icd10_label"],
                "category": d["category"],
                "biomarkers": list(d["biomarkers"].keys()),
            }
            for name, d in DISEASE_DB.items()
        ]
    }


@app.get("/biomarkers")
def list_biomarkers():
    seen = {}
    for disease_name, disease_def in DISEASE_DB.items():
        for bm_name, bm_def in disease_def["biomarkers"].items():
            if bm_name not in seen:
                seen[bm_name] = {"unit": bm_def["unit"], "reference": bm_def["reference"], "relevant_diseases": []}
            seen[bm_name]["relevant_diseases"].append(disease_name)
    return {"count": len(seen), "biomarkers": seen}


@app.post("/analyse")
def analyse(req: AnalyseRequest):
    if not req.biomarkers:
        raise HTTPException(status_code=400, detail="No biomarkers provided.")

    results = []
    for disease_name, disease_def in DISEASE_DB.items():
        scored = compute_disease_score(req.biomarkers, disease_def)
        if scored is None:
            continue

        confidence_pct = round(scored["score"] * 100)
        risk = get_risk_label(scored["score"])
        recommendation = get_recommendation(scored["score"], disease_def["urgency"])

        results.append({
            "disease": disease_name,
            "icd10": disease_def["icd10"],
            "icd10_label": disease_def["icd10_label"],
            "category": disease_def["category"],
            "specialist": disease_def["specialist"],
            "confidence": confidence_pct,
            "score": scored["score"],
            "risk": risk,
            "coverage": scored["coverage"],
            "markers_matched": scored["markers_matched"],
            "markers_total": scored["markers_total"],
            "markers": scored["markers"],
            "recommendation": recommendation,
        })

    if not results:
        raise HTTPException(status_code=422, detail="None of the submitted biomarkers matched any disease in the database.")

    results.sort(key=lambda x: x["score"], reverse=True)
    primary = results[0]

    return {
        "patient": {"name": req.patient_name, "id": req.patient_id, "age": req.age, "sex": req.sex},
        "primary_diagnosis": primary,
        "all_results": results,
        "database_version": "ICD-10-2024",
        "total_diseases_screened": len(DISEASE_DB),
    }
