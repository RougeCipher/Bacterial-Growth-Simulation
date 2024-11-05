import numpy as np 
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog

# Define base growth rates for bacteria (optimal conditions)
base_growth_rates = {
    "Escherichia coli": 0.25,              # 25% increase per hour
    "Staphylococcus aureus": 0.20,         # 20% increase per hour
    "Streptococcus pneumoniae": 0.15,      # 15% increase per hour
    "Mycobacterium tuberculosis": 0.03,    # 3% increase per hour
    "Pseudomonas aeruginosa": 0.20,        # 20% increase per hour
    "Neisseria meningitidis": 0.25,        # 25% increase per hour
    "Clostridium difficile": 0.13,         # 13% increase per hour
    "Helicobacter pylori": 0.08,           # 8% increase per hour
    "Salmonella enterica": 0.20,           # 20% increase per hour
    "Vibrio cholerae": 0.25,               # 25% increase per hour
    "Klebsiella pneumoniae": 0.20,         # 20% increase per hour
    "Acinetobacter baumannii": 0.15,       # 15% increase per hour
    "Bacillus anthracis": 0.13,            # 13% increase per hour
    "Brucella abortus": 0.03,              # 3% increase per hour
    "Listeria monocytogenes": 0.18,        # 18% increase per hour
    "Shigella dysenteriae": 0.20,          # 20% increase per hour
    "Campylobacter jejuni": 0.08,          # 8% increase per hour
    "Haemophilus influenzae": 0.18,        # 18% increase per hour
    "Bordetella pertussis": 0.08,          # 8% increase per hour
    "Legionella pneumophila": 0.08,        # 8% increase per hour
    "Rickettsia rickettsii": 0.03,         # 3% increase per hour
    "Borrelia burgdorferi": 0.02,          # 2% increase per hour
    "Enterococcus faecalis": 0.13,         # 13% increase per hour
    "Streptococcus pyogenes": 0.20,        # 20% increase per hour
    "Corynebacterium diphtheriae": 0.13,   # 13% increase per hour
    "Neisseria gonorrhoeae": 0.18,         # 18% increase per hour
    "Clostridium botulinum": 0.08,         # 8% increase per hour
    "Proteus mirabilis": 0.20,             # 20% increase per hour
    "Clostridium tetani": 0.08,            # 8% increase per hour
    "Bacteroides fragilis": 0.13,          # 13% increase per hour
    "Francisella tularensis": 0.03,        # 3% increase per hour
    "Serratia marcescens": 0.20,           # 20% increase per hour
    "Enterobacter cloacae": 0.20,          # 20% increase per hour
    "Chlamydophila pneumoniae": 0.02,      # 2% increase per hour
    "Bartonella henselae": 0.02,           # 2% increase per hour
    "Streptococcus agalactiae": 0.15,      # 15% increase per hour
    "Actinomyces israelii": 0.02,          # 2% increase per hour
    "Vibrio parahaemolyticus": 0.25,       # 25% increase per hour
    "Nocardia asteroides": 0.02,           # 2% increase per hour
    "Moraxella catarrhalis": 0.18,         # 18% increase per hour
    "Tropheryma whipplei": 0.005,          # 0.5% increase per hour
    "Coxiella burnetii": 0.02,             # 2% increase per hour
    "Leptospira interrogans": 0.02,        # 2% increase per hour
    "Staphylococcus epidermidis": 0.18,    # 18% increase per hour
    "Prevotella melaninogenica": 0.08,     # 8% increase per hour
    "Clostridium perfringens": 0.28,       # 28% increase per hour
    "Treponema pallidum": 0.01,            # 1% increase per hour
    "Chlamydia trachomatis": 0.02,         # 2% increase per hour
    "Yersinia pestis": 0.08                # 8% increase per hour
}


# Define default optimal temperature and pH for bacteria
optimal_conditions = {
    "Escherichia coli": {"temperature": 37, "pH": 7.0},
    "Staphylococcus aureus": {"temperature": 37, "pH": 7.0},
    "Streptococcus pneumoniae": {"temperature": 37, "pH": 7.4},
    "Mycobacterium tuberculosis": {"temperature": 37, "pH": 6.5},
    "Pseudomonas aeruginosa": {"temperature": 37, "pH": 7.0},
    "Neisseria meningitidis": {"temperature": 37, "pH": 7.4},
    "Clostridium difficile": {"temperature": 37, "pH": 6.5},
    "Helicobacter pylori": {"temperature": 37, "pH": 6.0},
    "Salmonella enterica": {"temperature": 37, "pH": 7.0},
    "Vibrio cholerae": {"temperature": 37, "pH": 8.5},
    "Klebsiella pneumoniae": {"temperature": 37, "pH": 7.0},
    "Acinetobacter baumannii": {"temperature": 37, "pH": 7.0},
    "Bacillus anthracis": {"temperature": 37, "pH": 7.0},
    "Brucella abortus": {"temperature": 37, "pH": 7.0},
    "Listeria monocytogenes": {"temperature": 30, "pH": 7.0},
    "Shigella dysenteriae": {"temperature": 37, "pH": 7.0},
    "Campylobacter jejuni": {"temperature": 42, "pH": 7.0},
    "Haemophilus influenzae": {"temperature": 37, "pH": 7.6},
    "Bordetella pertussis": {"temperature": 37, "pH": 7.6},
    "Legionella pneumophila": {"temperature": 35, "pH": 6.9},
    "Rickettsia rickettsii": {"temperature": 37, "pH": 7.0},
    "Borrelia burgdorferi": {"temperature": 33, "pH": 7.6},
    "Enterococcus faecalis": {"temperature": 37, "pH": 7.5},
    "Streptococcus pyogenes": {"temperature": 37, "pH": 7.4},
    "Corynebacterium diphtheriae": {"temperature": 37, "pH": 7.2},
    "Neisseria gonorrhoeae": {"temperature": 37, "pH": 7.4},
    "Clostridium botulinum": {"temperature": 35, "pH": 7.0},
    "Proteus mirabilis": {"temperature": 37, "pH": 7.0},
    "Clostridium tetani": {"temperature": 37, "pH": 7.0},
    "Bacteroides fragilis": {"temperature": 37, "pH": 7.0},
    "Francisella tularensis": {"temperature": 37, "pH": 6.9},
    "Serratia marcescens": {"temperature": 30, "pH": 7.0},
    "Enterobacter cloacae": {"temperature": 37, "pH": 7.0},
    "Chlamydophila pneumoniae": {"temperature": 37, "pH": 7.0},
    "Bartonella henselae": {"temperature": 37, "pH": 7.0},
    "Streptococcus agalactiae": {"temperature": 37, "pH": 7.4},
    "Actinomyces israelii": {"temperature": 37, "pH": 7.3},
    "Vibrio parahaemolyticus": {"temperature": 37, "pH": 8.0},
    "Nocardia asteroides": {"temperature": 37, "pH": 7.0},
    "Moraxella catarrhalis": {"temperature": 37, "pH": 7.0},
    "Tropheryma whipplei": {"temperature": 37, "pH": 7.0},
    "Coxiella burnetii": {"temperature": 37, "pH": 4.5},
    "Leptospira interrogans": {"temperature": 30, "pH": 7.4},
    "Staphylococcus epidermidis": {"temperature": 37, "pH": 7.0},
    "Prevotella melaninogenica": {"temperature": 37, "pH": 7.0},
    "Clostridium perfringens": {"temperature": 37, "pH": 7.0},
    "Treponema pallidum": {"temperature": 33, "pH": 7.2},
    "Chlamydia trachomatis": {"temperature": 37, "pH": 7.0},
    "Yersinia pestis": {"temperature": 28, "pH": 7.2}
}

# Define default parameter values for antibiotics
default_parameters = {
    "Escherichia coli": {
        "Ciprofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.2, "antibiotic_efficiency": 0.95},
        "Meropenem": {"mutation_rate": 0.001, "growth_rate": 1.2, "antibiotic_efficiency": 0.98},
        "Amoxicillin": {"mutation_rate": 0.005, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
    },
    "Staphylococcus aureus": {
        "Vancomycin": {"mutation_rate": 0.002, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Linezolid": {"mutation_rate": 0.0015, "growth_rate": 1.15, "antibiotic_efficiency": 0.92},
        "Clindamycin": {"mutation_rate": 0.002, "growth_rate": 1.1, "antibiotic_efficiency": 0.88},
    },
    "Streptococcus pneumoniae": {
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.96},
        "Levofloxacin": {"mutation_rate": 0.0015, "growth_rate": 1.05, "antibiotic_efficiency": 0.93},
        "Ceftriaxone": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.95},
    },
    "Mycobacterium tuberculosis": {
        "Rifampicin": {"mutation_rate": 0.0015, "growth_rate": 0.8, "antibiotic_efficiency": 0.93},
        "Isoniazid": {"mutation_rate": 0.001, "growth_rate": 0.9, "antibiotic_efficiency": 0.92},
        "Ethambutol": {"mutation_rate": 0.002, "growth_rate": 0.8, "antibiotic_efficiency": 0.85},
    },
    "Pseudomonas aeruginosa": {
        "Piperacillin": {"mutation_rate": 0.001, "growth_rate": 1.2, "antibiotic_efficiency": 0.86},
        "Tobramycin": {"mutation_rate": 0.0015, "growth_rate": 1.1, "antibiotic_efficiency": 0.88},
        "Ceftazidime": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
    },
    "Neisseria meningitidis": {
        "Ceftriaxone": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.96},
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.05, "antibiotic_efficiency": 0.90},
        "Moxifloxacin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.91},
    },
    "Clostridium difficile": {
        "Metronidazole": {"mutation_rate": 0.002, "growth_rate": 1.1, "antibiotic_efficiency": 0.85},
        "Fidaxomicin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.94},
        "Vancomycin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
    },
    "Helicobacter pylori": {
        "Clarithromycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
        "Amoxicillin": {"mutation_rate": 0.003, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Tetracycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.82},
    },
    "Salmonella enterica": {
        "Ciprofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.92},
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.05, "antibiotic_efficiency": 0.86},
        "Cefixime": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
    },
    "Vibrio cholerae": {
        "Doxycycline": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Furazolidone": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
    },
    "Klebsiella pneumoniae": {
        "Carbapenems": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Ceftriaxone": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Polymyxin B": {"mutation_rate": 0.003, "growth_rate": 1.0, "antibiotic_efficiency": 0.80},
    },
    "Acinetobacter baumannii": {
        "Colistin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Meropenem": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Tigecycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
    },
    "Bacillus anthracis": {
        "Ciprofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.2, "antibiotic_efficiency": 0.95},
        "Doxycycline": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
    },
    "Brucella abortus": {
        "Doxycycline": {"mutation_rate": 0.01, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Rifampicin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.93},
        "Streptomycin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
    },
    "Listeria monocytogenes": {
        "Ampicillin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Gentamicin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Trimethoprim": {"mutation_rate": 0.003, "growth_rate": 1.0, "antibiotic_efficiency": 0.80},
    },
    "Shigella dysenteriae": {
        "Ciprofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.93},
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Cefotaxime": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
    },
    "Campylobacter jejuni": {
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.05, "antibiotic_efficiency": 0.85},
        "Ciprofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Erythromycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.82},
    },
    "Haemophilus influenzae": {
        "Amoxicillin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Ceftriaxone": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Levofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
    },
    "Bordetella pertussis": {
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Clarithromycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
        "Erythromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
    },
    "Legionella pneumophila": {
        "Levofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.91},
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.05, "antibiotic_efficiency": 0.88},
        "Moxifloxacin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Rickettsia rickettsii": {
        "Doxycycline": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Chloramphenicol": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Tetracycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
    },
    "Borrelia burgdorferi": {
        "Doxycycline": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Ceftriaxone": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Amoxicillin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Enterococcus faecalis": {
        "Vancomycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Ampicillin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Daptomycin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
    },
    "Streptococcus pyogenes": {
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.98},
        "Clindamycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Cephalexin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
    },
    "Corynebacterium diphtheriae": {
        "Penicillin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Azithromycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Erythromycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Neisseria gonorrhoeae": {
        "Ceftriaxone": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.05, "antibiotic_efficiency": 0.90},
        "Spectinomycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
    },
    "Clostridium botulinum": {
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.93},
        "Metronidazole": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Antitoxin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.98},
    },
    "Proteus mirabilis": {
        "Ampicillin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Ceftriaxone": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Sulfamethoxazole": {"mutation_rate": 0.003, "growth_rate": 1.0, "antibiotic_efficiency": 0.80},
    },
    "Clostridium tetani": {
        "Metronidazole": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.93},
        "Erythromycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Bacteroides fragilis": {
        "Metronidazole": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Imipenem": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Clindamycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Francisella tularensis": {
        "Streptomycin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.95},
        "Gentamicin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Ciprofloxacin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Serratia marcescens": {
        "Cefepime": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Piperacillin-tazobactam": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Aztreonam": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
    },
    "Enterobacter cloacae": {
        "Carbapenems": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Cefepime": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
        "Amikacin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Chlamydophila pneumoniae": {
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Doxycycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Rifampin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Bartonella henselae": {
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Doxycycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Rifampin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Streptococcus agalactiae": {
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Ampicillin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Vancomycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
    },
    "Actinomyces israelii": {
        "Penicillin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Doxycycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Clindamycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Vibrio parahaemolyticus": {
        "Doxycycline": {"mutation_rate": 0.0015, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Ciprofloxacin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.92},
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
    },
    "Nocardia asteroides": {
        "Trimethoprim-sulfamethoxazole": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Amikacin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Imipenem": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Moraxella catarrhalis": {
        "Amoxicillin-clavulanate": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Cefuroxime": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Cefixime": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
    },
    "Tropheryma whipplei": {
        "Ceftriaxone": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.95},
        "Trimethoprim-sulfamethoxazole": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Penicillin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Coxiella burnetii": {
        "Doxycycline": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Hydroxychloroquine": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Chloramphenicol": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Leptospira interrogans": {
        "Doxycycline": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.93},
        "Ceftriaxone": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Staphylococcus epidermidis": {
        "Vancomycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
        "Linezolid": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Daptomycin": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Prevotella melaninogenica": {
        "Metronidazole": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Clindamycin": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
        "Amoxicillin-clavulanate": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
    "Clostridium perfringens": {
        "Penicillin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.93},
        "Clindamycin": {"mutation_rate": 0.002, "growth_rate": 1.1, "antibiotic_efficiency": 0.90},
        "Metronidazole": {"mutation_rate": 0.002, "growth_rate": 1.1, "antibiotic_efficiency": 0.88},
    },
    "Treponema pallidum": {
        "Penicillin G": {"mutation_rate": 0.001, "growth_rate": 1.0, "antibiotic_efficiency": 0.98},
        "Doxycycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Tetracycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.88},
    },
    "Chlamydia trachomatis": {
        "Azithromycin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.85},
        "Doxycycline": {"mutation_rate": 0.002, "growth_rate": 1.0, "antibiotic_efficiency": 0.92},
        "Ofloxacin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.87},
    },
    "Yersinia pestis": {
        "Streptomycin": {"mutation_rate": 0.001, "growth_rate": 1.1, "antibiotic_efficiency": 0.95},
        "Doxycycline": {"mutation_rate": 0.002, "growth_rate": 1.05, "antibiotic_efficiency": 0.92},
        "Ciprofloxacin": {"mutation_rate": 0.0015, "growth_rate": 1.0, "antibiotic_efficiency": 0.90},
    },
}


# Functions for antibiotic resistance mechanisms and environmental adjustments
def apply_efflux_pump_effect(antibiotic_efficiency, time_step, efflux_rate=0.01):
    return max(0, antibiotic_efficiency - (efflux_rate * time_step))

def apply_target_modification_effect(antibiotic_efficiency, mutation_rate, target_modification_factor=0.1):
    return max(0, antibiotic_efficiency * (1 - (mutation_rate * target_modification_factor)))

def adjust_antibiotic_efficiency(antibiotic_efficiency, time_step, mutation_rate, temperature, pH):
    antibiotic_efficiency = apply_efflux_pump_effect(antibiotic_efficiency, time_step)
    antibiotic_efficiency = apply_target_modification_effect(antibiotic_efficiency, mutation_rate)
    antibiotic_efficiency = adjust_antibiotic_efficiency_for_environment(antibiotic_efficiency, temperature, pH)
    antibiotic_efficiency = np.clip(antibiotic_efficiency, 0, 1)  # Ensure efficiency is within [0,1]
    return antibiotic_efficiency

def adjust_antibiotic_efficiency_for_environment(antibiotic_efficiency, temperature, pH):
    if temperature < 15 or temperature > 45:  # Arbitrary stress limits
        antibiotic_efficiency *= 0.7
    if pH < 5 or pH > 9:  # Arbitrary pH stress limits
        antibiotic_efficiency *= 0.8
    return antibiotic_efficiency

def adjust_growth_rate_for_environment(growth_rate, temperature, optimal_temp, pH, optimal_pH):
    # Adjust for temperature
    if temperature < optimal_temp - 5 or temperature > optimal_temp + 5:
        growth_rate *= 0.5  # Halve growth rate if temperature is too high or too low
    elif temperature < optimal_temp - 2 or temperature > optimal_temp + 2:
        growth_rate *= 0.8  # Slight reduction for suboptimal temperature
    # Adjust for pH
    if pH < optimal_pH - 1 or pH > optimal_pH + 1:
        growth_rate *= 0.7  # Reduce growth rate by 30% if pH is not optimal
    return growth_rate

# Functions for dynamic environmental conditions
def dynamic_temperature(t, base_temp, amplitude=5, period=24):
    """Sinusoidal temperature variation."""
    temp = base_temp + amplitude * np.sin(2 * np.pi * t / period)
    temp = np.clip(temp, 0, 50)  # Ensure temperature stays within realistic bounds
    return temp

def dynamic_pH(t, base_pH, amplitude=0.3, period=24):
    """Sinusoidal pH variation."""
    pH = base_pH + amplitude * np.sin(2 * np.pi * t / period)
    pH = np.clip(pH, 1, 14)  # Ensure pH stays within realistic bounds
    return pH

# Functions for growth rate adjustments
def adjust_growth_for_oxygen(growth_rate, oxygen_level, is_aerobic=True):
    oxygen_level = np.clip(oxygen_level, 0, 1)  # Ensure oxygen level is within [0,1]
    if is_aerobic:
        return growth_rate * oxygen_level  # Aerobic bacteria need oxygen
    else:
        return growth_rate * (1 - oxygen_level)  # Anaerobic bacteria grow better with less oxygen

def adjust_growth_for_moisture(growth_rate, moisture_level):
    moisture_level = np.clip(moisture_level, 0, 1)  # Ensure moisture level is within [0,1]
    return growth_rate * moisture_level

def adjust_growth_for_specific_nutrients(growth_rate, nutrients):
    nutrient_factors = {
        "glucose": 1.2,
        "amino_acids": 1.1,
        "vitamins": 1.05
    }
    for nutrient, level in nutrients.items():
        level = np.clip(level, 0, 1)  # Ensure nutrient level is within [0,1]
        growth_rate *= nutrient_factors.get(nutrient, 1) * level
    return growth_rate

def deplete_nutrients(nutrients, population_size, depletion_rate=0.00001):
    for nutrient in nutrients:
        nutrients[nutrient] -= depletion_rate * population_size
        nutrients[nutrient] = max(0, nutrients[nutrient])  # Prevent negative levels
    return nutrients

def adjust_growth_for_waste(growth_rate, waste_level, toxicity_threshold=0.5):
    if waste_level > toxicity_threshold:
        growth_rate *= max(0, (1 - (waste_level - toxicity_threshold)))
    return growth_rate

def accumulate_waste(waste_level, population_size, production_rate=0.00001):
    return waste_level + production_rate * population_size

def dynamic_mutation_rate(base_mutation_rate, antibiotic_efficiency, stress_factor=0.1):
    return base_mutation_rate * (1 + antibiotic_efficiency * stress_factor)

def combined_stress_effect(growth_rate, temperature, pH, nutrients, oxygen_level, moisture_level, waste_level, is_aerobic, optimal_temp, optimal_pH):
    growth_rate = adjust_growth_rate_for_environment(growth_rate, temperature, optimal_temp, pH, optimal_pH)
    growth_rate = adjust_growth_for_specific_nutrients(growth_rate, nutrients)
    growth_rate = adjust_growth_for_oxygen(growth_rate, oxygen_level, is_aerobic)
    growth_rate = adjust_growth_for_moisture(growth_rate, moisture_level)
    growth_rate = adjust_growth_for_waste(growth_rate, waste_level)
    growth_rate = max(growth_rate, 0)  # Ensure growth rate is non-negative
    return growth_rate

def decay_antibiotic_concentration(antibiotic_efficiency, time_step, decay_rate=0.02):
    return antibiotic_efficiency * np.exp(-decay_rate * time_step)

def get_growth_phase(t, lag_time=2, exponential_time=12, stationary_time=18):
    if t < lag_time:
        return "lag"
    elif t < exponential_time:
        return "exponential"
    elif t < stationary_time:
        return "stationary"
    else:
        return "death"

def adjust_growth_rate_for_phase(growth_rate, phase):
    if phase == "lag":
        return 0.1 * growth_rate  # Slow growth in lag phase
    elif phase == "exponential":
        return growth_rate  # Full growth rate in exponential phase
    elif phase == "stationary":
        return 0.05 * growth_rate  # Reduced growth in stationary phase
    elif phase == "death":
        return -0.1 * abs(growth_rate)  # Negative growth in death phase

def plot_antibiotic_decay(time, antibiotic_efficiencies):
    plt.figure(figsize=(10, 5))
    plt.plot(time, antibiotic_efficiencies, label="Antibiotic Efficiency")
    plt.xlabel("Time (hours)")
    plt.ylabel("Antibiotic Efficiency")
    plt.title("Antibiotic Decay Over Time")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_population_diversity(time, non_resistant_population, resistant_population):
    plt.figure(figsize=(10, 5))
    plt.plot(time, non_resistant_population, label="Non-Resistant Population")
    plt.plot(time, resistant_population, label="Resistant Population", linestyle='--')
    plt.xlabel("Time (hours)")
    plt.ylabel("Population Size")
    plt.title("Population Diversity Over Time")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_environmental_conditions(time, temperature_series, pH_series):
    fig, ax1 = plt.subplots(figsize=(10, 5))

    ax1.set_xlabel("Time (hours)")
    ax1.set_ylabel("Temperature (°C)", color='red')
    ax1.plot(time, temperature_series, label="Temperature (°C)", color='red')
    ax1.tick_params(axis='y', labelcolor='red')

    ax2 = ax1.twinx()
    ax2.set_ylabel("pH Level", color='blue')
    ax2.plot(time, pH_series, label="pH Level", color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')

    plt.title("Environmental Conditions Over Time")
    fig.tight_layout()
    plt.grid(True)
    plt.show()

def plot_waste_accumulation(time, waste_levels):
    plt.figure(figsize=(10, 5))
    plt.plot(time, waste_levels, label="Waste Level")
    plt.xlabel("Time (hours)")
    plt.ylabel("Waste Level")
    plt.title("Waste Accumulation Over Time")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_nutrient_depletion(time, nutrients_over_time):
    plt.figure(figsize=(10, 5))
    for nutrient, levels in nutrients_over_time.items():
        plt.plot(time, levels, label=f"{nutrient.capitalize()} Level")
    plt.xlabel("Time (hours)")
    plt.ylabel("Nutrient Level")
    plt.title("Nutrient Depletion Over Time")
    plt.legend()
    plt.grid(True)
    plt.show()

def calculate_phase_metrics(population_data, time):
    lag_phase_duration = time[np.argmax(population_data > population_data[0]*1.01)]  # Time until 1% growth
    peak_population = max(population_data)
    return lag_phase_duration, peak_population

def get_user_inputs():
    root = tk.Tk()
    root.withdraw()
    try:
        species = simpledialog.askstring("Input", "Select bacterial species (Escherichia coli, Staphylococcus aureus, Yersinia pestis): ")
        while species not in base_growth_rates:
            species = simpledialog.askstring("Error", "Invalid species. Please select from (Escherichia coli, Staphylococcus aureus, Yersinia pestis): ")

        antibiotic_options = list(default_parameters.get(species, {}).keys())
        antibiotic = simpledialog.askstring("Input", f"Select antibiotic (if any) from {antibiotic_options} (or leave blank): ")
        if antibiotic not in antibiotic_options:
            antibiotic = ""

        nutrient_environment = simpledialog.askstring("Input", "Nutrient Environment (Rich/Deficient/Neutral): ")
        selected_nutrients = simpledialog.askstring("Input", "Enter nutrients present (comma-separated, e.g., glucose, amino_acids, vitamins): ").split(',')
        nutrient_levels = {}
        for nutrient in selected_nutrients:
            nutrient = nutrient.strip()
            level = float(simpledialog.askstring("Input", f"Enter level for {nutrient} (0 to 1): "))
            while not (0 <= level <= 1):
                level = float(simpledialog.askstring("Error", f"Invalid level. Enter level for {nutrient} (0 to 1): "))
            nutrient_levels[nutrient] = level

        initial_temperature = float(simpledialog.askstring("Input", "Enter initial temperature in °C (0 to 50): "))
        while not (0 <= initial_temperature <= 50):
            initial_temperature = float(simpledialog.askstring("Error", "Invalid temperature. Enter initial temperature in °C (0 to 50): "))

        initial_pH = float(simpledialog.askstring("Input", "Enter initial pH level (1 to 14): "))
        while not (1 <= initial_pH <= 14):
            initial_pH = float(simpledialog.askstring("Error", "Invalid pH. Enter initial pH level (1 to 14): "))

        oxygen_level = float(simpledialog.askstring("Input", "Enter oxygen level (0 to 1): "))
        while not (0 <= oxygen_level <= 1):
            oxygen_level = float(simpledialog.askstring("Error", "Invalid level. Enter oxygen level (0 to 1): "))

        moisture_level = float(simpledialog.askstring("Input", "Enter moisture level (0 to 1): "))
        while not (0 <= moisture_level <= 1):
            moisture_level = float(simpledialog.askstring("Error", "Invalid level. Enter moisture level (0 to 1): "))

        is_aerobic_input = simpledialog.askstring("Input", "Is the bacterium aerobic? (Yes/No): ").lower()
        while is_aerobic_input not in ['yes', 'no']:
            is_aerobic_input = simpledialog.askstring("Error", "Invalid input. Is the bacterium aerobic? (Yes/No): ").lower()
        is_aerobic = is_aerobic_input == 'yes'

        initial_population = simpledialog.askstring("Input", "Enter initial population size [default 1000]: ")
        initial_population = int(initial_population) if initial_population else 1000

        resistant_population = simpledialog.askstring("Input", "Enter initial resistant population size [default 0]: ")
        resistant_population = int(resistant_population) if resistant_population else 0

        time_hours = simpledialog.askstring("Input", "Enter the timeframe for simulation in hours [default 24]: ")
        time_hours = int(time_hours) if time_hours else 24

    except Exception as e:
        print("Error in input:", e)
        return None

    return {
        "species": species.strip(),
        "antibiotic": antibiotic.strip(),
        "nutrient_environment": nutrient_environment.strip(),
        "nutrients": nutrient_levels,
        "initial_temperature": initial_temperature,
        "initial_pH": initial_pH,
        "oxygen_level": oxygen_level,
        "moisture_level": moisture_level,
        "is_aerobic": is_aerobic,
        "initial_population": initial_population,
        "resistant_population": resistant_population,
        "time_hours": time_hours
    }

def simulate_bacterial_growth(parameters):
    species = parameters["species"]
    time_hours = parameters["time_hours"]
    time = np.linspace(0, time_hours, num=time_hours + 1)
    non_resistant_population = np.zeros(len(time))
    resistant_population = np.zeros(len(time))
    non_resistant_population[0] = parameters["initial_population"] - parameters["resistant_population"]
    resistant_population[0] = parameters["resistant_population"]
    antibiotic_efficiencies = np.zeros(len(time))

    # Initialize lists to store environmental data
    temperature_series = []
    pH_series = []
    waste_levels = []
    # Initialize dictionary to store nutrient levels over time
    nutrients_over_time = {nutrient: [] for nutrient in parameters["nutrients"].keys()}

    # Get optimal temperature and pH for the species
    optimal_temp = optimal_conditions.get(species, {}).get("temperature", 37)
    optimal_pH = optimal_conditions.get(species, {}).get("pH", 7.0)

    # Base growth rate
    base_growth_rate = base_growth_rates.get(species, 0.1)

    # Waste level
    waste_level = 0

    # Nutrients
    nutrients = parameters["nutrients"].copy()  # Copy to avoid modifying original

    for t_idx, t in enumerate(time):
        # Dynamic environmental conditions
        temperature = dynamic_temperature(t, parameters["initial_temperature"])
        pH = dynamic_pH(t, parameters["initial_pH"])

        # Record environmental conditions
        temperature_series.append(temperature)
        pH_series.append(pH)
        waste_levels.append(waste_level)
        for nutrient in nutrients:
            nutrients_over_time[nutrient].append(nutrients[nutrient])

        # Adjust growth rate based on combined stressors
        current_growth_rate = combined_stress_effect(
            base_growth_rate,
            temperature,
            pH,
            nutrients,
            parameters["oxygen_level"],
            parameters["moisture_level"],
            waste_level,
            parameters["is_aerobic"],
            optimal_temp,
            optimal_pH
        )

        # Adjust growth rate based on growth phase
        phase = get_growth_phase(t)
        current_growth_rate = adjust_growth_rate_for_phase(current_growth_rate, phase)

        # Antibiotic parameters
        if parameters["antibiotic"]:
            antibiotic_params = default_parameters.get(species, {}).get(parameters["antibiotic"], {})
            base_mutation_rate = antibiotic_params.get("mutation_rate", 0.001)
            initial_antibiotic_efficiency = antibiotic_params.get("antibiotic_efficiency", 0.9)
        else:
            base_mutation_rate = 0
            initial_antibiotic_efficiency = 0

        # Decay antibiotic concentration over time
        if parameters["antibiotic"]:
            antibiotic_efficiency = decay_antibiotic_concentration(initial_antibiotic_efficiency, t)
            # Adjust antibiotic efficiency due to resistance mechanisms and environmental factors
            antibiotic_efficiency = adjust_antibiotic_efficiency(antibiotic_efficiency, t_idx, base_mutation_rate, temperature, pH)
            antibiotic_efficiency = np.clip(antibiotic_efficiency, 0, 1)
            antibiotic_efficiencies[t_idx] = antibiotic_efficiency
        else:
            antibiotic_efficiency = 0

        # Dynamic mutation rate
        mutation_rate = dynamic_mutation_rate(base_mutation_rate, antibiotic_efficiency)

        # Mutation leading to resistance
        if parameters["antibiotic"] and t_idx > 0:
            new_resistant_bacteria = non_resistant_population[t_idx - 1] * mutation_rate * antibiotic_efficiency
            non_resistant_population[t_idx - 1] -= new_resistant_bacteria
            resistant_population[t_idx - 1] += new_resistant_bacteria

        # Non-resistant population growth
        if t_idx > 0:
            if parameters["antibiotic"]:
                growth_factor = 1 + current_growth_rate * (1 - antibiotic_efficiency)
            else:
                growth_factor = 1 + current_growth_rate
            non_resistant_population[t_idx] = non_resistant_population[t_idx - 1] * growth_factor

        # Resistant population growth
        if t_idx > 0:
            resistant_population[t_idx] = resistant_population[t_idx - 1] * (1 + current_growth_rate)

        # Ensure population values are non-negative
        non_resistant_population[t_idx] = max(non_resistant_population[t_idx], 0)
        resistant_population[t_idx] = max(resistant_population[t_idx], 0)

        # Deplete nutrients and accumulate waste
        total_population = non_resistant_population[t_idx] + resistant_population[t_idx]
        nutrients = deplete_nutrients(nutrients, total_population)
        waste_level = accumulate_waste(waste_level, total_population)

    total_population = non_resistant_population + resistant_population

    # Calculate phase metrics
    lag_phase_duration, peak_population = calculate_phase_metrics(total_population, time)

    # Compute growth rates
    with np.errstate(divide='ignore', invalid='ignore'):
        growth_rates = np.divide(np.diff(total_population), total_population[:-1])
        growth_rates = np.nan_to_num(growth_rates)  # Replace nan with 0

    # Determine phases based on growth rates
    phases = np.array([''] * len(time))
    phases[0] = 'lag'  # Start with lag phase
    for i in range(1, len(time)):
        gr = growth_rates[i-1]
        if gr >= 0.02:
            phases[i] = 'exponential'
        elif 0 <= gr < 0.02:
            phases[i] = 'lag'
        elif -0.01 <= gr < 0:
            phases[i] = 'stationary'
        elif gr < -0.01:
            phases[i] = 'death'
        else:
            phases[i] = 'stationary'  # default

    # Calculate phase durations
    phase_durations = {'lag': 0, 'exponential': 0, 'stationary': 0, 'death': 0}
    current_phase = phases[0]
    phase_start_time = time[0]
    for i in range(1, len(time)):
        if phases[i] != current_phase:
            phase_end_time = time[i-1]
            duration = phase_end_time - phase_start_time
            phase_durations[current_phase] += duration
            current_phase = phases[i]
            phase_start_time = time[i]
    # Add the last phase duration
    phase_end_time = time[-1]
    duration = phase_end_time - phase_start_time
    phase_durations[current_phase] += duration

    # Detect waste inhibition time
    toxicity_threshold = 0.5  # from adjust_growth_for_waste function
    waste_level_array = np.array(waste_levels)
    waste_inhibition_time = None
    if np.any(waste_level_array > toxicity_threshold):
        waste_inhibition_index = np.argmax(waste_level_array > toxicity_threshold)
        waste_inhibition_time = time[waste_inhibition_index]
    else:
        waste_inhibition_time = 'Not reached'

    # Plotting
    plot_growth(time, total_population, non_resistant_population, resistant_population, species)
    if parameters["antibiotic"]:
        plot_antibiotic_decay(time, antibiotic_efficiencies)
    plot_population_diversity(time, non_resistant_population, resistant_population)
    plot_environmental_conditions(time, temperature_series, pH_series)
    plot_waste_accumulation(time, waste_levels)
    plot_nutrient_depletion(time, nutrients_over_time)

    # Print simulation summary
    print("---------------------------------------")
    print("Simulation Summary")
    print("---------------------------------------")
    print(f"Lag Phase Duration: {lag_phase_duration:.1f} hours")
    peak_population_formatted = "{:.1f} x 10^{}".format(
        peak_population / (10 ** int(np.log10(peak_population))),
        int(np.log10(peak_population))
    )
    print(f"Peak Population: {peak_population_formatted} CFU")
    print("Phase Durations:")
    for phase in ['lag', 'exponential', 'stationary', 'death']:
        duration = phase_durations.get(phase, 0)
        print(f"  - {phase.capitalize()}: {duration:.1f} hours")
    print("Environmental Summary:")
    # Temperature variation
    base_temp = parameters["initial_temperature"]
    temp_amplitude = 5  # from dynamic_temperature function
    print(f"  - Temperature variation: ±{temp_amplitude}°C around {base_temp}°C")
    # pH variation
    base_pH = parameters["initial_pH"]
    pH_amplitude = 0.3  # from dynamic_pH function
    print(f"  - pH variation: ±{pH_amplitude} around pH {base_pH}")
    # Nutrient Depletion Rate
    depletion_rate = 0.00001  # from deplete_nutrients function
    print(f"  - Nutrient Depletion Rate: {depletion_rate}")
    # Waste Accumulation
    if waste_inhibition_time != 'Not reached':
        print(f"  - Waste Accumulation: Growth-inhibiting after {waste_inhibition_time:.1f} hours")
    else:
        print(f"  - Waste Accumulation: Toxicity threshold not reached")

def plot_growth(time, total_population, non_resistant_population, resistant_population, species):
    plt.figure(figsize=(12, 7))
    plt.plot(time, total_population, label=f"Total {species} Population")
    plt.plot(time, non_resistant_population, label=f"Non-Resistant {species} Population", linestyle='--')
    plt.plot(time, resistant_population, label=f"Resistant {species} Population", linestyle='-.')
    plt.xlabel("Time (hours)")
    plt.ylabel("Bacterial Population")
    plt.title(f"Bacterial Growth Simulation for {species}")
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    parameters = get_user_inputs()
    if parameters:
        simulate_bacterial_growth(parameters)
    else:
        print("Simulation aborted due to input errors.")

if __name__ == "__main__":
    main()
