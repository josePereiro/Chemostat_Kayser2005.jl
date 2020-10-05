import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import Chemostat_Kayser2005: KayserData, Chemostat
const Kd = KayserData
const ChU = Chemostat.Utils

## ------------------------------------------------------------------
# INFO
# Data from Kayser, Anke, Jan Weber, Volker Hecht, and Ursula Rinas. “Metabolic Flux Analysis of 
# Escherichia Coli in Glucose-Limited Continuous Culture. I. Growth-Rate-Dependent Metabolic 
# Efficiency at Steady State.” Microbiology, 151, no. 3 (2005): 693–706. https://doi.org/10.1099/mic.0.27481-0.
#

## ------------------------------------------------------------------
# TABLE1
"""
    TABLE1","Concentrations of biomass, glucose, acetate 
    and ammonium and volumetric CTR and OTR during steady-state 
    growth of E. coli TG1 at various dilution rates in 
    glucose-limited continuous cultures. 
    Data consistency was analysed by determination of the carbon 
    and nitrogen recoveries by applying reactor mass balances."
    Units: D => (1/ hr), Xv => (g/ L), Conc => (g/ L), 
    CTR and OTR => (g/ L hr), Recoveries => (%)
"""
table1 = Dict();
table1["D"] = ["Dilution", "1/ hr", 
    0.044, 0.066, 0.134, 0.150, 0.170, 0.203, 0.265,
    0.280, 0.300, 0.347, 0.375, 0.388, 0.397, 0.410, 0.415];
table1["Xv"] = ["Cell density", "gDW/ L",
    5.07, 5.05, 5.29, 5.24, 5.23, 5.41, 5.28, 
    5.53, 5.53, 5.61, 5.69, 5.88, 5.27, 3.82, 1.05];
table1["sGLC"] = ["Glucose effluent conc", "g/ L",  
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
    0.000, 0.000, 0.229, 0.295, 0.259, 0.398, 1.822, 6.048];
table1["sAC"] = ["Acetic acid effluent conc", "g/ L", 
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.048, 0.051, 0.062, 1.140, 1.500];
table1["sNH4"] = ["Ammoniate effluent conc", "g/ L", 
    1.709, 1.688, 1.652, 1.650, 1.656, 1.629, 1.632,
    1.620, 1.612, 1.617, 1.590, 1.556, 1.635, 1.835, 2.245];
# CTR: Carbon dioxide transfer rate
table1["CTR"] = ["Carbon dioxide transfer rate", "g/ L hr", 
    0.286, 0.440, 0.924, 0.915, 1.113, 1.276, 1.536,
    1.676, 1.795, 1.905, 1.888, 1.971, 2.015, 0.462, 0.207];
# Oxygen transfer rate
table1["OTR"] = ["Oxygen transfer rate", "g/ L hr", 
    0.222, 0.288, 0.576, 0.615, 0.896, 0.950, 1.044,
    1.226, 1.312, 1.357, 1.373, 1.434, 1.387, 0.330, 0.147];
# Recovery % carbon
table1["RC"] = ["Carbon recovery", "%", 
    96, 97, 101, 95, 98, 98, 93, 
    97, 97, 96, 95, 97, 93, 76, 89];
# Recovery % nitrogen
table1["RN"] = ["Nitrogen recovery", "%", 
    99, 98, 98, 98, 98, 98, 97,
    98, 98, 99, 98, 98, 97, 98, 99];


## ------------------------------------------------------------------
# TABLE1 CONVERTED
"""
    table1_converted, Converted concentrations of biomass, glucose, acetate 
    and ammonium and volumetric CTR and OTR during steady-state 
    growth of E. coli TG1 at various dilution rates in 
    glucose-limited continuous cultures. 
    Data consistency was analysed by determination of the carbon 
    and nitrogen recoveries by applying reactor mass balances."
    Met ids was set to be compatible with EColi BiGG model."
    Units: D => (1/ hr), Xv => (g/ L), Conc => (mM), 
    CTR and OTR => (g/ L hr), Recoveries => (%)
"""

table1_converted = deepcopy(table1);
# C(g/L) * 1000/ MM(g/mol) = mM
table1_converted["sGLC"][2:end] = ["mM"; table1["sGLC"][3:end] * 1e3 / 180.156];
table1_converted["sAC"][2:end] = ["mM"; table1["sAC"][3:end] * 1e3 / 60.02];
table1_converted["sNH4"][2:end] = ["mM"; table1["sNH4"][3:end] * 1e3 / 18];

## ------------------------------------------------------------------
# MEDIUM
"""
    MEDIUM, The E. coli K-12 strain TG1 
    (=DSM 6056) was grown on a defined medium with (per litre) 
    11 g glucose monohydrate (corresponding to 10 g glucose), 8 g 
    (NH4)2SO4, 0.8 g (NH4)2HPO4, 2.7 g KH2PO4, 0.35 g citric acid 
    monohydrate, 1 g MgSO4.7H2O, 12 mg ferric citrate, 0.5 mg 
    CoCl2.6H2O, 3 mg MnCl2.4H2O, 0.3 mg CuCl2.2H2O, 0.6 mg 
    H3BO3, 0?5 mg Na2MoO4.2H2O, 1?6 mg Zn(CH3COO)2.2H2O, 
    1.7 mg EDTA and '4 mg thiamine'. Glucose and magnesium sulfate 
    were autoclaved separately. Thiamine was added by sterile filtration
"""
# I take into consideration only the typical limiting metabolites
medium_converted = Dict(
    # 2 * C(8 g/L) * 1000/ MM(132.14 g/mol) + 2 * C(0.8 g/L) * 1000/ MM(132.07 g/mol) = mM
    "cNH4"  => Dict("name" => "Ammonium medium conc", "unit" => "mM", "c" => 2 * (8 * 1e3 / 132.14) + 2 * (0.8 * 1e3 / 132.07)),
    # C(10 g/L) * 1000/ MM(180.156 g/mol) = mM
    "cGLC" => Dict("name" => "glucose", "unit" => "mM", "c" => 10 * 1e3 / 180.156),
    # C(4 mg/L)/ MM(265 g/mol) = mM
    "cTHM" => Dict("name" => "Thiamine medium conc", "unit" => "mM", "c" => 4 / 265),
    # C(12 mg/L)/ MM(245 g/mol) = mM
    # "fe3dcit_e" => 12/ 245,
    # C(0.35 g/L) * 1000/ MM(192 + 18 g/mol) = mM
    "cCIT" => Dict("name" => "Citric acid medium conc", "unit" => "mM", "c" => 0.35 / (192 + 18))
);

## ------------------------------------------------------------------
# ATP DEMAND

"""
    ATPM, the energy requirement for maintenance in the 
    absence of growth and determined as mATP=2.81 mmol/ g hr
"""
atpm = 2.81;

## ------------------------------------------------------------------
# SAVING
# Table1
ChU.save_data(Kd.KAYSER_CONV_TABLE1_FILE, table1_converted)

# Medium
ChU.save_data(Kd.KAYSER_CONV_MEDIUM_FILE, medium_converted)
