## ------------------------------------------------------------------
# maps between Kayser2005 https://doi.org/10.1099/mic.0.27481-0.
# data and the model

function load_mets_map()
    Kd_mets_map = Dict()
    Kd_mets_map["GLC"] = "glc_DASH_D_e"
    Kd_mets_map["THM"] = "thm_e"
    Kd_mets_map["NH4"] = "nh4_e"
    Kd_mets_map["CIT"] = "cit_e"
    Kd_mets_map["AC"] = "ac_e"
    return Kd_mets_map
end

load_exch_met_map() = load_data(EXCH_MET_MAP_FILE; verbose = false)

function load_iders_to_plot_map()
    iders_to_plot_map = Dict()
    mets_map = load_mets_map()
    for Kd_ider in Kd.iders_to_plot
        model_ider = 
        iders_to_plot_map[]
    end
end

## ------------------------------------------------------------------
# The intakes bounds of the network are determined by the 
# medium concentration in the Chemostat model (see Cossios paper)
# This is a base medium for modeling
function load_base_intake_info()
    base_intake_info = Dict()
    # form fba
    base_intake_info["EX_o2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_pi_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_ca2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_cl_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_cobalt2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_cu2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_fe2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_k_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_mg2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_mn2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_mobd_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_ni2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_so4_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_zn2_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_cbl1_e"] = Dict("c" => MAX_CONC, "lb" => -ABS_MAX_BOUND)

    # from kayser
    base_intake_info["EX_glc__D_e"] = Dict("c" => Kd.medium["cGLC"]["c"], "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_nh4_e"] = Dict("c" => Kd.medium["cNH4"]["c"], "lb" => -ABS_MAX_BOUND)
    base_intake_info["EX_cit_e"] = Dict("c" => Kd.medium["cCIT"]["c"], "lb" => -ABS_MAX_BOUND)
    # cTHM: thiamine missing from model. TODO: I can add an ne exchange

    return base_intake_info
end

## ------------------------------------------------------------------
function load_beg_rxns_map()
    beg_rxns_map = Dict()
    beg_rxns_map["carbonyl reductase (NADPH)"] = ["P5CR"]
    beg_rxns_map["alcohol dehydrogenase (NADP+)"] = ["ALCD19", "ALCD2x"]
    beg_rxns_map["quinate/shikimate dehydrogenase"] = ["SHK3Dr"]
    beg_rxns_map["malate dehydrogenase (decarboxylating)"] = ["MDH","MDH2", "MDH3"]
    beg_rxns_map["3alpha-hydroxysteroid dehydrogenase (B-specific)"] = [""]
    beg_rxns_map["2-hydroxy-3-oxopropionate reductase"] = [""]
    beg_rxns_map["glucose dehydrogenase (acceptor)"] = ["G6PDH2r", "UDPGD"]
    beg_rxns_map["cellobiose dehydrogenase (acceptor)"] = [""]
    beg_rxns_map["peroxidase"] = ["THIORDXi"]
    beg_rxns_map["catechol 2,3-dioxygenase"] = [""]
    beg_rxns_map["arachidonate 8-lipoxygenase"] = [""]
    beg_rxns_map["calcidiol 1-monooxygenase"] = [""]
    beg_rxns_map["nitric-oxide synthase"] = [""]
    beg_rxns_map["phenylalanine 4-monooxygenase"] = [""]
    beg_rxns_map["tryptophan 5-monooxygenase"] = ["MTRPOX"]
    beg_rxns_map["Carboxylate reductase"] = ["P5CR"]
    beg_rxns_map["arsenate reductase (donor)"] = ["ASR"]
    beg_rxns_map["biliverdin reductase"] = [""]
    beg_rxns_map["15-oxoprostaglandin 13-oxidase"] = [""]
    beg_rxns_map["coproporphyrinogen oxidase"] = ["CPPPGO","CPPPGO2","PPPGO","PPPGO3"]
    beg_rxns_map["long-chain-acyl-CoA dehydrogenase"] = ["HACD8","HADPCOADH3","OXCOAHDH"]
    beg_rxns_map["butyryl-CoA dehydrogenase"] = ["ACOAD1f"]
    beg_rxns_map["acyl-CoA dehydrogenase"] = ["ACOAD2f","ACOAD3f","ACOAD4f",
                                    "ACOAD5f","ACOAD6f","ACOAD7f","ACOAD8f",
                                    "HACD1","HACD2","HACD3","HACD4","HACD5","HACD6",
                                    "HACD7"]
    beg_rxns_map["L-amino-acid oxidase"] = ["ASPO3","ASPO4","ASPO5","ASPO6"]
    beg_rxns_map["amine oxidase (flavin-containing)"] = ["PYAM5PO"]
    beg_rxns_map["methylenetetrahydrofolate reductase [NAD(P)H]"] = ["MTHFR2"]
    beg_rxns_map["formyltetrahydrofolate dehydrogenase"] = ["MTHFD"]
    beg_rxns_map["sarcosine oxidase"] = ["SARCOX"]
    beg_rxns_map["nitrate reductase (NADH)"] = [""]
    beg_rxns_map["nitrite reductase (NO-forming)"] = [""]
    beg_rxns_map["nitrate reductase"] = ["NO3R1bpp","NO3R2bpp","NO3R2pp","NO3R2bpp"]
    beg_rxns_map["trypanothione-disulfide reductase"] = ["TDSR1","TDSR2"]
    beg_rxns_map["glutathione-disulfide reductase"] = [""]
    beg_rxns_map["thioredoxin-disulfide reductase"] = [""]
    beg_rxns_map["thiol oxidase"] = [""]
    beg_rxns_map["nitrate reductase (cytochrome)"] = [""]
    beg_rxns_map["aspartate carbamoyltransferase"] = ["ASPCT"]
    beg_rxns_map["serine O-acetyltransferase"] = ["SERAT"]
    beg_rxns_map["protein-glutamine gamma-glutamyltransferase"] = [""]
    beg_rxns_map["gamma-glutamyltransferase"] = ["CRNBTCT"]
    beg_rxns_map["citrate (Si)-synthase"] = ["CS"]
    beg_rxns_map["kaempferol 3-O-galactosyltransferase"] = ["GALT1"]
    beg_rxns_map["NAD+ ADP-ribosyltransferase"] = ["NNDMBRT"]
    beg_rxns_map["di-trans,poly-cis-decaprenylcistransferase"] = ["ACGAMT"]
    beg_rxns_map["cystathionine gamma-synthase"] = [""]
    beg_rxns_map["adenosine kinase"] = ["ADNK1"]
    beg_rxns_map["glycerate kinase"] = ["GLYCK","GLYCK2"]
    beg_rxns_map["galactokinase"] = ["GALKr"]
    beg_rxns_map["[pyruvate dehydrogenase (acetyl-transferring)] kinase"] = [""]
    beg_rxns_map["guanylate kinase"] = ["GK1"]
    beg_rxns_map["FMN adenylyltransferase"] = ["FMNAT"]
    beg_rxns_map["tRNA adenylyltransferase"] = [""]
    beg_rxns_map["aryl sulfotransferase"] = [""]
    beg_rxns_map["aminoacyl-tRNA hydrolase"] = [""]
    beg_rxns_map["carboxymethylenebutenolidase"] = [""]
    beg_rxns_map["ubiquitin thiolesterase"] = [""]
    beg_rxns_map["fructose-bisphosphatase"] = ["FBP","FBA"]
    beg_rxns_map["[phosphorylase] phosphatase"] = [""]
    beg_rxns_map["phosphoglycolate phosphatase"] = ["PGLYCP"]
    beg_rxns_map["protein-tyrosine-phosphatase"] = ["TYRPpp"]
    beg_rxns_map["inositol-polyphosphate 5-phosphatase"] = ["MI1PP"]
    beg_rxns_map["3',5'-cyclic-GMP phosphodiesterase"] = ["23PDE2pp","23PDE5pp",
                                                        "23PDE7pp","23PDE9pp",
                                                        "PDE1","PDE4"]
    beg_rxns_map["beta-glucosidase"] = ["MLTG1","MLTG2","MLTG3","MLTG4","MLTG5"]
    beg_rxns_map["beta-glucuronidase"] = [""]
    beg_rxns_map["glucosylceramidase"] = [""]
    beg_rxns_map["cyclomaltodextrinase"] = [""]
    beg_rxns_map["alpha-N-arabinofuranosidase"] = [""]
    beg_rxns_map["purine nucleosidase"] = ["AMPN","AHCYSNS","CMPN","MTAN","NMNN"]
    beg_rxns_map["rRNA N-glycosylase"] = [""]
    beg_rxns_map["NAD+ nucleosidase"] = ["NADN"]
    beg_rxns_map["Xaa-Pro aminopeptidase"] = ["AMPTASEPG","AMPTASECG"]
    beg_rxns_map["dipeptidyl-peptidase I"] = ["ALAALAD","UM4PCP"]
    beg_rxns_map["peptidyl-dipeptidase A"] = [""]
    beg_rxns_map["coagulation factor Xa"] = [""]
    beg_rxns_map["t-Plasminogen activator"] = [""]
    beg_rxns_map["cathepsin B"] = [""]
    beg_rxns_map["envelysin"] = [""]
    beg_rxns_map["amidase"] = ["AGM3PA","AGM3PApp","AGM4PA","AGM4PApp","AM3PA",
                                "AM4PA","GSPMDA","NMNDA","NNAM"]
    beg_rxns_map["formamidase"] = [""]
    beg_rxns_map["arginase"] = [""]
    beg_rxns_map["guanidinoacetase"] = ["GUAD"]
    beg_rxns_map["apyrase"] = [""]
    beg_rxns_map["phloretin hydrolase"] = [""]
    beg_rxns_map["Orotidine-5'-phosphate decarboxylase"] = ["OMPDC"]
    beg_rxns_map["4-Hydroxybenzoate decarboxylase"] = ["OPHBDC"]
    beg_rxns_map["Threonine aldolase"] = ["THRAi","THRA2i"]
    beg_rxns_map["enoyl-CoA hydratase"] = ["DHACOAH","ECOAH1","ECOAH2","ECOAH3",
                        "ECOAH4","ECOAH5","ECOAH6","ECOAH7","ECOAH8"]
    beg_rxns_map["Uroporphyrinogen-III synthase"] = ["UPP3S"]
    beg_rxns_map["dihydroxy-acid dehydratase"] = ["DHAD1","DHAD2"]
    beg_rxns_map["pectin lyase"] = [""]
    beg_rxns_map["DNA-(apurinic or apyrimidinic site) lyase"] = [""]
    beg_rxns_map["lactoylglutathione lyase"] = ["LGTHL"]
    beg_rxns_map["guanylate cyclase"] = ["GUACYC"]
    beg_rxns_map["dTDP-4-dehydrorhamnose 3,5-epimerase"] = ["TDPDRE"]
    beg_rxns_map["UDP-glucose 4-epimerase"] = ["UDPG4E"]
    beg_rxns_map["Triose-phosphate isomerase"] = ["TPI"]
    beg_rxns_map["steroid DELTA-isomerase"] = [""]
    beg_rxns_map["dodecenoyl-CoA isomerase"] = ["CTECOAI6","CTECOAI7","CTECOAI8"]
    beg_rxns_map["Glutamate-1-semialdehyde 2,1-aminomutase"] = [""]
    beg_rxns_map["Chalcone isomerase"] = [""]
    beg_rxns_map["Chloromuconate cycloisomerase"] = [""]
    beg_rxns_map["Tyrosine-tRNA ligase"] = [""]
    beg_rxns_map["Threonine-tRNA ligase"] = [""]
    beg_rxns_map["Isoleucine-tRNA ligase"] = [""]
    beg_rxns_map["Lysine-tRNA ligase"] = [""]
    beg_rxns_map["formate-tetrahydrofolate ligase"] = ["FTHFLi"]
    beg_rxns_map["Adenylosuccinate synthase"] = ["ADSS"]
    beg_rxns_map["DNA ligase (NAD+)"] = [""]
    return beg_rxns_map
end
