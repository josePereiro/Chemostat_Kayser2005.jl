
function add_kayser_biomass(model::ChU.MetNet; UB = maximum(model.ub))

    model = ChU.expanded_model(model, size(model, 1) + 9, size(model, 2) + 10)

    ## ------------------------------------------------------------------
    # DNA
    # r88) 0.00247 ATP' + 0.00247 UTP + 0.00254 GTP + 0.00254 CTP + 0.015 ATP = DNA
    DNA = ChU.Met("DNA")
    DNA_RXN = ChU.Rxn("DNA_RXN"; 
        mets = ["atp_c", "utp_c", "gtp_c", "ctp_c", DNA.id, 
            "pi_c", "h_c"],
        S = [-(0.00247 + 0.015), -0.00247, -0.00254, -0.00254, 1.0, 
            (0.00247 + 0.015) + 0.00247 + 0.00254 + 0.00254,
            (0.00247 + 0.015) + 0.00247 + 0.00254 + 0.00254,
        ],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, DNA)
    ChU.set_rxn!(model, DNA_RXN)

    ## ------------------------------------------------------------------
    # One-carbon units and polyamine
    # r95) 0.00485 SER = C1
    C1 = ChU.Met("C1")
    C1_RXN = ChU.Rxn("C1_RXN"; 
        mets = ["ser_DASH_L_c", C1.id],
        S = [-0.00485, 1.0],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, C1)
    ChU.set_rxn!(model, C1_RXN)

    # r96) 0.0119 ATP + 0.01779 NADPH + 0.0119 GLUT = PA + 0.01397 AKG
    PA = ChU.Met("PA")
    PA_RXN = ChU.Rxn("PA_RXN"; 
        mets = ["atp_c", "nadph_c", "glu_DASH_L_c", PA.id, "akg_c", 
            "pi_c", "h_c"],
        S = [-0.0119, -0.01779, -0.0119, 1.0, 0.01397, 
            0.0119, 0.0119],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, PA)
    ChU.set_rxn!(model, PA_RXN)

    ## ------------------------------------------------------------------
    # r94) 0.0154 GLC6P + 0.0154 ATP = GLYC
    GLYC = ChU.Met("GLYC")
    GLYC_RXN = ChU.Rxn("GLYC_RXN"; 
        mets = ["g6p_c", "atp_c", GLYC.id, 
            "adp_c", "pi_c", "h_c"],
        S = [-0.0154, -0.0154, 1.0, 
            0.0154, 0.0154, 0.0154], 
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, GLYC)
    ChU.set_rxn!(model, GLYC_RXN)

    ## ------------------------------------------------------------------
    # r92) 0.00509 GLC6P + 0.0481 ATP + 0.0376 NADPH + 0.033 ACCOA + 
    # 0.0023 RIB5P + 0.0023 PEP + 0.00392 GLUT + 0.00235 G3P = 
    # LPS + 0.0023 NADH + 0.00392 AKG
    LPS = ChU.Met("LPS")
    LPS_RXN = ChU.Rxn("LPS_RXN"; # Lipopolysaccharider
        mets = [
            "g6p_c", "atp_c", "nadph_c", "coa_c", 
            "r5p_c", "pep_c", "glu_DASH_L_c", "g3p_c", 
            LPS.id, "nadh_c", "akg_c", 
            "adp_c", "pi_c", "h_c"
        ],
        S = [
            -0.00509, -0.0481, -0.0376, -0.033, 
            -0.0023, -0.0023, -0.00392, -0.00235, 
            1.0, 0.0023, 0.00392, 
            0.0481, 0.0481, 0.0481
        ],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, LPS)
    ChU.set_rxn!(model, LPS_RXN)

    ## ------------------------------------------------------------------
    # Peptidoglycan
    # r93) 0.00276 FRU6P + 0.0055 ACCOA + 0.00276 PEP + 0.00276 PYR + 0.00276 OAA + 
    # 0.02484 ATP + 0.0193 GLUT + 0.0193 NADPH = PG + 0.0138 AKG
    PG = ChU.Met("PG")
    PG_RXN = ChU.Rxn("PG_RXN";
        mets = [
            "f6p_c", "coa_c", "pep_c", "pyr_c", "oaa_c", 
            "atp_c", "glu_DASH_L_c", "nadph_c", PG.id, "akg_c", 
            "adp_c", "pi_c", "h_c"
        ],
        S = [
            -0.00276, -0.0055, -0.00276, -0.00276, -0.00276, 
            -0.02484, -0.0193, -0.0193, 1.0, 0.0138, 
            0.02484, 0.02484, 0.02484
        ],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, PG)
    ChU.set_rxn!(model, PG_RXN)

    ## ------------------------------------------------------------------
    # r89) 0.0129 PAL + 0.0129 OL + 0.0129 GAP + 0.0129 SER + 0.0258 ATP = LIPID
    # r90) 8 ACCOA + 7 ATP + 13 NADPH = PAL
    # r91) 9 ACCOA + 9 ATP + 15 NADPH = OL
    # -------------------------------------------------------------------
    # Total 
    # 0.0129 (8 ACCOA + 7 ATP + 13 NADPH) + 
    # 0.0129 (9 ACCOA + 9 ATP + 15 NADPH) + 
    # 0.0129 GAP + 0.0129 SER + 0.0258 ATP = LIPID
    # -------------------------------------------------------------------
    # (0.0129 * (8 + 9)) ACCOA + (0.0129 * (8 + 9)) ATP + 
    # (0.0129 * (13 + 15)) NADHP + 0.0129 GAP + 0.0129 SER + 0.0258 ATP = LIPID
    LIPID = ChU.Met("LIPID")
    LIPID_RXN = ChU.Rxn("LIPID_RXN"; 
        mets = ["coa_c", "atp_c", "nadph_c", "g3p_c", "ser_DASH_L_c", LIPID.id, 
            "adp_c", "pi_c", "h_c"],
        S = [-(0.0129 * (8 + 9)), -(0.0129 * (8 + 9)), -(0.0129 * (13 + 15)), -0.0129, -0.0129, 1.0, 
            (0.0129 * (8 + 9)), (0.0129 * (8 + 9)), (0.0129 * (8 + 9))],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, LIPID)
    ChU.set_rxn!(model, LIPID_RXN)

    ## ------------------------------------------------------------------
    # RNA
    # r87) 0.0165 ATP' + 0.0203 GTP + 0.0136 UTP + 0.0126 CTP + 0.0256 ATP = RNA
    RNA = ChU.Met("RNA")
    RNA_RXN = ChU.Rxn("RNA_RXN"; 
        mets = ["atp_c", "gtp_c", "utp_c", "ctp_c", RNA.id, 
            "pi_c", "h_c"],
        S =    [-(0.0165 + 0.0256), -0.0203, -0.0136, -0.0126, 1.0, 
            0.0165 + 0.0256, 0.0165 + 0.0256],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, RNA)
    ChU.set_rxn!(model, RNA_RXN)

    ## ------------------------------------------------------------------
    #  PROTEIN
    # 1.18 ATP  + 0.042 VAL + 0.005 TRP + 0.013 TYR + 0.024 THR + 
    # 0.021 SER + 0.021 PRO + 0.018 PHE + 0.015 MET + 0.033 LYS + 
    # 0.043 LEU + 0.028 ILE + 0.009 HIS + 0.059 GLY + 0.025 GLUM + 0.025 GLUT + 
    # 0.009 CYS + 0.023 ASP + 0.023 ASN + 0.028 ARG + 0.049 ALA = PROTEIN
    PROTEIN = ChU.Met("PROTEIN")
    PROTEIN_RXN = ChU.Rxn("PROTEIN_RXN"; 
        mets = [
            "atp_c"       , "val_DASH_L_c", "trp_DASH_L_c", "tyr_DASH_L_c", "thr_DASH_L_c", 
            "ser_DASH_L_c", "pro_DASH_L_c", "phe_DASH_L_c", "met_DASH_L_c", "lys_DASH_L_c",
            "leu_DASH_L_c", "ile_DASH_L_c", "his_DASH_L_c", "gly_c"       , "gln_DASH_L_c", "glu_DASH_L_c",
            "cys_DASH_L_c", "asp_DASH_L_c", "asn_DASH_L_c", "arg_DASH_L_c", "ala_DASH_L_c", PROTEIN.id,
            "adp_c", "pi_c", "h_c"],
        S = [
            -1.18 , -0.042, -0.005, -0.013, -0.024,
            -0.021, -0.021, -0.018, -0.015, -0.033,
            -0.043, -0.028, -0.009, -0.059, -0.025, -0.025,
            -0.009, -0.023, -0.023, -0.028, -0.049, 1.0, 
            1.18, 1.18, 1.18
        ],
        lb = 0.0, ub = UB
    )
    ChU.set_met!(model, PROTEIN)
    ChU.set_rxn!(model, PROTEIN_RXN)

    ## ------------------------------------------------------------------
    # Based on this simple stoichiometric model,the ATP consumption for 
    # maintenance requirements inthe absence of growthmATPwas determined as 
    # 2.81 mmol/ g h and the growth-associated energy consumption YX/ATP as 11.6 g / mol. 
    Y_X_ATP = 11.6

    ## ------------------------------------------------------------------
    # BIOMASS
    # r97) 1.12 PROTEIN + 0.56 RNA + LIPID + LPS + GLYC + PG + DNA + C1 + PA + 1/Y_X_ATP ATP = BIOMASS 
    # + 1/Y_X_ATP adp_c + 1/Y_X_ATP h_c + 1/Y_X_ATP pi_c + 1/Y_X_ATP ppi_c
    KAYSER_BIOMASS_RXN = ChU.Rxn(BIOMASS_IDER; 
        mets = [PROTEIN.id, RNA.id, LIPID.id, LPS.id, GLYC.id, PG.id, DNA.id, C1.id, PA.id, "atp_c", 
            "adp_c", "h_c", "pi_c"],
        S =    [ -1.12    , -0.56 , -1.0    , -1.0  , -1.0   , -1.0 , -1.0  , -1.0 ,  -1.0, -(1.0/Y_X_ATP) * 1e3,
            (1.0/Y_X_ATP) * 1e3, (1.0/Y_X_ATP) * 1e3, (1.0/Y_X_ATP) * 1e3],
        lb = 0.0, ub = 1000.
    );
    ChU.set_rxn!(model, KAYSER_BIOMASS_RXN)

    return model
end