# DIRS
const KAYSER_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "Hayser2013___data")
const KAYSER_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(KAYSER_RAW_DATA_DIR))
mkpath(KAYSER_PROCESSED_DATA_DIR)    

# FILES
const KAYSER_CONV_TABLE1_FILE = joinpath(PROCESSED_DATA_DIR, "Hayser2013___table1_conv.bson")
const KAYSER_CONV_MEDIUM_FILE = joinpath(PROCESSED_DATA_DIR, "Hayser2013___medium_conv.bson")
