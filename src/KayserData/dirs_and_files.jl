# DIRS
const HAYSER_RAW_DATA_DIR = joinpath(RAW_DATA_DIR, "Hayser2013___data")
const HAYSER_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, basename(HAYSER_RAW_DATA_DIR))
mkpath(HAYSER_PROCESSED_DATA_DIR)    

# FILES
const HAYSER_CONV_TABLE1_FILE = joinpath(PROCESSED_DATA_DIR, "Hayser2013___table1_conv.bson")
const HAYSER_CONV_MEDIUM_FILE = joinpath(PROCESSED_DATA_DIR, "Hayser2013___medium_conv.bson")
