
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import MAT
import Chemostat_Kayser2005: iJO1366, KayserData
const Kd = KayserData
const iJO = iJO1366
import Chemostat
const Ch = Chemostat
const ChU = Chemostat.Utils

##
mat_model = MAT.matread(iJO.MODEL_RAW_MAT_FILE)["iJO1366"]
model = ChU.MetNet(mat_model; reshape = true);

##


##
ChU.summary(model, 886)
