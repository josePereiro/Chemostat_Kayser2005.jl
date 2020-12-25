import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using MAT

import Chemostat_Kayser2005: Chemostat
const Ch = Chemostat
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP

import Chemostat_Kayser2005: KayserData, iJR904, BegData, iJO1366
const iJR = iJR904
const iJO = iJO1366
const Kd = KayserData # experimental data
const Bd = BegData    # cost data

## ------------------------------------------------------------------
# LOAD RAW MODEL
src_file = iJR.MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["model"]

## ------------------------------------------------------------------
const MetNet = ChU.MetNet
const EMPTY_SPOT = ""
# This just prepare the metnet to hold more elements by making 
# all the numerical fields larger
function expanded_model(metnet::MetNet{T}, newM::Int, newN::Int) where T
    M, N = size(metnet)
    @assert all((newM, newN) .>= (M, N))

    function _similar_copy(col, fill, newdim)
        L = length(col)
        newcol = similar(col, newdim)
        newcol[L:end] .= fill
        newcol[1:L] .= col[1:L]
        return newcol
    end
    
    net = Dict()

    net[:S] = similar(metnet.S, newM, newN)
    @show findall(net[:S][762, :] .!= 0.0)

    net[:S][M:end, :] .= zero(T)
    net[:S][:, N:end] .= zero(T)
    net[:S][1:M, 1:N] .= metnet.S

    @show findall(net[:S][762, :] .!= 0.0)

    net[:b] = _similar_copy(metnet.b, 0, newM)
    net[:c] = _similar_copy(metnet.c, 0, newN)
    net[:lb] = _similar_copy(metnet.lb, 0, newN)
    net[:ub] = _similar_copy(metnet.ub, 0, newN)
    net[:rxns] = _similar_copy(metnet.rxns, EMPTY_SPOT, newN)
    net[:mets] = _similar_copy(metnet.mets, EMPTY_SPOT, newM)
    
    return MetNet(metnet; reshape = false, net...)
end

findempty(metnet::MetNet, field::Symbol) = findfirst(isequal(EMPTY_SPOT), getfield(metnet, field))
## ------------------------------------------------------------------
let 
    model = MetNet(mat_model; reshape=true)
    @show size(model)
    model = expanded_model(model, size(model, 1) + 1, size(model, 2) + 2)
    @show size(model)
    meti = findempty(model, :mets)
    @show meti
    @show rxnis = findall(model.S[meti, :] .!= 0.0)
    # model.S[rxnis, meti]
end;