#TODO: find julia native alternative for it
#=
    This script just run all the other scripts in the correct order.
=#

## ------------------------------------------------------------------------
using ArgParse

set = ArgParseSettings()
@add_arg_table! set begin
    "--force-pull"
        help = "make a sync with the remote"
        action = :store_true
    "--install"
        help = "run an installation script"
        action = :store_true

end

parsed_args = parse_args(set)
force_pull_flag = parsed_args["force-pull"]
install_flag = parsed_args["install"]
clear_cache_flag = parsed_args["clear-cache"]
clear_fva_models_flag = parsed_args["clear-fva-models"]

## ------------------------------------------------------------------------
using Pkg
try import DrWatson
catch
    import Pkg
    pkg"add DrWatson"
end
import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

## ------------------------------------------------------------------------
function force_pull()
    here = pwd()
    cd(@__DIR__)
    run(`git reset -q HEAD -- .`)
    run(` git checkout -- .`)
    try;
        run(`git pull`)
    catch err
        @warn string("fail to pull: error ", err)
    end
    cd(here)
end
force_pull_flag && force_pull()

## ------------------------------------------------------------------------
if install_flag
    # sync
    force_pull()
    # Install unregistered packages
    try
        pkg"rm Chemostat"
        pkg"rm UtilsJL"
    catch; end
    pkg"add https://github.com/josePereiro/UtilsJL.git#master"
    pkg"add https://github.com/josePereiro/Chemostat#adbeb2f"
    pkg"instantiate"
    pkg"build"
    pkg"test Chemostat"
end