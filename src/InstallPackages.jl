## InstallPackages.jl
# ------------------------------------------------------------------------------
# Installs all required packages automatically
# ---
# created by Sebastian Schmitt, 04.03.2022
# ------------------------------------------------------------------------------

import Pkg

function InstallPackages()
    global no_procs = 1
    all_installed = false
    while !all_installed
        try
            include("Init.jl")
            all_installed = true
        catch e
            if typeof(e) == LoadError
                pkg_name = split(e.error.msg,"\"")[2]
                Pkg.add(pkg_name)
                println("---------- Installed Package \"$(pkg_name)\"! ----------")
            else
                throw(e)
            end
        end
    end
end

InstallPackages()
