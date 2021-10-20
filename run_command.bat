@echo off
set /p var= "Path of input file: "
set var=%var:\=\\%

julia MD_EVAL_Master.jl %var%

@pause