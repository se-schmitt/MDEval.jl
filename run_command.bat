@echo off
set /p var1= "Path of input file:   "
set var1=%var1:\=\\%
set /p var2= "Number of processors: "

julia MD_EVAL_Master.jl %var1% %var2%

@pause