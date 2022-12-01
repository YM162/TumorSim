#NOT WORKING RIGHT NOW
using DrWatson
@quickactivate "TumorSim"

using PackageCompiler
using Pkg 

tumorsim_path = projectdir()
build_dir_path = projectdir("build")
binary_path = projectdir("build","TumorSim.so")
precompile_path = projectdir("scripts","precompile_script.jl")

cd(build_dir_path)

Pkg.activate("")
Pkg.add("Agents")
Pkg.add("DrWatson")
Pkg.add("DataFrames")
Pkg.add(path=tumorsim_path)

mkpath(build_dir_path)
try rm(binary_path) catch end

PackageCompiler.create_sysimage(["TumorSim","Agents","DrWatson","DataFrames"]; sysimage_path=binary_path,precompile_execution_file=precompile_path)

rm("Manifest.toml")
rm("Project.toml")