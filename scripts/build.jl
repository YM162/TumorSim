using DrWatson
@quickactivate "TumorSim"

tumorsim_path = projectdir()
build_dir_path = projectdir("build")
binary_path = projectdir("build","TumorSim.so")
precompile_path = projectdir("test","runtests.jl")

rm(build_dir_path,force=true,recursive=true)
mkpath(build_dir_path)

using PackageCompiler
#using Pkg 

#cd(build_dir_path)

#Pkg.activate("")
#Pkg.add("Agents")
#Pkg.add("DrWatson")
#Pkg.add("DataFrames")
#Pkg.add(path=tumorsim_path)



PackageCompiler.create_sysimage(["TumorSim","Agents","DrWatson","DataFrames"]; sysimage_path=binary_path,precompile_execution_file=precompile_path)

#rm("Manifest.toml")
#rm("Project.toml")