cd(@__DIR__)
mkpath("build")
cd("build")

using PackageCompiler
using Pkg
Pkg.activate("")
Pkg.add("Agents")
Pkg.add("DataFrames")
Pkg.add(url="https://github.com/YM162/TumorSim.git")
PackageCompiler.create_sysimage(["TumorSim","Agents","DataFrames"]; sysimage_path="TumorSim.so",precompile_execution_file="../precompile.jl")

rm("Manifest.toml")
rm("Project.toml")
