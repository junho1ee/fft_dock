
## 
using BioStructures
using Bio3DView
using Blink

##
# Stored in the current working directory by default
downloadpdb("1EN2", dir="./data/pdb/")

##
# Read the PDB file
struc = read("./data/pdb/1EN2.pdb", PDBFormat)

##
# Print the structure
println(struc)

##
# visualize the structure
# https://biojulia.dev/BioStructures.jl/stable/documentation/
# viewpdb("1EN2")
viewstruc(struc)
