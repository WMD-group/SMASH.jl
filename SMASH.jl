# S*M*A*S*H - module and XDATCAR Reader
# Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics.
# VASP (electronic structure package) XDATCAR (ab-initio molecular dynamics file format) reader
#
# Jarvist Moore Frost, University of Bath
# File begun 2014-07-12

module SMASH 

export atomic, Trajectory, read

atomic=["H", "He", 
"Li", "Be", "B", "C", "N", "O", "F", "Ne", 
"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te","I", "Xe", 
"Cs","Ba",
"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
"Hf","Ta","W","Re","Os","Ir","Pt", "Au", "Hg", "Tl", "Pb", "Bi","Po","At","Rn","Fr","Ra"]
# indices are atomic number...

# Print our (Atomic Number) table
#for Z in 1:length(atomic)
#    println(Z," ",atomic[Z])
#end

type Trajectory #NB: need to read moar on constructors...
   cell
   natoms::Int
   frames
   atomlookup::Array
end

function readnlines(f,n)
    local lines=""
    local i=1
    for i=1:n
        lines=lines*readline(f)
    end
    return (lines)
end

#readmatrix(f, nlines) = readdlm(IOBuffer(string([readline(f) for i in 1:nlines])))
readmatrix(f, nlines) = readdlm(IOBuffer(readnlines(f,nlines)))

# XDATCAR LOOKS LIKE THIS: 

# Perovskite_MD
#           1
#    12.580170    0.000078   -0.039648
#     0.000077   12.547782   -0.000152
#    -0.039643   -0.000153   12.594017
#   C    N    H    Pb   I
#   8   8  48   8  24
#Direct configuration=    50
#   0.43568603  0.50228894  0.53798652
#   0.03842720  0.49679247  0.48113604

function read_XDATCAR(f::IOStream) 
    t=Trajectory([],0,[],[]) 
    
    l=readline(f) #Title
    l=readline(f) #Always a '1' ?
    
    t.cell=readdlm(IOBuffer(readnlines(f,3))) #Unit cell spec (3x3 matrix of basis vectors)
    println(t.cell)

#    atomlookup=readdlm(IOBuffer(readnlines(f,2)))
    atomcrossref=readmatrix(f,2) # Ref to POTCAR; AtomName and #ofatoms

#    println(atomlookup)

    atoms=Int(sum(atomcrossref[2,1:end])) #Total atoms in supercell; quite ugly; but works
    t.natoms=atoms

    for atomtype = 1:length(atomcrossref[2,1:end]) # Each Atom... 
        for i in 1:atomcrossref[2,atomtype] # For i number of each atoms
            push!(t.atomlookup,indexin([atomcrossref[1,atomtype]],atomic)[1])
        end
    end
    print(t.atomlookup)

    println("$atoms atoms in XDATCAT frames")
    #frames=readdlm(f , dlm=(r"\r?\n?",r"Direct configuration=?"))
    #print(frames)
    
    nframe=0
    while !eof(f) 
        nframe=nframe+1

        stepsizeline=readline(f)
        push!(t.frames,readmatrix(f,atoms))   # Fractional coordinates!
#        print(frame)
    end
    println("read_XDATCAR: $nframe Green Bottles...")

    return t
end

# Print titles...
function print_titles()
    SMASH=Any[ Any["System","Systematic","Sub","Simulated","Standard","Symbiotic"],
            Any["Method","Mash","Martian","Metrication","Molecular","Mutual"],
            Any["Analysis","Analytic","Ability","And","Atomic","Aristotype"],
            Any["Subtype","Suitable","Sublime","Subtle"],
            Any["Holonomy","Homeotype","Hypothetic"]]
    TILT=   Any["(+)","(-)","(0)"]

    println("Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics")
    print("S*M*A*S*H: ")
    for WORDS in SMASH
        print(  WORDS[1+rand(Uint32)%length(WORDS)]," ",
                 TILT[1+rand(Uint32)%length( TILT)]," " )
    end
    println()
end
# Quite enough of that; let's get on with the real work...

print_titles()

end # Module
