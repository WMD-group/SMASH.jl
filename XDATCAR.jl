# S*M*A*S*H
# Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics.
# VASP (electronic structure package) XDATCAR (ab-initio molecular dynamics file format) reader
#
# Jarvist Moore Frost, University of Bath
# File begun 2014-07-12

atomic={"H", "He", 
"Li", "Be", "B", "C", "N", "O", "F", "Ne", 
"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te","I", "Xe", 
"Cs","Ba",
"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
"Hf","Ta","W","Re","Os","Ir","Pt", "Au", "Hg", "Tl", "Pb", "Bi","Po","At","Rn","Fr","Ra"}
# indices are atomic number...

# Print our (Atomic Number) table
#for Z in 1:length(atomic)
#    println(Z," ",atomic[Z])
#end

type Trajectory #NB: need to read moar on constructors...
   cell
   natoms::Int
   frames
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

function read_XDATCAR(f::IOStream, t::Trajectory)
    l=readline(f) #Title
    l=readline(f) #Always a '1' ?
    
    t.cell=readdlm(IOBuffer(readnlines(f,3)))
    println(t.cell)

#    atomlookup=readdlm(IOBuffer(readnlines(f,2)))
    atomlookup=readmatrix(f,2)

    println(atomlookup)

    atoms=int(sum(atomlookup[2,1:end])) #quite ugly; but works

    println("$atoms atoms in XDATCAT frames")
    #frames=readdlm(f , dlm=(r"\r?\n?",r"Direct configuration=?"))
    #print(frames)

    nframe=0
    while !eof(f) 
        nframe=nframe+1

        stepsizeline=readline(f)
        push!(t.frames,readmatrix(f,atoms))
#        print(frame)
    end
    println("read_XDATCAR: $nframe Green Bottles...")
end


