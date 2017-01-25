push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using SMASH

# Test routines...
f=open("testmd2-nonselective_XDATCAR","r")
t=SMASH.read_XDATCAR(f) #Returns type XDATCAR.Trajcetory

"minimimum distance image convention; using a unitcell=3x3, with a and b vectors in fractional coords"
function minimd(a, b, unitcell; verbose::Bool=false)
    
    if verbose # pretty-print the call signature
        println("minimd(a,b,unitcell) with variables: \n")
        display(a)
        display(b)
        display(unitcell)
    end

    d=unitcell.*(b-a) # Nb: changed in Julia 0.5, need to be explicit about per-element broadcast (.)
    d=d-round(d/unitcell)*unitcell
    d # returns in units of the unitcell; i.e. Angstrom
end


"iterate over frames, calculate distnaces between Pb and I. Uses minimd PBCs!"
function PbIdistance(t)
    for i in 1
        Pbs=t.frames[i][t.atomlookup .=="Pb",:]
        Is=t.frames[i][t.atomlookup .=="I",:]

        for j=1:size(Pbs,1)
            Pb=Pbs[j,:]
            #display(Pb)
            print("\nPb: at ",Pb," Fractional ")
            for k=1:size(Is,1)
                I=Is[k,:]
                #display(I)

                #println("Pb: ",Pb," I: ",I," Diff: ",Pb-I, " Norm: ",norm(Pb-I))
                if (norm(minimd(Pb,I,t.cell))<4) # MAGIC NUMBER; Pb-I distance angstroms
                    #@printf(" %0.3f",norm(minimd(Pb,I,t.cell)))
                    myd=minimd(Pb,I,t.cell)
                    @printf("\n d=%0.3f\n\t%0.3f x %0.3f y %0.3f z",norm(myd),myd[1],myd[2],myd[3] )
                end
            end
            println()
        end
    end
end

PbIdistance(t)
