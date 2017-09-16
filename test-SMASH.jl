push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using SMASH

# Test routines...
f=open("testmd2-nonselective_XDATCAR","r")
t=SMASH.read_XDATCAR(f) #Returns type XDATCAR.Trajcetory

"Iterate over frames, calculate distances between Pb and I. Uses minimd PBCs!"
function PbIdistance(t)
    grandsum=0.0
    octahedra=0

    for i in 1:t.nframes
        Pbs=t.frames[i][t.atomlookup .=="Pb",:]
        Is=t.frames[i][t.atomlookup .=="I",:]

        for j=1:size(Pbs,1)
            Pb=Pbs[j,:]
            #display(Pb)
#            @printf("\nPb %d: at [%f,%f,%f] Fractional ",j,Pb[1],Pb[2],Pb[3])

            sumd=0.0
            octahedrapoints=Matrix(0,3)

            for k=1:size(Is,1)
                I=Is[k,:]
                #display(I)

                d=minimd(Pb,I,t.cell)
                #println("Pb: ",Pb," I: ",I," Diff: ",Pb-I, " Norm: ",norm(Pb-I))
                if (norm(d)<4) # MAGIC NUMBER; Pb-I distance angstroms
                    #@printf(" %0.3f",norm(minimd(Pb,I,t.cell)))
#                    @printf("\n I %d at \td=%0.3f \t[%0.3f,%0.3f,%0.3f,]",k,norm(d),d[1],d[2],d[3] )
                    octahedrapoints=[octahedrapoints; d']
                    sumd+=d
                end
            end

            shapeparam=minimumVolumeEllipsoid(octahedrapoints,verbose=false)
            @printf("Minimum volume ellipsoid shape param: %f \n",shapeparam)
            @printf("Pb-I6 'sumd' vector: \td=%f \t[%0.3f,%0.3f,%0.3f] \n",norm(sumd),sumd[1],sumd[2],sumd[3])

            grandsum+=sumd
            octahedra+=1
            #println()
        end
    end
    println("Grand sum: ",grandsum)
    println("Grandsum / number octahedra: ",grandsum/octahedra)
end

PbIdistance(t)


