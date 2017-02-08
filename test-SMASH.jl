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

    # Rewrote in a more sane order, via short note of W.Smith; and some googling
    d=b-a # b and a are fractional
    d=d-round(d) # Minimum image convention in fractional space
    d=unitcell*d # Project back to real space; i.e. Angstrom
    if verbose
        display(d)
    end
    d # returns in units of the unitcell; i.e. Angstrom
end

# This method from 
# In particular, reading J.Cumby's source code for PIEFACE:
# https://github.com/jcumby/PIEFACE/blob/372b6ff6166e4996d86084a3116a8b606c25acfa/pieface/ellipsoid.py#L54-L71
# Nb: Python/Numpy has crazy definiton of dotproduct = matrix multiplication. 8-[
#
# Khachiyan 1996 is pretty impenetrable:
# http://dx.doi.org/10.1006/jagm.1996.0062
# But these slides (see page 17) are fairly understandable:
# https://people.orie.cornell.edu/miketodd/ublndeta.pdf
#
# Cumby's Python implementation seems heavily influenced by Nina Moshtagh's Matlab impelementation:
# https://uk.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
# https://doi.org/10.1.1.116.7691
function minimumVolumeEllipsoid(points; MaxCycles=20, tolerance=1e-6)
    println("minimumVolumeEllipsoid( $points, MaxCycles=$MaxCycles")

    N=6 #hard coded for now
    D=3 # Number of dimensions? seems to affect step size; might be leading to numerical instability and solution collapse

    X=Array{Float64}(points)
    X=hcat(X,ones(X)) # I don't really understand what this is doing - padding Q with '1.0's
    println("X: ")
    display(X)

    err = 1e7
    u = zeros(N).+(1/N)
    display(u)

    while err>tolerance
        V=X*(diagm(u)*X')
        display(V)
        M=diag(X'*(inv(V)*X))
        display(M)

		j=indmax(M)
        maximum = M[j]
        println("\nindmax(M)=$j, maximum=M[j]=$maximum")

        delta = (maximum - D - 1) / ((D + 1) * (maximum - 1))
        println("\nDelta: $delta (should be >0)")
        new_u = (1.0 - delta) * u
        new_u[j] += delta
        display(new_u)
        err = norm(new_u - u)
        println("Err: $err")

		u=new_u
    end

end

"iterate over frames, calculate distnaces between Pb and I. Uses minimd PBCs!"
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

            minimumVolumeEllipsoid(octahedrapoints)

            @printf("\nPb-I6 'sumd' vector: \td=%f \t[%0.3f,%0.3f,%0.3f]",norm(sumd),sumd[1],sumd[2],sumd[3])
            grandsum+=sumd
            octahedra+=1
            #println()
        end
    end
    println("Grand sum: ",grandsum)
    println("Grandsum / number octahedra: ",grandsum/octahedra)
end

PbIdistance(t)


