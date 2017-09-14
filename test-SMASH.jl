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

"""
 minimumVolumeEllipsoid(points; tolerance=1e-3, verbose::Bool=true)

Calculates the minimum volume ellipsoid for the submitted point cloud.
Dimension of problem and Number of points discovered by inspection of
point cloud.

This method uses the Khachiyan algorithm, solving the dual problem in N+1 dimensions.

Khachiyan 1996 is pretty impenetrable:
http://dx.doi.org/10.1006/jagm.1996.0062
But these slides (see page 17) are fairly understandable:
https://people.orie.cornell.edu/miketodd/ublndeta.pdf
   - here I use variable names as in the above talk.

It was implemented with considerable assistance by reading J.Cumby's source code for PIEFACE:
https://github.com/jcumby/PIEFACE/blob/372b6ff6166e4996d86084a3116a8b606c25acfa/pieface/ellipsoid.py#L54-L71
Nb: Python/Numpy has crazy definiton of dotproduct = matrix multiplication. 8-[

Cumby's Python implementation seems heavily influenced by Nina Moshtagh's Matlab implementation:
https://uk.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
https://doi.org/10.1.1.116.7691
   - My implementation follows Moshtagh closely, as I find the code  well documented & the Matlab syntax is very close to Julia.
"""
function minimumVolumeEllipsoid(points; tolerance=1e-3, verbose::Bool=true)
    # N - number of points; D - dimension of problem, by inspection of point cloud
    (N,D)=size(points)

    X=Array{Float64}(points) # Forces into Float64 repr (from Any)
    X=hcat(X,ones(N))' # Lift D-dimension vectors into a higher dimensional space; i.e. pad each D-dimension tuple with 1.0
    # And then transpose to match the expected ordering for the linear algebra below
    if verbose
        println("Padded set of vectors lifted into higher dimension (X): ")
        display(X)
    end

    err = 1e7
    u = zeros(N).+(1/N) # Starts with u=eye(N)/N ; the uniform distribution

    count=1

    while err>tolerance
        V=X*(diagm(u)*X')
        #display(V)
        M=diag(X'*(inv(V)*X))
        #display(M)

		(maximum,j)=findmax(M)
        delta = (maximum - D - 1) / ((D + 1) * (maximum - 1))
        
        new_u = (1.0 - delta) * u
        new_u[j] += delta
        err = norm(new_u - u)
         
        count=count+1
		u=new_u

        if verbose
            println("indmax(M)=$j, maximum=M[j]=$maximum")
            println("Delta: $delta (should be >0)")
            println("Err: $err")
            println("Loops: $count")
        end 
    end

    if verbose
        println("CONVERGED TO MINIMUM VOLUME")
    end
	
# Comment from Moshtagh's code:
# %%%%%%%%%%%%%%%%%%% Computing the Ellipse parameters%%%%%%%%%%%%%%%%%%%%%%
# % Finds the ellipse equation in the 'center form': 
# % (x-c)' * A * (x-c) = 1
# % It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
# % of the ellipse. 
    P=points'
    U=diagm(u)

	centre = P*u
	A = 1/D * inv(P*U*P' - (P*u)*(P*u)')
	U, s, rotation = svd(A)

    if verbose
        println("Centre of ellipse: $centre")
        println("The 'A' matrix for ellipse, centre-form (x-c)'*A*(x-c)=1:")
        display(A)
        println("\nSVD decomposition, svd(A):")

        println("U:") 
        display(U)
        println("\nSigma:")
        display(s)
        println("\nRotation matrix:")
        display(rotation)
    end

    radii = 1.0./sqrt(s)

    vol=4/3*pi*prod(radii)
    sphereradius=prod(radii)^(1/3)
    
    # Ellipsoid shape measure r3/r2 - r2/r1
    # Nb: As haven't crosschecked SVD definition of radii + comparison to numpy, could be back to front here
    shapeparam=radii[3]/radii[2] - radii[2]/radii[1]

    if verbose
        println("\nradii: $radii")
        println("Volume ellipsoid: $vol Equivalent Sphere radius: $sphereradius")
        println("Ellipsoid shape measure r3/r2 - r2/r1 (Care! definitons.) $shapeparam")
    end

    return(shapeparam)
end

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


