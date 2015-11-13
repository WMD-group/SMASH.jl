# S*M*A*S*H
# Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics.
#
# Jarvist Moore Frost, University of Bath
# File begun 2014-07-07

include("XDATCAR.jl") #XDATCAR reader...

# Packages to pull in...
#using Gaston    # interface to plot using GnuPlot
                # I would use the MatPlotlib hooks; but this actually works
#set_terminal("x11") # Rubbish installation on Jarv's Mac @ Work... - Keep it old school!

# Print titles...
function print_titles()
    SMASH={ {"System","Systematic","Sub","Simulated","Standard","Symbiotic"},
            {"Method","Mash","Martian","Metrication","Molecular","Mutual"},
            {"Analysis","Analytic","Ability","And","Atomic","Aristotype"},
            {"Subtype","Suitable","Sublime","Subtle"},
            {"Holonomy","Homeotype","Hypothetic"}}
    TILT=   {"(+)","(-)","(0)"}

    println("Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics")
    print("S*M*A*S*H: ") 
    for WORDS in SMASH
        print(  WORDS[1+rand(Uint32)%length(WORDS)]," ",
                 TILT[1+rand(Uint32)%length( TILT)]," " )
    end
    println()
end
# Quite enough of that; let's get on with the real work... 

function plot_octahedra() #Doesn't work currently...
    figure(1)
    c = Gaston.CurveConf()
    c.plotstyle="pm3d"
    for i=1:4
        Xa,Ya,Za=0.0,0.0,0.0
        Xb=i%3
        Yb=i%2
        Zb=i
        gnuplot_send("set arrow $i from $Xa,$Ya,$Za to $Xb,$Yb,$Zb nohead") 
    end
    plot()
end


# MAIN
print_titles()
#plot_octahedra() #currently crashes out....

# Test routines...
t=Trajectory({},0,{},{})
f=open("testmd2-nonselective_XDATCAR","r")
read_XDATCAR(f,t)

println("OK; frames read...")
println("Fractional coordinates; t.frames[123]")
println(t.frames[123])

println("Expanding to real coordinates; t.frames[123]*t.cell")
println(t.frames[123]*t.cell)

println("Dividing to fractional coordinates; t.frames[123]/t.cell")
println(t.frames[123]/t.cell)

Pb={} # ToDo: This code doesn't work :) FIXME 
for (i,v) in enumerate(t.frames[123][:,1])
    println(i,v)
    if (t.atomlookup[i]==82)
        push!(Pb,v)
    end
end
#plot(Pb[:,1],Pb[:,2])

plot(t.frames[123][:,1],t.frames[123][:,2])

read(STDIN,Char) # wait for character before ending (and thus closing the GnuPlots)
