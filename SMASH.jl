# S*M*A*S*H
# Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics.
#
# Jarvist Moore Frost, University of Bath
# File begun 2014-07-07

#include("XDATCAR.jl") #XDATCAR reader...
push!(LOAD_PATH,"./")
using XDATCAR

# Packages to pull in...
#using Gaston    # interface to plot using GnuPlot
                # I would use the MatPlotlib hooks; but this actually works
#set_terminal("x11") # Rubbish installation on Jarv's Mac @ Work... - Keep it old school!

using Plots

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
t=Trajectory([],0,[],[])
f=open("testmd2-nonselective_XDATCAR","r")
XDATCAR.read(f,t)

println("OK; frames read...")
println("Fractional coordinates; t.frames[123]")
println(t.frames[123])

println("Expanding to real coordinates; t.frames[123]*t.cell")
println(t.frames[123]*t.cell)

println("Dividing to fractional coordinates; t.frames[123]/t.cell")
println(t.frames[123]/t.cell)

SUPER=2

Pb={}
for i in 1:123
    push!(Pb,mod((t.frames[i][65:72,:]*SUPER),1)*t.cell/SUPER)
end
println(Pb)

println("STATS")
println(mean(Pb))

plot(t.frames[123][:,1],t.frames[123][:,2]) #,"plotstyle","linespoints")

read(STDIN,Char) # wait for character before ending (and thus closing the GnuPlots)

end

println("Pb Extractor...")
Pb={} # ToDo: This code doesn't work :) FIXME 
v=hcat(t.frames[123],t.atomlookup) # pastes coords with atom #

for a in v[:]
    println(v," ") #,t.atomlookup[indexin(v,t.frames[123])])
    if (a[4]==82)
        push!(Pb,v)
    end
end

println(Pb)
#plot(Pb[:,1],Pb[:,2])

plot(t.frames[123][:,1],t.frames[123][:,2],"plotstyle","points")

read(STDIN,Char) # wait for character before ending (and thus closing the GnuPlots)

