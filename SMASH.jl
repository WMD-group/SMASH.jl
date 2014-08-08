# S*M*A*S*H
# Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics.
#
# Jarvist Moore Frost, University of Bath
# File begun 2014-07-07

# Print titles...

SMASH={ {"System","Systematic","Sub","Simulated","Standard","Symbiotic"},
        {"Method","Mash","Martian","Metrication","Molecular","Mutual"},
        {"Analysis","Analytic","Ability","And","Atomic"},
        {"Subtype","Suitable","Sublime","Subtle"},
        {"Holonomy","Homeotype","Hypothetic"}}
TILT=   {"(+)","(-)","(0)"}

println("Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics")
print("S*M*A*S*H: ") 
for WORDS in SMASH
    print(  WORDS[1+rand(Uint32)%length(WORDS)]," ",
            TILT[1+rand(Uint32)%length(TILT)]," ")
end
println()

# Quite enough of that; let's get on with the real work... 
