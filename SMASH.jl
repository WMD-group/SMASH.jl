# S*M*A*S*H
# Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics.
#
# Jarvist Moore Frost, University of Bath
# File begun 2014-07-07

# Print titles...

SMASH={ {"System","Systematic","Sub","Simulated","Standard"},
    {"Method","Mash","Metrication"},
    {"Analysis","Analytic","Ability"},
    {"Subtype","Suitable",},
    {"Holonomy","Homeotype"}}
TILT={"+","-","0"}

println("Reconstituting the Glazer tilt notation for Perovskites from sampling molecular dynamics")
print("S*M*A*S*H: ") 
for i in SMASH
    word=i[1+rand(Uint32)%length(i)]
    print(word," ",TILT[1+rand(Uint32)%length(TILT)]," ")
end
println()

# Quite enough of that; let's get on with the real work... 
