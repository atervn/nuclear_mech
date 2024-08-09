
include("setup_functions.jl")

if !(@isdefined envelopeType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end

ipar = inputParametersType();
ipar = read_parameters(ipar,"parameters.txt")

chroStiff = [10e-4 14e-4 18e-4]
lamStiff = [10e-4 13e-4 16e-4];
areaComp = [0.1e3 0.5e3 1e3 2e3]
bulk = [1000 5000 10000 20000]

ind = 1
for i1 = lamStiff
    for i2 = areaComp
        for i3 = bulk
            for i4 = chroStiff
                ipar.laminaStiffness = i1;
                ipar.areaCompressionStiffness = i2;
                ipar.bulkModulus = i3;
                ipar.chromatinStiffness = i4;

                f = open(".\\parameters\\parameters_"*string(ind)*".txt", "w")

                fields = propertynames(ipar)

                for i = eachindex(fields)

                    write(f,string(fields[i]) *","* string(getfield(ipar,fields[i])) * "\n" )

                end

                close(f)
                global ind += 1
            end
        end
    end
end
