function print_error(error)

    if !(error isa InterruptException)
        printstyled("Simulation failed...\n"; color=:red)
        rethrow(error)
    else
        printstyled("Simulation stopped\n"; color=:red)
    end

end