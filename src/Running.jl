function run_model(filename)
    println(string("Running ",filename))
    pt=ENV["LUMERICAL_RUNNING_PATH"]
    run(`$pt $filename`)
    println(string("Done running ",filename))
end
