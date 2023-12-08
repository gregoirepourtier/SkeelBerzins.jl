using Test
using ExampleJuggler

ExampleJuggler.verbose!(true)

# Credit: https://github.com/j-fu/VoronoiFVM.jl/blob/master/test/runtests.jl
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method runtests() returning true
# or false depending on success.
function run_tests_from_directory(testdir, prefix)
    @info "Directory $(testdir):"
    examples = filter(ex -> length(ex) >= length(prefix) && ex[1:length(prefix)] == prefix,
                      basename.(readdir(testdir)))
    @info examples
    @testmodules(testdir, examples)
end

function run_all_tests()
    example_dir = joinpath(@__DIR__, "..", "examples")

    @testset "Test solve parameters" begin
        run_tests_from_directory(@__DIR__, "test_")
    end

    @testset "PDE Examples" begin
        run_tests_from_directory(example_dir, "Example1")
    end
    @testset "Interpolation Methods" begin
        run_tests_from_directory(example_dir, "Example2")
    end
    @testset "Integration SciML Ecosystem" begin
        run_tests_from_directory(example_dir, "Example3")
    end
end

run_all_tests()
