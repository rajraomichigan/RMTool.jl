# test/readme_examples.jl
#
# Extracts every ```julia fenced block from README.md, writes them to a sandbox
# file, and runs them, so the examples shown in the README cannot drift from
# what the package actually does. Blocks that are illustrative rather than
# runnable can be skipped by fencing them as ```julia-norun instead of ```julia.

using Test

@testset "README examples run" begin
    readme = joinpath(@__DIR__, "..", "README.md")
    if !isfile(readme)
        @info "README.md not found next to the package; skipping example test."
        return
    end

    text = read(readme, String)
    # capture ```julia ... ``` blocks (but not ```julia-norun ... ```)
    blocks = String[]
    for mblock in eachmatch(r"```julia\n(.*?)```"s, text)
        push!(blocks, mblock.captures[1])
    end
    @test !isempty(blocks)   # the README should contain at least one runnable example

    mktempdir() do dir
        script = joinpath(dir, "readme_block.jl")
        for (i, code) in enumerate(blocks)
            # Skip blocks that are REPL transcripts (lines beginning with "julia>")
            # or that reference Plots, which is not a dependency of the test env.
            occursin(r"^\s*julia>", code)  && continue
            occursin("Plots", code)        && continue
            occursin("using Plots", code)  && continue
            write(script, "using RMTool\nusing LinearAlgebra\n" * code)
            @testset "README block $i" begin
                @test (include(script); true)
            end
        end
    end
end
