using Documenter
using RMTool

makedocs(
    sitename = "RMTool",
    modules = [RMTool],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "1. Getting started with RMTool" => "gettingstarted.md",
        "2. Computing the spectral measure of random matrices" => "spectral.md",
        "Porting notes: MATLAB → Julia" => "porting.md",
        "Function reference" => "reference.md",
    ],
    warnonly = true,
)
