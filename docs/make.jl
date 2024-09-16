using TestPackage
using Documenter

DocMeta.setdocmeta!(TestPackage, :DocTestSetup, :(using TestPackage); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
  file for
  file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [TestPackage],
    authors = "Paul Schrimpf <paul.schrimpf@ubc.ca>",
    repo = "https://github.com/ubcecon/TestPackage.jl/blob/{commit}{path}#{line}",
    sitename = "TestPackage.jl",
    format = Documenter.HTML(; canonical = "https://ubcecon.github.io/TestPackage.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/ubcecon/TestPackage.jl")
