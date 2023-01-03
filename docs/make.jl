using Documenter
using SkeelBerzins

makedocs(
    sitename = "SkeelBerzins",
    format = Documenter.HTML(),
    modules = [SkeelBerzins]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
)
