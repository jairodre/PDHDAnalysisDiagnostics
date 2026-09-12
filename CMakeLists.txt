# PDHDAnalysisDiagnostics package scaffold.
#
# The source templates are intentionally neither built nor registered as art
# plugins. Add a target only after its InputTags, output schema, missing-product
# policy, and data/MC behavior have been implemented and reviewed.

add_subdirectory(DataProducts)
add_subdirectory(Alg)
add_subdirectory(Modules)
add_subdirectory(job)

install_source()
