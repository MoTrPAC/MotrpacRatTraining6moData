# ---- resave the .rds files in data/ with compressed format ----
#
# details:
#   The type of compression is chosen automatically by tools
#   This command need only be run once and may take several minutes.

tools::resaveRdaFiles(paths = 'data')
