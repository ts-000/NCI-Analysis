$ErrorActionPreference = "Stop"
$DsspBin = if ($env:DSSP_BIN) { $env:DSSP_BIN } else { "mkdssp" }

python verify_install.py --dssp $DsspBin
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
python tests\selftest.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
python stage1_collect.py --output stage1_gold
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
python stage2_annotate.py `
  --stage1 stage1_gold `
  --output stage2_annotations `
  --dssp $DsspBin `
  --workers 2
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
python stage3_prepare.py `
  --stage2 stage2_annotations `
  --output stage3_fragments `
  --max-per-motif 0
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
