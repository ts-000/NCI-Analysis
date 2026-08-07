$ErrorActionPreference = "Stop"
$DsspBin = if ($env:DSSP_BIN) { $env:DSSP_BIN } else { "mkdssp" }

python verify_install.py --dssp $DsspBin
python tests\selftest.py
python stage1_collect.py --output stage1_gold
python stage2_annotate.py `
  --stage1 stage1_gold `
  --output stage2_annotations `
  --dssp $DsspBin `
  --workers 2
python stage3_prepare.py `
  --stage2 stage2_annotations `
  --output stage3_fragments `
  --max-per-motif 0
