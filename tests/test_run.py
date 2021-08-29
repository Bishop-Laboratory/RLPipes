# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build, check, run


def test_check():
  runner = CliRunner()
  buildres = runner.invoke(build, ['rseq_out_bams/', 'test_data/bam_test_samples_1.csv'])
  checkres = runner.invoke(check, ['rseq_out_bams/'])
  # runres = runner.invoke(run, ['rseq_out_bams/'])
  assert checkres.exit_code == 0
