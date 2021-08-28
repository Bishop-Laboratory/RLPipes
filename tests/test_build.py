# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build


def test_build():
  runner = CliRunner()
  result = runner.invoke(build, ['rseqOutdir/', 'test_data/samples.csv'])
  assert result.exit_code == 0

def test_build_fail_mode():
  runner = CliRunner()
  result = runner.invoke(build, ['rseqOutdir/', 'test_data/samples.csv', "-m", "bisDRIP"])
  assert result.exit_code == 2

def test_build_fail_mkdir():
  runner = CliRunner()
  result = runner.invoke(build, ['/asd/as/d/sad/asd/rseqOutdir/', 'test_data/samples.csv', "-m", "bisDRIP"])
  assert result.exit_code == 2

def test_build_bam():
  runner = CliRunner()
  buildres = runner.invoke(build, ['rseq_out_bams/', 'test_data/bam_test_samples_1.csv'])
  assert buildres.exit_code == 0

