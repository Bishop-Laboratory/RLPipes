# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build, check


def test_check():
  runner = CliRunner()
  buildres = runner.invoke(build, ['rseq_out/', 'test_data/samples.csv'])
  checkres = runner.invoke(check, ['rseq_out/'])
  assert checkres.exit_code == 0

