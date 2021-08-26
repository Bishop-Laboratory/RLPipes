import click
import ast
from .run_workflow import make_snakes

# Allows passing strings to CLI and eval as python objects
# From https://stackoverflow.com/questions/47631914/how-to-pass-several-list-of-arguments-to-click-option
# Had to add the "str(value)" part or default would throw ValueErrors.
class PythonLiteralOption(click.Option):
  def type_cast_value(self, ctx, value):
    try:
      return ast.literal_eval(str(value))
    except ValueError:
      raise click.BadParameter(value)
          

# Help text
snakeHelp = """
      A dctionary of arguments to pass to the snakemake python API. Default: "{'use_conda': True}".
      Please visit the snakemake API reference for the full list of available options.
"""


# Get the shared options
# From https://stackoverflow.com/questions/40182157/shared-options-and-flags-between-commands
verify_run_options = [
    click.option("--smargs", "-s", cls=PythonLiteralOption, help=snakeHelp, default="{'use_conda': True}")
]


# Function for addint these options to the click command
def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options


@click.group()
def cli(**kwargs):
  """
  RSeq: An R-loop mapping pipeline with built-in QC.
  """
  pass


@cli.command("build")
@click.option("--genome", "-g", help="UCSC genome for samples (e.g., 'hg38')")
def build(genome):
    """Configures an RSeq workflow."""
    click.echo(f"{genome}")
    

@cli.command("check")
@click.argument("config")
@add_options(verify_run_options)
def check(config, **kwargs):
    """
    Validates an RSeq workflow.
    
    CONFIG: required. Filepath to config.json generated by RSeqCLI build.
    """
    click.echo(f"{config}")
    print(kwargs)
    smargs=kwargs['smargs']
    click.echo("smargs, type: {}  value: {}".format(
        type(smargs), smargs))
    make_snakes(config, smargs, verify=True)


@cli.command("run")
@click.argument("config")
@add_options(verify_run_options)
def run(config, **kwargs):
    """
    Executes an RSeq workflow.
    
    CONFIG: required. Filepath to config.json generated by RSeqCLI build.
    """
    click.echo(f"{config}")
    print(kwargs)
    smargs=kwargs['smargs']
    click.echo("smargs, type: {}  value: {}".format(
        type(smargs), smargs))

