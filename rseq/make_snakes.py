import snakemake as snk
import sys
import json


def make_snakes(config_files):
    """Python wrapper for snakemake api for RSeq

    :param config_files: list of snakemake config files in .json format
    :return: None
    """
    for config_file in config_files:
        config = json.load(open(config_file))
        snk.snakemake('rseq/rseq.smk', unlock=True,
                      config=config)
        cores = config['cores'][0]
        ret_bool = snk.snakemake('rseq/rseq.smk', force_incomplete=True,
                                 config=config, cores=cores)
        print(ret_bool)
        if not ret_bool:
            break


if __name__ == "__main__":
    make_snakes(sys.argv[1:])
