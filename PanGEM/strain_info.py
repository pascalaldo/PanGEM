import pandas as pd
import argparse
from pathlib import Path

def initialize_parser(parser):
    parser.description = "Generate strain info."
    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="File containing a list of filenames for gbk files to analyze.",
    )
    parser.add_argument(
        "--output1", "-o1",
        type=str,
        required=True,
        help="Output file that contains an overview of all strains (StrainInformation.csv).",
    )
    parser.add_argument(
        "--output2", "-o2",
        type=str,
        required=True,
        help="Strain info output, containing all strains and their source features.",
    )


def collect_strain_info(
    gbk_list_path,
    output_strains_path,
    output_info_path,
):
    with open(gbk_list_path, "r") as f:
        gbk_files = [Path(gbk_path.strip()) for gbk_path in f.readlines()]

    strains = pd.DataFrame({
        "strain": [p.stem for p in gbk_files],
        "Pathotype": list(range(len(gbk_files)))
    })
    strains["NCBI ID"] = strains["strains"]

    strain_info = []
    for gbk_path in gbk_files:
        with open(gbk_path, "r") as f:
            records = SeqIO.parse(f, "genbank")
            for record in records:
                for feature in record.features:
                    if feature.type == "source":
                        info = {}
                        info["file"] = gbk_path.name
                        info["id"] = gbk_path.stem
                        info["genome_size"] = len(record.seq) / 1000000
                        for q in feature.qualifiers.keys():
                            info[q] = "|".join(f.qualifiers[q])
                        strain_info.append(info)
    strain_info = pd.DataFrame(strain_info)
    strain_info.drop_duplicates(
        subset="id", keep="first", inplace=True, ignore_index=False
    )
    strain_info.to_csv(output_info_path)

def run(args):
    collect_strain_info(
        args.input,
        args.output1,
        args.output2,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
