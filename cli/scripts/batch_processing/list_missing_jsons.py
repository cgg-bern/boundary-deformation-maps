import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Find .ovm files missing a corresponding result JSON file."
    )
    parser.add_argument("ovm_dir", type=Path, help="Directory containing .ovm files")
    parser.add_argument("json_dir", type=Path, help="Directory containing result JSON files")
    parser.add_argument(
        "suffix",
        type=str,
        help="Suffix appended to the .ovm stem to form the result filename "
             "(e.g. '_boundary_deformation_map_output.json')",
    )
    parser.add_argument("output_txt", type=Path, help="Path to write the list of missing names")
    args = parser.parse_args()

    ovm_files = {f.stem for f in args.ovm_dir.glob("*.ovm")}
    result_stems = {
        f.name[: -len(args.suffix)]
        for f in args.json_dir.glob(f"*{args.suffix}")
    }

    missing = sorted(ovm_files - result_stems)

    print(f"Total .ovm files: {len(ovm_files)}")
    print(f"Total result files: {len(result_stems)}")
    print(f"Missing results: {len(missing)}")

    args.output_txt.write_text("\n".join(missing) + ("\n" if missing else ""))
    print(f"Wrote missing list to {args.output_txt}")


if __name__ == "__main__":
    main()
