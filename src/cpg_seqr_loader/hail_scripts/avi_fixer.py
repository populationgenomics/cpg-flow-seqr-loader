import argparse
import hail as hl


def main():
    parser = argparse.ArgumentParser(description='Fix a tuple-wrapped row field in a Hail MatrixTable')
    parser.add_argument('--input', '-i', required=True, help='Path to the input MatrixTable (.mt)')
    parser.add_argument('--output', '-o', required=True, help='Path to write the fixed MatrixTable (.mt)')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite the output path if it exists')

    args = parser.parse_args()

    print(f'Reading MatrixTable from {args.input}...')
    mt = hl.read_matrix_table(args.input)

    # The Fix: Access the first element of the tuple to get the raw value back
    # We use annotate_rows to overwrite the existing 'avis' field
    print("Unwrapping the 'avis' field...")
    mt = mt.annotate_rows(avis=mt.avis[0])

    print(f'Writing fixed MatrixTable to {args.output}...')
    mt.write(args.output, overwrite=args.overwrite)

    print("Done! Field 'avis' has been flattened.")


if __name__ == '__main__':
    main()
