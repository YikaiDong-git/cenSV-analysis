import os
import csv
import json
import argparse

def main():
    parser = argparse.ArgumentParser(description='Generate CSV table from benchmark data')
    parser.add_argument('--input-json', required=True)
    parser.add_argument('--output-csv', required=True)
    parser.add_argument('--output-html', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.input_json):
        print("Input JSON not found")
        return

    with open(args.input_json, 'r') as f:
        data = json.load(f)

    details = data.get('details', [])
    if details:
        with open(args.output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['comparison','chrom','svtype','ac_ref','ac_tgt','svlen_ref','svlen_tgt','diff_type'])
            writer.writeheader()
            for d in details:
                row = {
                    'comparison': f"Q{d.get('mapq')}_L{d.get('minsvlen')}",
                    'chrom': d.get('chrom'),
                    'svtype': d.get('svtype'),
                    'ac_ref': d.get('ac_ref',''),
                    'ac_tgt': d.get('ac_tgt',''),
                    'svlen_ref': d.get('svlen_ref',''),
                    'svlen_tgt': d.get('svlen_tgt',''),
                    'diff_type': ','.join(d.get('diff_cols')) if d.get('diff_cols') else d.get('type')
                }
                writer.writerow(row)

    html = "<html><head><script src=\"https://cdn.plot.ly/plotly-2.27.0.min.js\"></script></head><body><h1>Diff visualization (open CSV for details)</h1><p>CSV saved at: %s</p></body></html>" % args.output_csv
    with open(args.output_html,'w') as f:
        f.write(html)

    print('Saved CSV and minimal HTML at', args.output_csv, args.output_html)

if __name__ == '__main__':
    main()
