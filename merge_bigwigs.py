import pyBigWig
import argparse

parser=argparse.ArgumentParser(description="Normalize bigwigs")
parser.add_argument("-bigwig","--input_bigwig", help="input bigwig to normalize")
parser.add_argument("-o","--output_path", help="output path")
parser.add_argument("-s","--scale", type=float,  help="scale")
args = parser.parse_args()
