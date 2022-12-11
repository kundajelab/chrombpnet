from weasyprint import HTML, CSS
import argparse

def main():

	parser = argparse.ArgumentParser(description='Convert html to pdf')
	parser.add_argument('-i','--input_html', required=True, type=str,  help='input file path to html')
	parser.add_argument('-o','--output_pdf', required=True, type=str,  help='output file path to pdf')
	args = parser.parse_args()

	css = CSS(string='''
		@page {
    		size: 1800mm 1300mm;
    		margin: 0in 0in 0in 0in;
		}
	''')
	HTML(args.input_html).write_pdf(args.output_pdf, stylesheets=[css])

if __name__=="__main__":
	main()
