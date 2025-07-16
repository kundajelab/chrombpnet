from weasyprint import HTML, CSS
import argparse

def main(input_html,output_pdf):
	css = CSS(string='''
		@page {
    		size: 1800mm 1300mm;
    		margin: 0in 0in 0in 0in;
		}
	''')
	HTML(input_html).write_pdf(output_pdf, stylesheets=[css])

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Convert html to pdf')
	parser.add_argument('-html','--input_html', required=True, type=str,  help='input file path to html')
	parser.add_argument('-pdf','--output_pdf', required=True, type=str,  help='output file path to pdf')
	args = parser.parse_args()

	main(args.input_html,args.output_pdf)
