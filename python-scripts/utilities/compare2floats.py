import argparse

def compareTwoFloats(number1,number2):
	if float(number1) > float(number2):
		print True
	elif float(number1) < float(number2):
		print False
	else:	
		print "Numbers provided are the same!"

parser=argparse.ArgumentParser(description='This is a simple script to compare two floats and return true if first number is higher than the second.')
parser.add_argument(dest='number1', help='first number to compare', type=float)
parser.add_argument(dest='number2', help='second number to compare', type=float)

args=parser.parse_args()

compareTwoFloats(args.number1,args.number2)
