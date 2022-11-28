#!/usr/bin/python 

import sys, re, os
import hashlib
import argparse as ap 
from collections import Counter


def parse_argument(args):
	parser = ap.ArgumentParser(description='transfer data from source directory to destination directory')
	parser.add_argument(
		'-s', '--source', metavar='<source directory>', type=str,
		help='directory for source files', required=True)
	parser.add_argument(
		'-d', '--dest', metavar='<destination directory>', type=str,
		 help='directory for destination files', required=True)
	parser.add_argument(
		'-n', '--suffix', metavar='<suffix>', type=str,
		 help='the prefix of source files', required=True)		 
	parser.add_argument(
		'-l', '--log', metavar='<log file>', type=str,
		help='log file for success', required=True)

	return parser.parse_args()

 
def md5sum(filename):
	
	hash_md5 = hashlib.md5()
	with open(filename, "rb") as f:
		for chunck in iter(lambda: f.read(4096), b""):
			hash_md5.update(chunck)
	return hash_md5.hexdigest()


def sync_file(file_name, file_path, dest_path):
	
	if os.path.isfile(file_path):
		dest_file = os.path.join(dest_path, file_name)
		if not os.path.isfile(dest_file):
			sync_cmd = "rsync -atvP  " + file_path + "  "+ dest_path
			os.system(sync_cmd)

	return dest_file
	

def check_sum(sync, dest):
	
	sync_md5 = md5sum(sync)
	dest_md5 = md5sum(dest)

	if Counter(sync_md5) == Counter(dest_md5):
		log = "Successed"
	else:
		log = "Failed"

	return log


def recursion(file_name, file_path, dest):
	sync_res = sync_file(file_name, file_path, dest)
	if os.path.isfile(sync_res):
		log_res = check_sum(file_path, sync_res)
	if 	log_res != "Successed":
		return recursion(file_name, file_path, dest)
	else:
		return log_res


def main():
	args = parse_argument(sys.argv)
	sfx = "." + args.suffix

	log_f = open(args.log, "w")
	for subdir, dirs, files in os.walk(args.source):
		for file in files:
			filepath = subdir + os.sep + file
			if filepath.endswith(sfx):
				log_result = recursion(file, filepath, args.dest)
				log_f.write(file + "\t" + log_result + "\n")

	log_f.close()

if __name__ == "__main__":
	main()
