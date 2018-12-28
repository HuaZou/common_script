#!/usr/bin/python3

#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This program abtains jobs status by user                  		         #
# Args :  2 args                                                             #
#   username:		username                                    			 #
#   Output: 		prefix              									 #
#                                                                            #
# Output: one result		                                                 #
#	result:         jobs status with MEN ect                     			 #
#----------------------------------------------------------------------------#

import sys
import os
import re
try:
    import argparse
except ImportError:
    sys.exit("Unable to find the argparse module")


def parse_arguments(args):
    parser = argparse.ArgumentParser(
        description="jobs status")
    parser.add_argument('-u', '--user', metavar='<username>',
        help="username in linux(centos)", required=True)
    parser.add_argument('-o', '--output', metavar='<prefix>',
        help="jobs status", required=True)
    return parser.parse_args()


def Status(username, output):
    command = "qstat -u " + str(username)
    qstat = os.popen(command, "r")
    head = qstat.readline
    lines = qstat.readline()
    out = open(output, "w")
    out.write("jobID\tqsubOrder\tusername\thard\tstatus\tinformation")
    print("jobID\tqsubOrder\tusername\thard\tstatus\tinformation")
    while lines:
        lines = qstat.readline()
        sample = lines.strip().split()
        if len(sample) != 0:
            command_jobs = "qstat -j " + str(sample[0])
            jobs = os.popen(command_jobs, "r")
            jobs_lines = jobs.readlines()
            for jobs_line in jobs_lines:
                jobs_sample = jobs_line.strip().split("\t")
                if jobs_sample[0] == "hard resource_list":
                    hard = jobs_sample[1]
                if jobs_sample[0] == "usage":
                    usage = jobs_sample[1]
                else:
                    usage = "None"
                res = str(sample[0]) + str(sample[0]) + str(sample[1]) + str(hard) + str(sample[4]) + str(usage)
                out.write("\n" + res)
                print("\n" + res)
            jobs.close()
    qstat.clost()
    out.close()


def main():
    args = parse_arguments(sys.argv)
    Status(args.user, args.output)


main()
