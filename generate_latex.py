import os
import argparse


def parse_list(txt):
    main_list = txt.split(",")
    run_id_list = []
    for sublist in main_list:
        if ":" in sublist:
            left = sublist.split(":")[0]
            right = sublist.split(":")[1]
            run_id_list.extend(range(int(left), int(right)+1))
        else:
            run_id_list.append(int(sublist))


    return run_id_list  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Generate latex reports")
    parser.add_argument('--run_range', type=str, action = 'store')
    parser.add_argument('--run_id',type=int, action ='store')
    parser.add_argument('--file',type=str,action = 'store')
    args = parser.parse_args()

    if args.run_id is not None:
        run_id_list = [args.run_id]

    elif args.run_range is not None:
        run_id_list = parse_list(args.run_range)

    elif args.file is not None:
        f = open(args.file)
        a = f.readline()
        a = a[:-1]
        
        run_id_list = parse_list(a)
        

    for run_id in run_id_list:
        print "Generating report #%s..." % run_id
        os.system("pdflatex %s_report.tex" % run_id)
        os.system("pdflatex %s_report.tex" % run_id)
