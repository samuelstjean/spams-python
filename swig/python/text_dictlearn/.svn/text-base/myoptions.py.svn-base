import re, os, sys

global cur_doc
cur_doc = "";

def sphinxdoc_of_options(options):
    """
    Build a doc string in sphinx format from options list.
      - input : options = list of option definitions (see parse_args below)
      - output: tuple made of usage line and string of detailed options help
    """
    args = []
    argsdoc = []
    for (opt,nbargs,required,argname,comm) in options:
        if comm != "":
            if opt == None:
                argsdoc.append("\t* **%s**: %s" %(argname,comm))
            else:
                if (nbargs > 0):
                    argsdoc.append("\t* **%s** %s: %s" %(opt,argname,comm))
                else:
                    argsdoc.append("\t* **%s**: %s" %(opt,comm))
        if opt == None:
            s = "**%s**" %argname
        else:
            if (nbargs > 0):
                s = "**%s** *%s*" %(opt,argname)
            else:
                s = "**%s**" %opt
        if required:
            args.append(s)
        else:
            args.append("[%s]" %s)
    return (" ".join(args),"\n".join(argsdoc) + "\n")

def doc_of_options(options):
    """
    Make doc from options list.

      - input : options = list of option definitions (see parse_args below)
      - output:

        - txt : string of args list
        - doc : detailed help of options
    
    """
    (txt,doc) = sphinxdoc_of_options(options)
    rexp = re.compile('\*')
    return (re.sub(rexp,'',txt),re.sub(rexp,'',doc))

def prepare_doc(options,before = "\n",after= "",name = None):
    """
    Return a string made of options documentation with optional strings 
    prepended and appended.
    The resulting string is stored in a global variable for use by usage()
    """
    global cur_doc
    (txt,ldoc) = sphinxdoc_of_options(options)
    sname = None
    if name == '__main__' or name == None:
        sname = os.path.basename(sys.argv[0])
    if sname == '' or sname == None:
        sname = name
    if sname == None:
        doc = "%s\n" %txt
    else:
        doc = "**%s** %s\n" %(sname,txt)
    doc = doc + before + ldoc + after
    cur_doc = doc
    return doc

def usage(msg = None):
    """
    Prints the help message, that must have been built with prepare_doc
    """
    if(msg != None):
        print msg
    rexp = re.compile('\*')
    rexp2 = re.compile('\n\|')
    s = re.sub(rexp,'',cur_doc)
    s = re.sub(rexp2,'\n',s)
    print "Usage: ",  s
    sys.exit(1)
    
def parse_args(options,argv,usage = usage):
    """
    Process command line arguments (argv) according to a list of tuples (options):
     - the 1st element of a tuple is the option string,
     - the 2nd is the number of arguments of the option,
     - the 3rd indicates if option is mandatory,
     - the 4th is argname(s) and the last is a description string.

    Arguments starting with '-' are considered options.

    Returns a list of (option-strings,arg-value) and a list of non options arguments.
    """
    args = []
    opts = []
    h_opts = {}
    found_required = {}
    for (opt,nbargs,required,argname,comm) in options:
        if opt == None:
            continue
        h_opts[opt] = nbargs
        if required:
            found_required[opt] = False
    n = len(argv)
    i = 0
    while i < n:
        opt = argv[i]
        i += 1
        if opt[0] == '-':
            if opt in h_opts:
                if opt in found_required:
                    found_required[opt] = True
                nbargs = h_opts[opt]
                if(nbargs > 0):
                    i2 = i + nbargs - 1
                    if i2 >= n:
                        usage("Missing arg for option %s" %opt)
                    i2 += 1
                    if(nbargs == 1):
                        arg = argv[i]
                    else:
                        arg = argv[i:i2]
                    i = i2
                else:
                    arg = None
            else:
                usage('Bad option %s' %opt)
            opts.append((opt,arg))
        else:
            args.append(opt)
    l = []
    for opt in found_required.keys():
        if (not found_required[opt]):
            l.append(opt)
    if (len(l) != 0):
        usage("The following options are required : %s" %str(l))
    return opts, args

