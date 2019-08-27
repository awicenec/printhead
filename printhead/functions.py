####
# All of this is for the interactive version...
###
__version__ = "5.0"
from printhead.printhead import FitsHead, HeadDict

def usage():
        """
        Prints out a short help.
        """
        msg = ("Script dumps headers of FITS files to stdout or creates header",
               "files. It supports compressed files (.gz and .Z)",
               "",
               "Synopsis: printhead.py [-s <KEYWORD> -H <number> -M <extnum> -S -x <type> -e -h]",
               " fname1 [fname2]...",
               "",
               "If only a file name is given, the primary header of that file is printed.",
               "",
               "--extract|-e:   All the headers of the files found are then",
               "                extracted to directories with the same name as the last",
               "                directory found in the path-names of the files. The",
               "                header files will have the same base name as the file, but",
               "                with the extension '.hdr'.",
               "--skey|-s:      if the given KEYWORD is found",
               "                in the header only the matching lines will be printed.",
               "",
               "--header|-H:    <number> specifies the number of the header to be printed.",
               "                If 99 is given, all headers are printed. If <number> is",
               "                negative only the structure of the file is printed.",
               "",
               "--xml|-x:       <type> is either 'vo' or 'xf'. All the headers of the files found",
               "                are then extracted to directories with the same name as the last",
               "                directory found in the path-names of the files. The",
               "                header files will have the same base name as the file, but",
               "                with the extension '.xml'. The files use XFits as a format",
               "                if 'xf' is specified and VOTable format if 'vo' is specified",
               "                for <type>.",
               "--Struct|-S     Show the structure of the FITS file",
               "",
               "--tsv|-t        Print keywords as a tab-separated-value list. The list contains",
               "                the items: keyword name, keyword value, comment, keyword type, index",
               "                The index item is the running number of the keyword within the header.",
               "",
               "--check|-c      If this flag is set the CRC32 checksum of the data part of the",
               "                extensions is calculated.",
               "--mode-m        [1]|0: If set to 0 the program does not try to skip the data part",
               "                between headers. This is useful for interpreting header files.",
               "--parse|-p      Switch the full parsing of the header on",
               "                extensions is calculated.",
               "--help|-h:      print this help and exit.",
               "",
               "Version: " + __version__)
        print('\n'.join(msg))


def run(args, skey='END', header=0, mode=1, struct=0, check=0):
        """
        Implements the loop around several files and opens either a
        pipe (compressed files) or the file directly.
        """
        for name in args:
          try:
            pH = FitsHead(name, skey=skey, show=header,
                          struct=struct, check=check, mode=mode)
            if skey != 'END':
                if header == 99:
                    heads = range(len(pH.HEAD))
                else:
                    heads = [header]
                for h in heads:
                    if list(pH.Extension[h]['index'].values()).count(skey) == 0:
                        print('%s\t%3d\t%s\t*not found*' % (name, h, skey))
                    else:
                       print("%s\t%3d\t%s\t%s" % (name, h, skey,
                                                  pH.Extension[h].getKeyword(skey)[1]))
            else:
                print(pH.HEAD[header])
          except Exception as e:
            pH = ''
            print(e)
#            sys.exit('<ERROR> unable to open file:' +name+' <ERROR>')
        return pH


def tsvFunc(args, skey='END', header=0, mode=1):
        """
        Implements the loop around several files and opens either a
        pipe (compressed files) or the file directly.

        INPUT:     string list, file name to process
                   string attribute skey, keyword to parse, default 'END', optional
                   int attribute header, >=0 number of header to return, default 0, optional
        OUTPUT:    tuple, (<FitsHead instance>, <list of tsv formatted lines>)
        """

        lines = []
        for name in args:
          try:
            pH = FitsHead(name, skey=skey, show=header, struct=1, mode=mode)
            tupleList = pH.parseFitsHead2TupleList(forceString=1)
            if header == 99:
                    hrange = range(len(tupleList))
            else:
                    hrange = [header]
            for hind in hrange:
                if skey != 'END':
                    if pH.Extension[hind].getKeyPos(skey) == -1:
                        lines += ['%s\t%s\t*not found*' % (name, skey)]
                    else:
                        ind = pH.Extension[hind].getKeyPos(skey)
                        # print skey, ind, tupleList[ind]
                        lines += ascii_load_lines([tupleList[hind]
                                                   [ind]], '\t', '\n')
                else:
                    lines += ascii_load_lines(tupleList[hind], '\t', '\n')
          except Exception as e:
              print(e)
              return
        return (pH, lines)


def ascii_load_lines(res, TABsep, RETsep):
    """
    Helper function takes a list of tuples, [(1,2,3,3,),(4,5,6,7)],
    where each tuple represents a record to be loaded, and
    returns a string formated according to the syntax used by the IQ load command.
    """
#    lines = RETsep.join(map(lambda x:TABsep.join(x),res))
    lines = map(lambda x: TABsep.join(x) + RETsep, res)

#    lines = ""
#    for row in res:
#       line =""
#       for column in row[:-1]:
#          # columns + Tab separator
#          line = line + str(column) + TABsep
#
#       # Columns + Last column + Enter Separator
#       line = line  + str(row[-1]) + RETsep
#       lines = lines + line
    return(lines)


def hdrExtract(name, xmlfl='', xtract=0, skey='END', show=0, struct=1, check=0, mode=1):
    """
    Extracts headers of all files found by glob(name) into
    header file <file_id>.hdr or <file_id>.xml. The last directory
    in the path defined by <name> is maintained also for the
    header files.
    """
    file_list = glob(name)
    if xmlfl >= 1:
        oext = '.xml'
    else:
        oext = '.hdr'

    if len(file_list) == 0:
        return -1
    for file in file_list:
        (path, base) = os.path.split(file)
        (fileb, ext) = os.path.splitext(base)
        if path:
            #last directory of orig-files will be used to order the
            #extracted headers

            night = os.path.split(path)[1]
        else:
            night = ''

        pH = FitsHead(file, skey=skey, show=show, struct=struct,
                      check=check, mode=mode)
        pH.fd.close()

        if ext == '.Z' or ext == '.gz':
            (file_id, ext) = os.path.splitext(fileb)
        else:
            file_id = fileb

        if night:
            if not os.path.isdir(night):
                    os.mkdir(night)
            ofnm = night + '/' + file_id + oext
        else:
            ofnm = file_id + oext
        if xtract == 1:
                #            print 'extracting header of file ',file,' to ',ofnm
            o = open(ofnm, 'w')
            o.write(pH.HEAD[0])
            o.close()
        elif xmlfl != '':
                #            print 'extracting header of file ',file,' to ',ofnm
            pH.parseFitsHead()
            XmlHead = pH.xmlHead(format=xmlfl, head=show)

            # if outfile is specified write the XML to it

            if len(ofnm) > 0:
                try:
                    o = open(ofnm, 'w')
                except:
                    print("ERROR: Unable to open ", outfile)
                    return 1

                for xml in XmlHead:
                    if type(xml) == type(''):
                        o.write(xml + "\n")
                    elif type(xml) == type([]):
                        o.write('\n'.join(xml))
                o.close()

        else:
            pH.parseFitsHead()

    fh = pH.Extension[0].Serialize()

    return pH


def mergeExtPrimary(file, extnum=1, outf=1, verb=1):
    """
    Merge Extension <extnum> (default 1) with primary header and attach the
    data of extension <extnum> as the primary data part.

    This is only possible if there is no original primary data part (NAXIS = 0)
    and if the data part of the extension is an image (XTENSION = 'IMAGE')
    """

    pk = {'SIMPLE': 0, 'XTENSION': 0, 'BITPIX': 1, 'NAXIS': 2, 'NAXIS1': 3,
          'NAXIS2': 4, 'NAXIS3': 5, 'NAXIS4': 6}

    pH = FitsHead(file, struct=1, show=99)
    pH.parseFitsHead()

    if pH.SIZE[0] != 0:
        print('There is a primary data part already! Size: ',
              pH.SIZE[0], ' bytes')
        print('Bailing out!')
        return pH

    print(len(pH.Extension))
    if pH.Extension[extnum].getKeyword('XTENSION')[1] != 'IMAGE':
        print('The extension is not an IMAGE but ',
              pH.Extension[extnum].getKeyword('XTENSION')[1])
        print('Bailing out!')
        return pH

    extHead = pH.Extension[extnum]

    maxind = pH.Extension[0]['index'].keys()[-1]
#    del(pH.Extension[0]['index'][maxind])   # get rid of the END keyword

    for k in list(extHead['index'].values())[:-1]:

        keyDict = extHead.getKeyDict(k)
        keyDict['index'] = {}
        if k[0:8] != 'XTENSION' and k[:3] != 'END':

            ind = pH.Extension[0].getKeyIndex(k)
            if ind == -1:
                ind = max(pH.Extension[0]['index'].keys())
                if k in pk:
                    ind = pk[k]
            keyDict['index'].update({ind: k})

            print(keyDict)
            pH.Extension[0].updateKeyword(keyDict, force=1)

#            pH.Extension[0]['nodes'].update({k:extHead['nodes'][k]})
#            print {k:extHead['nodes'][k]}
#        elif k[0:8] == 'HIERARCH':
#
#            hind = "['nodes']"
#            hkeys = k.split()
#
#            fullInd = ''
#            for hk in hkeys:
#                fullInd += "['"+hk+"']"
#
#            fitsLine = k
#            fitsLine = fitsLine + (29 - len(fitsLine)) * ' ' + '= '
#            for hk in hkeys:
#                if eval("'hk' in pH.Extension[0]"+hind):
#                    hind = hind + "['" + hk + "']"
#                    if eval("'Value' in extHead"+hind+"):
#                        vale = eval("extHead"+hind+"['Value']")
#                        valo = eval("pH.Extension[0]"+hind+"['Value']")
#                        if vale != valo:
#                            eval("pH.Extension[0]"+hind+".update({'Value':" + vale +"})")
#                            com = eval("extHead"+hind+"['Comment']")
#                            eval("pH.Extension[0]"+hind+".update({'Comment':'" + com +"'})")
#                        maxind = pH.Extension[0]['index'].keys()[-1]
#                        pH.Extension[0]['index'].update({maxind+1:k})
#                else:
#                    val = eval("extHead"+hind + "['" + hk + "']")
#                    eval("pH.Extension[0]"+hind+".update({hk:val})")
#                    del(val)
#                    hind = hind + "['" + hk + "']"

    maxind = pH.Extension[0]['index'].keys()[-1]
    pH.Extension[0]['index'].update({maxind+1: 'END'})

    if outf != 0:
        (path, base) = os.path.split(file)
        (fileb, ext) = os.path.splitext(base)
        outf = fileb + ".new" + ext
        outFd = open(outf, 'w')

    for dd in range(len(pH.Extension)):
        if dd != extnum:      # extnum header and data are with primary
                              # but keep other extensions.
            if (verb):
                    print("Extracting FITS header for extension number ", dd)
            if dd == 0:
                pH.Extension[dd].sortKeys()
                datapart = extnum
            else:
                datapart = dd
            FitsHd = pH.Extension[dd].Serialize()


# if outfile is specified write the FHead to it

        if outf != 0:
                #            print "Header cards:",len(FHead)
            totLen = 0
            for fitsLine in FitsHd:
                totLen += len(fitsLine)
                if type(fitsLine) == type(''):
                    outFd.write(fitsLine)
                elif type(fitsLine) == type([]):
                    for x in fitsLine:
                        outFd.write(x)

#           append the binary part to the header
            if (verb):
                    print("Extracting datapart number ", datapart)
            pH.fd.seek(pH.POS[datapart][0], 0)
            dataBlocks = int(ceil(pH.SIZE[datapart]/2880.))
            data = ''
            for d in range(dataBlocks):
                data = pH.fd.read(2880)
                outFd.write(data)
            del(data)

    if outf != 0:
            outFd.close()

    del(FitsHd)
    return pH
