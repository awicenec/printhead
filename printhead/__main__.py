#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2012
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
import sys
import getopt
from printhead.functions import *

def main(args=sys.argv[1:]):
        opts, args = getopt.getopt(args, "s:H:x:M:m:peSctqh",
                                   ["parse", "extract", "skey=", "header=", "xml=", "struct", "merge=",
                                    "mode=", "check", "tsv", "quiet", "help"])
        _VERBOSE_ = 1

        xtract = 0
        parse = 0
        xmlfl = ''
        skeyfl = 0
        skey = 'END'
        show = -1
        struct = 0
        tsv = 0
        check = 0
        mergefl = 0
        hfl = 0
        breakfl = 0
        mode = 1

        while True:
            if len(args) == 0:
                usage()
                break
    #            sys.exit()
            try:

                for o, v in opts:
                    if o in ("-e", "--extract"):
                        xtract = 1
                    if o in ("-p", "--parse"):
                        parse = 1
                    if o in ("-s", "--skey"):
                        skey = v.strip()
                        skeyfl = 1
                    if o in ("-t", "--tsv"):
                        if hfl == 0:
                                show = 99
                        struct = 1
                        tsv = 1
                    if o in ("-H", "--header"):
                        hfl = 1
                        show = int(v)
                        struct = 1
                    if o in ("-x", "--xml"):
                        if v[:2].lower() == 'xf':
                            xmlfl = 'xfits'
                        elif v.lower() == 'vo':
                            xmlfl = 'vo'
                    if o in ("-S", "--Struct"):
                        show = -99
                        struct = 1
                    if o in ("-m", "--mode"):
                        mode = int(v)
                    if o in ("-M", "--merge"):
                        mergefl = 1
                        show = -1
                        struct = 1
                        breakfl = 1
                        print(">>> Careful this does not work correctly!!!!")
                        for f in args:
                            pH = mergeExtPrimary(f, extnum=int(v), verb=1)
                    if o in ("-q", "--quiet"):
                        usage()
                        breakfl = 1
                    if o in ("-c", "--check"):
                        # makes only sense with showing the structure.
                        show = -99
                        struct = 1
                        check = 1
                    if o in ("-h", "--help"):
                        usage()
                        breakfl = 1
            except Exception as e:
                errMsg = "Problem parsing command line options: %s" % str(e)
                print(errMsg)
                break
            try:
                if tsv == 1:
                    head = int(show)
                    if head < 0:
                            head = 0
                    (pH, lines) = tsvFunc(args, skey=skey, header=head, mode=mode)
                    for l in lines:
                        print(l[:-1])  # don't print the \n

                elif xtract == 1:
                    if xmlfl != '':
                        xtract = 0
                    for f in args:
                        pH = hdrExtract(f, xmlfl=xmlfl, show=show,
                                        xtract=xtract, mode=mode)
                elif skeyfl == 1:
                    for f in args:
                        head = int(show)
                        if head < 0:
                                head = 0
                        pH = run([f], skey=skey, header=head, mode=mode, struct=struct, check=check)
                elif xmlfl != '':
                    struct = 1
                    for f in args:
                        pH = FitsHead(f, skey=skey, show=show, struct=struct,
                                      check=check, mode=mode)
                        pH.fd.close()
                        pH.parseFitsHead()
                        XmlHead = pH.xmlHead(format=xmlfl, head=show)
                        for xml in XmlHead:
                            if type(xml) == type(''):
                                print(xml + "\n")
                            elif type(xml) == type([]):
                                print('\n'.join(xml))

                elif struct > 0:
                    if mergefl == 0:
                        for f in args:
                            pH = FitsHead(f, struct=struct, check=check, verbose=0,
                                          show=show, mode=mode)
                            if show == -99:
                                output = '\n'.join(pH.STRUCT)
                            elif show == 99:
                                output = ''.join(pH.HEAD)
                            elif show >= 0 and show <= len(pH.HEAD):
                                output = pH.HEAD[show]
                            else:
                                output = "Invalid header number specified. Should be: [0-%d,99]" % \
                                    (len(pH.HEAD)-1)
                            print(output)
                elif breakfl == 1:
                    break
                else:
                   pH = run(args)
                break
            except Exception as e:
               errMsg = "Problem extracting headers: %s" % str(e)
               print(errMsg)
               break


if __name__ == '__main__':

        args = []
        args.append(sys.argv[1:])
        main(args=args)
