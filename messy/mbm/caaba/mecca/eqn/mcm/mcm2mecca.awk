# Authors: Sebastian Tauer, Hartwig Harder, Rolf Sander (MPI Mainz, 2017)

# This script converts the input file extracted from MCM in mecca format,
# by separating species and equations and changing name of photolysis rates

# Example usage:
# gawk -f mcm2mecca.awk -v unkfile=tmp_unknown.kpp -v jsubstfile=Jsubstitution.txt -v spcfile=tmp_allnew.spc -v eqnfile=$eqnfile -v logfile=$logfile mcmfile.kpp
# Normally, mcm2mecca.awk is called via xmcm2mecca and not directly.

BEGIN {
    print "executing mcm2mecca.awk" >> logfile
    dontedit = "created automatically by mcm2mecca.awk, do not edit!"
    printf "// %s\n", dontedit > eqnfile
    # declare MCM constants:
    print "#INLINE F90_GLOBAL"                                      >> eqnfile
    print "REAL(dp) :: &"                                           >> eqnfile
    print "  KRO2NO, KRO2HO2, KAPHO2, KAPNO, KRO2NO3, KNO3AL, &"    >> eqnfile
    print "  KDEC, KROPRIM, KROSEC, KCH3O2, K298CH3O2, K14ISOM1, &" >> eqnfile
    print "  KD0, KDI, KRD, FCD, NCD, FD, KBPAN, KC0, KCI, KRC, &"  >> eqnfile
    print "  FCC, NC, FC, KFPAN, K10, K1I, KR1, FC1, NC1, F1, &"    >> eqnfile
    print "  KMT01, K20, K2I, KR2, FC2, NC2, F2, KMT02, K30, &"     >> eqnfile
    print "  K3I, KR3, FC3, NC3, F3, KMT03, K40, K4I, KR4, FC4, &"  >> eqnfile
    print "  NC4, F4, KMT04, KMT05, KMT06, K70, K7I, KR7, FC7, &"   >> eqnfile
    print "  NC7, F7, KMT07, K80, K8I, KR8, FC8, NC8, F8, KMT08, &" >> eqnfile
    print "  K90, K9I, KR9, FC9, NC9, F9, KMT09, K100, K10I, &"     >> eqnfile
    print "  KR10, FC10, NC10, F10, KMT10, K1, K3, K4, K2, &"       >> eqnfile
    print "  KMT11, K120, K12I, KR12, FC12, NC12, F12, KMT12, &"    >> eqnfile
    print "  K130, K13I, KR13, FC13, NC13, F13, KMT13, K140, &"     >> eqnfile
    print "  K14I, KR14, FC14, NC14, F14, KMT14, K150, K15I, &"     >> eqnfile
    print "  KR15, FC15, NC15, F15, KMT15, K160, K16I, KR16, &"     >> eqnfile
    print "  FC16, NC16, F16, KMT16, K170, K17I, KR17, FC17, &"     >> eqnfile
    print "  NC17, F17, KMT17, KMT18, KPPN0, KPPNI, KRPPN, &"       >> eqnfile
    print "  FCPPN, NCPPN, FPPN, KBPPN"                             >> eqnfile
    print "#ENDINLINE {above lines go into MODULE KPP_ROOT_Global}" >> eqnfile
    ERR=0
    EQUATION=0
    DEFVAR=0
    INLINE=0
    RO2CALCULATION=0
    RO2END=0
    # read in J-Substitution table into array JMatrix
    # first line is header, skip
    getline line < jsubstfile
    while ( (getline line < jsubstfile) > 0) {
        split(line, Jfeld);
        # make use of associative arrays in awk, index of array can be string
        # will be now MCM number, stored in second column
        JMatrix[Jfeld[2]] = Jfeld[1]
    }
}

{
    ERR=0
    PHOTO=0
    outfile=unkfile;
    # remove leading white space in case of preprocessor command:
    if (match($1,"#")) $1=$1
    # check if reached section of def var:
    if (match($0,"#DEFVAR")) {
        INLINE=0
        DEFVAR=1
        EQUATION=0
    }
    if (DEFVAR) {
        # change chloride from CL in MCM to Cl:
        if (match($0,/\s*CL\s*/)) {
            gsub("CL","Cl")
        }
        # change bromide from BR in MCM to Br:
        if (match($0,/\s*BR\s*/)) {
            gsub("BR","Br")
        }
        outfile=spcfile
        # only write if second element is ' = ':
        if (match($2,"=")) {
        } else {
            ERR=1
        }
    }
    # check if inline global or rconst has been reached:
    if (match($0,"#INLINE F90_GLOBAL")) {
        INLINE=1
        DEFVAR=0
        EQUATION=0
        ERR=0
    }
    if (match($0,"#INLINE F90_RCONST")) {
        INLINE=1
        DEFVAR=0
        EQUATION=0
        ERR=0

    }
    if (INLINE) {
        outfile=eqnfile;
        # find start of RO2 calculation
        if (match($0,"RO2 = &")) {
            sub("&","0.",$0)
            RO2CALCULATION=1
            print $0 >> outfile
        }

        # find lines of RO2 calculation
        if (RO2CALCULATION) {
            for (i=1; i<NF;i++) {
                if (match($i,/.{2}ind\_[0-z]+\)?/)) {
                    cindex=$i
                    # get rid of C( and ), an idea for combined removel would be appreciated
                    split(cindex,a,"(")
                    split(a[2],a,")")
                    cindex2 = "  IF ("a[1]">0)     RO2 = RO2 + C("a[1]")"
                    print cindex2 >> outfile
                }
            }
            if (!match($NF,"&") && !match($0,"RO2 = 0.")) {
                RO2CALCULATION=0
                RO2END=1
            }
        }
            
        
        # find reaction rate of HO2 + HO2 and exchange H2O with C(ind_H2O):
        if (match($1,"KMT06")) {
            sub("*H2O","*C(ind_H2O)",$0)
        }
        # find any other equations where H2O concentration is used:
        # Needs attention by user!
        if (match($0,/.*\={1}.*H2O.*/) && !match($1,"KMT06")) {
            print "Unknown water sink found."
            printf $0
        }
        # deleting lines which cannot be read by mecca:
        if (match($NF, "&")){
            multiline=1
        }
        if (multiline) {
            ERR=0
            multiline=0
        }
        if (match($2,"=")) {
            ERR=0
        }
        if (match($1,"=")) {
            ERR=1
        }
        if (match($1,"!")) {
            ERR=1
        }
    }
    # check if equation section has been reached:
    if (match($0,"#EQUATIONS")) {
        DEFVAR=0
        INLINE=0
        EQUATION=1
        outfile=eqnfile;
    }
    if (EQUATION) {
        if (match($1,"[0-9]+",Rnumber)) {
            outfile=eqnfile;
            # get reaction number from the number given by mcm
            gnum="<G_MCM_"Rnumber[0]">"
            # change chloride from CL in MCM to Cl:
            if (match($0,/\s*CL\s*/)) {
                gsub("CL","Cl")
            }
            # change bromide from BR in MCM to Br:
            if (match($0,/\s*BR\s*/)) {
                gsub("BR","Br")
            }
            # add DUMMY to equation to avoid multiple equation of same kind:
            if (match($0,/\yCClCONO3\y\s*\=\s*\yCHOCl\y\s*\+\s*\yCO\y\s*\+\s*\yHO2\y\s*\+\s*\yNO2\y/)) {
                print $0
                sub("=","= Dummy +")
            }
            # add significant sources of water:
            if (match($0,/OH\s*\+\s*HO2\s*\=/)) {
                sub("=","= H2O")
            }
            if (match($0,/OH\s*\+\s*CH4\s*\=/)) {
                sub("=","= H2O +")
            }
            if (match($0,/OH\s*\+\s*\HCHO\s*\=/)) {
                sub("=","= H2O +")
            }
            if (match($0,/OH\s*\+\s*\H2O2\s*\=/)) {
                sub("=","= H2O +")
            }
            # check for H2O sink in equations and add water as agent in chemical reaction:
            while (match($0,/\yH2O\y\*|\*\yH2O\y/,arr)) {
                sub(/\yH2O\y\*|\*\yH2O\y/,"")
                sub("=","+ H2O =")
            }
            for (i=1; i<NF;i++) {
                # add dummy if products are missing:
                if (match($0," = *:")) {
                    gsub("=","= Dummy",$0)
                    }
                # in case it's a photolysis
                if (match($i,/J\([0-9]*\)/)) {
                    PHOTO=1 # 1=This is a photolysis reaction.
                    gnum="<J_MCM_"Rnumber[0]">"
                }
            }
            # replace the inital eq. number:
            $1=gnum
            # append MCM as reference:
            sub(";","; {\\&MCM}",$0)
            if (PHOTO>0) {
                # photolysis reaction gets "hv" as reaction partner:
                sub("="," + hv =",$0)
                # and become J and Tr G as a eq attribute:
                sub(":",": {%TrGJ}",$0)
            } else {
                # normal gas phase reactions become Tr and G as a eq attribute:
                sub(":",": {%TrG}",$0)
            }
        }
    }
    if (match($0,"#INCLUDE atoms")) {
        ERR=1 # suppress standard output, atoms need manual modification
        print "#INCLUDE atoms" > spcfile
        print "Min;" > spcfile
        print "Pls;" > spcfile
    }
    # if no error, print the line:
    if ((ERR==0) && (RO2CALCULATION==0) && (RO2END==0)) {
        print $0 >> outfile
        # find start of INLINE F90_RCONST
        if (match($0,"#INLINE F90_RCONST")) {
            print "M = cair" >> outfile
            print "N2 = c(ind_N2)" >> outfile
            print "O2 = c(ind_O2)" >> outfile
            print "H2O = c(ind_H2O)" >> outfile
        }
        if (match($0,"#EQUATIONS")) {
            print "<G_N2> N2 = N2: {%TrG} 0. ;" >> outfile
            print "<G_O2> O2 = O2: {%TrG} 0. ;" >> outfile
            print "<G_H2O> H2O = H2O: {%TrG} 0. ;" >> outfile
        }
    } else if ((ERR==1) && (RO2CALCULATION==0) && (RO2END==0)) {  
        outfile=unkfile
        print $0 >> outfile
    }
    RO2END=0
}

END { }
