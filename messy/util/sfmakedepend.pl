#!/usr/bin/perl
#
# Fortran 90/77 dependency checker, 2002 version.
# http://www.arsc.edu/~kate/Perl/sfmakedepend
#
use 5.6.0;
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Getopt::Long;
use File::Copy;

my ($opt_help, $opt_man, $opt_file, @incdirs, @srcdirs, @hsrcdirs,
    $opt_modext, $opt_case, $compiler, $opt_depend, $drop_circ,
    @modsearch, $opt_modcase); #  mz_pj_20041124
our ($cpp, $add_ext, $mod_dir, $libdeps);
our $obj_ext = 'o';
our $ext = 'f';

# Parse the arguments, do the right thing for --help, --man.
Getopt::Long::Configure( "bundling" );
GetOptions("help" => \$opt_help,    "man" => \$opt_man,
        "file=s" => \$opt_file,     "I=s@" => \@incdirs,
	"srcdir=s@" => \@hsrcdirs,   "moddir=s" => \$mod_dir,
        "fext=s" => \$ext,          "objext=s" => \$obj_ext,
	"modext=s" => \$opt_modext, "addext=s" => \$add_ext,
	"case=s" => \$opt_case,     "compiler=s" => \$compiler,
	"depend=s" => \$opt_depend, "cpp" => \$cpp,
	"libdeps" => \$libdeps, "drop" => \$drop_circ,
	"modsearch=s@" => \@modsearch, "modcase=s" => \$opt_modcase ) #  mz_pj_20041124
  or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

our @suffixes = qw( .c .C .cc .cxx .cpp .f .F .f90 .F90 .f95 .F95 .for);
our @mod_files = ();
@incdirs = (@hsrcdirs, @incdirs) ; # mz_pj_20040420

my $mf = 'Makefile';
if ($opt_file) {
    $mf = $opt_file;
} elsif (-f "makefile") {
    $mf = 'makefile';
}
if ( !(-f $mf)) {
    system "touch $mf";
}

# extension used for compiler's private module information
our $mod_ext = "mod";
our $depend = "obj";
our $case = "lower";
our $obj_dep_flag;
my $ll = 79;             # maximum line length in Makefile
my $cray;
my $parasoft;
# Hernan hack my $nothing = "\n";

# Check the compiler first, then override the compiler-specific defaults
if ($compiler) {
    if ($compiler eq "crayold") {
        $cray = 1;
        $case = "lower";
        $depend = "obj";
        $obj_dep_flag = "-p";
    } elsif ($compiler eq "cray") {
        $case = "upper";
    } elsif ($compiler eq "parasoft") {
        $parasoft = 1;
        $case = "lower";
        $depend = "obj";
        $obj_dep_flag = "-module";
    } elsif ($compiler eq "sgiold") {
        $mod_ext = "kmo";
        $case = "upper";
    } elsif ($compiler eq "sgi" or $compiler eq "hp" or
	     $compiler eq "absoft") {
        $case = "upper";
    } elsif ($compiler eq "nag" or $compiler eq "ibm" or
	     $compiler eq "sun") {
        $case = "lower";
#  mz_pj_20080124+
    } elsif ($compiler eq "pathscale") {
	$case = "upper";
#  mz_pj_20080124-
    } else {
        warn "Unknown compiler: $compiler\n";
    }
}

$depend = $opt_depend if defined($opt_depend);
if ($depend eq "obj") {
    $drop_circ = 1;
}

$case = $opt_case if defined($opt_case);

# extension used for compiler's private module information
if ($opt_modext) {
    $mod_ext = $opt_modext;
}

# need to add some more dependencies so the .f file gets created
our $need_f;
if ($cpp and $depend eq "obj") {
    $need_f = 1;
}

my $mystring = '# DO NOT DELETE THIS LINE - used by make depend';

# Search for the includes in all the files
my $file;
my %sources;
# mz_pj_20040420+
my ( @flist, @hlist, $dir, @hsrcs, @srcs, $i) ;
if ($#ARGV + 1 > 0) {
    @flist = @ARGV
}
else {
    foreach $dir (@hsrcdirs,".") {
	opendir(DIR,$dir) ; @hlist =  readdir(DIR) ; closedir(DIR);
	@hsrcs = grep (/\.(f90|f|F|f95)$/, @hlist) ;
	### # exclude F*.f files
	### @srcs  = grep (!/F.*\.f/, @hsrcs) ;
	# exclude F*.f* files
	@srcs  = grep (!/F.*\.f.*/, @hsrcs) ;
	# prepend directory
        if ($dir ne ".") {
	    foreach $i (0 .. $#srcs) { $srcs[$i] = $dir."/".$srcs[$i]} } ;
	@flist = (@flist, @srcs) ;
    }
}
# mz_pj_20040420-

#foreach $file (@ARGV) { # mz_pj_20040420
foreach $file (@flist) { # mz_pj_20040420
    my $filesrc = findsrc($file);
    $sources{$file} = new Source_File($file, $filesrc);
    $sources{$file}->find_includes();
}

# Create new Makefile with new dependencies.

if ($mf ne "-") {
    copy( "$mf", "$mf.old");
    open(MFILE, "$mf.old") || die "can't read Makefile $mf.old: $!\n"; 
    open(NMFILE, "> $mf") || die "can't write $mf: $!\n";
    select(NMFILE);

    while (<MFILE>) {
        if (!/$mystring/) {
            print;
        } else {
            last;
        }
    }
    print $mystring, "\n";
}

# Now print out include and use dependencies in sorted order.
my $target;
foreach $target (sort keys(%sources)) {
    $sources{$target}->print();
# Hernan hack     print $nothing;
}

# print out module dependencies
if ( !( $cray || $parasoft) ) {
    my $modname;
    foreach $modname (sort keys(%main::mod_files)) {
        my ($name, $path, $suffix) = fileparse(
	    $sources{$main::mod_files{$modname}}->{'filepath'}, @suffixes);
	my $object = $path . $name . "." . $obj_ext;
	if (!( $drop_circ && lc($modname) eq lc($name)) ) {
	    $object =~ s#^\./##;
	    print "$modname.$mod_ext: $object\n";
	}
    }
}

#
# End of main
#

sub findfile {
# Let's see if we can find the included file. Look in current
# directory first, then in directories from -I arguments.
    my $file = shift;
    my ($found, $i, $filepath);

    $found = 0;

    if ( -f $file ) {
        $found = 1;
        $file =~ s#^\./##;          # convert ./foo.h to foo.h
        return $file;
    }
    foreach $i (0 .. $#incdirs) {
        $filepath = $incdirs[$i]."/".$file;
        if ( -f $filepath ) {
	    $found = 1;
            $filepath =~ s#^\./##;          # convert ./foo.h to foo.h
            return $filepath;
        }
    }
    if ( ! $found ) {
	$filepath = "";
    }
    $filepath;
}
#-----------------------------------------------------------------------
sub findsrc {
# Let's see if we can find the source-file.  Look in current
# directory first, then in directories from --srcdir arguments.
    my $src = shift;
    my($found, $i, $srcpath);

    $found = 0;

    if ( -f $src ) {
        $found = 1;
        $src =~ s#^\./##;          # convert ./foo.h to foo.h
        return $src;
    }
    foreach $i (0 .. $#srcdirs) {
        $srcpath = $srcdirs[$i]."/".$src;
        if ( -f $srcpath ) {
	    $found = 1;
            $srcpath =~ s#^\./##;          # convert ./foo.h to foo.h
            return $srcpath;
        }
    }
    if ( ! $found ) {
	$srcpath = "";
    }
    $srcpath;
}

#################################################################
package Source_File;

# hash containing names of included files
my %inc_files = ();

# Constructor
sub new {
    my $type = shift;
    my $filename = shift;
    my $path = shift;
    my $self = {};
    $self->{'Source_File'} = $filename;
    $self->{'filepath'} = $path;
    $self->{'includes'} = {};
    $self->{'uses'} = {};
    $self->{'modules'} = {};
    bless $self;
}

sub find_includes {
    my $self = shift;
    my $file = $self->{'filepath'};
    my($after, $filepath, $ref, $included, $use, $modname);
    local(*FILE);
    local($_);

    if (-f $file) {
        open(FILE, $file) || warn "Can't open $file: $!\n";
    } else {
	return;
    }
    while (<FILE>) {
	$included = "";
	$use = "";
        # look for Fortran style includes
        if (/^\s*include\s*['"]([^"']*)["']/i) {
            $included = $1;
	    $after = $';
	# C preprocessor style includes
        } elsif (/^#\s*include\s*["<]([^">]*)[">]/) {
            $included = $1;
	    $after = $';
        # Fortran 90 "use"
        } elsif (/^\s*use\s+(\w+)/i) {   
	    $use = $1;
# Change the case as needed - compiler dependent.
	    if ($main::case eq "upper") {
		$use = uc($use);
	    } elsif ($main::case eq "lower") {
		$use = lc($use);
	    }
	    $self->{'uses'}{$use} = 1;
	# Fortran 90 module
	} elsif (/^\s*module\s+(\w+)/i) {
	    $modname = $1;
	    if ($main::case eq "upper") {
		$modname = uc($modname);
	    } elsif ($main::case eq "lower") {
		$modname = lc($modname);
	    }
	    # Skip "module procedure" in interface blocks
	    next if (lc($modname) eq "procedure");

	    $main::mod_files{$modname} = $file;
	    $self->{'modules'}{$modname} = 1;
	}
	if ($included) {
	    # See if we've already searched this file
	    if ( $inc_files{$included} ) {
		$filepath = $inc_files{$included}{'filepath'};
	    } else {
                $filepath = main::findfile($included);
		$ref = new Source_File($included, $filepath);
                $inc_files{$included} = $ref;
# Search included file for includes
		$ref->find_includes();
	    }
            if ( $filepath ) {
		$self->{'includes'}{$included} = 1;
	    } else {
                if ($after !~ /bogus/i) {
#  mz_pj_20041124+
#		    warn "Can't find file: $included\n";
		    warn "... external include file: $included\n";
#  mz_pj_20041124-
		}
	    }
	}
    }
    close FILE;
}

sub print_includes {
    my $self = shift;
    my $target = shift;
    my $len_sum = shift;
    my $file;

    foreach $file (keys %{$self->{'includes'}}) {
	my $ref = $inc_files{$file};
	my $len = length($ref->{'filepath'}) + 1;
	if (($len_sum + $len > $ll) &&
	       (length($target) + 1 < $len_sum)) {
	    print "\n$target:";
	    $len_sum = length($target) + 1;
	}
	print " " . $ref->{'filepath'};
	$len_sum += $len;
	$len_sum = $ref->print_includes($target, $len_sum);
    }
    $len_sum;
}

# return list of modules used by included files
sub inc_mods {
    my $self = shift;
    my($file, $ref, $mod, @sub_list);
    my @list = ();

    foreach $mod (keys %{$self->{'uses'}}) {
	push(@list, $mod);
    }

    foreach $file (keys %{$self->{'includes'}}) {
	$ref = $inc_files{$file};
	@sub_list = $ref->inc_mods();
	@list = (@list, @sub_list);
    }
    @list;
}

# filenames containing the modules used by file and all its includes
sub find_mods {
    my $self = shift;
    my($modname, $file);
    my @module_files = ();
    my @mod_list = ();

# find modules used by include files
    if (%{$self->{'includes'}}) {
	foreach $file (keys %{$self->{'includes'}}) {
	    my $ref = $inc_files{$file};
	    my @list = $ref->inc_mods();
	    @mod_list = (@mod_list, @list);
	}
    }

# add them to the uses list (hash ensures uniqueness)
    foreach $modname (@mod_list) {
	$self->{'uses'}{$modname} = 1;
    }
    
# now find the filename that contains the module information
    foreach $modname (keys %{$self->{'uses'}}) {
	if ($main::depend eq "obj") {
	    if ($file = $main::mod_files{$modname}) {
                my $base = main::basename($file, @main::suffixes);
		$file = $base . "." . $main::obj_ext;
		push(@module_files, $file);
	    } else {
#  mz_pj_20041124+
#		warn "Don't know where module $modname lives.\n";
                ### The following is based on the assumption that
                ### f90 files are named as the modules they contain
		my $found = 0;
		my $foundmod;
		foreach $dir (@modsearch) {
		    if ( -f "$dir/$modname.f90" ) {
			$found = 1;
			$foundmod = "$dir/$modname.$mod_ext";
#  mz_pj_20080124+
		    if ($opt_modcase eq "upper") {
			$foundmod = "$dir/" . uc($modname) . ".$mod_ext";
		    } elsif ($main::case eq "lower") {
			$foundmod = "$dir/" . lc($modname) . ".$mod_ext";
		    } ;
#  mz_pj_20080124-
		    }
		}
		if ( $found ) {
		    warn "... dependency on $foundmod\n";
		    push(@module_files, "$foundmod");		    
		} else {
		    warn "... external module: $modname\n";
		}
#  mz_pj_20041124-
	    }
	} else {
	    if ($main::libdeps or defined($main::mod_files{$modname})) {
		$modname .= "." . $main::mod_ext;
		push(@module_files, $modname);
	    } else {
		warn "Couldn't locate source for module $modname\n";
	    }
	}
    }
    sort(@module_files);
}

sub print {
    my $self = shift;
    my $source = $self->{'Source_File'};
    my $compile_string = "\t" . '$(F90) $(F90FLAGS) -c';
    my($base, $object, $modname, $flag, $target, $ftarget);

    $base = main::basename($source, @main::suffixes);
    $target = $base . "." . $main::obj_ext;
    if ($main::cpp) {
	$ftarget = $base . "." . $main::ext;
    }

    $flag = $main::obj_dep_flag;
    
# print out "include" dependencies
    if (%{$self->{'includes'}}) {
	my $len_sum = length($target) + 1;
	if ($main::add_ext) {
	    print "$base.$main::add_ext ";
	    $len_sum += length($base) + length($main::add_ext) + 2;
	}
        print "$target:";
        $self->print_includes($target, $len_sum);
        print "\n";
        if ($main::cpp) {
	    $len_sum = length($ftarget) + 1;
            print "$ftarget:";
            $self->print_includes($ftarget, $len_sum);
            print "\n";
        }
    }

# clean out "use" of modules in own file
    my $mod;
    foreach $mod ( keys %{$self->{'uses'}} ) {
	if ( ${$self->{'modules'}}{$mod} ) {
	    delete ${$self->{'uses'}}{$mod};
	}
    }

# print out "use" dependencies
    if (%{$self->{'uses'}} || %{$self->{'includes'}}) {
	my @module_files = $self->find_mods();
#rch	my $len_sum = length($target) + 1;
#rch	print "$target:";
        my $len_sum = 0;
	my $file;
	foreach $file (@module_files) {
	    $file = $main::mod_dir . '/' . $file if $main::mod_dir;
            if( $len_sum < 1 ) {
 	        $len_sum = length($target) + 1;
 	        print "$target:";
            }
	    my $len = length($file) + 1;
	    if (($len_sum + $len > $ll) &&
		     (length($target) + 1 < $len_sum)) {
		print "\n$target:";
		$len_sum = length($target) + 1;
	    }
	    $len_sum += $len;
	    print " " . $file;
	}
	if ($main::need_f) {
	    my $len = length($ftarget) + 1;
	    if (($len_sum + $len > $ll) &&
		     (length($target) + 1 < $len_sum)) {
		print "\n$target:";
		$len_sum = length($target) + 1;
	    }
	    print " " . $ftarget if $len_sum;
	}
        print "\n" if $len_sum;
# extra Cray / Parasoft stuff
        if ($flag) {
	    print $compile_string;
	    foreach $file (@module_files) {
		print $flag . $file;
	    }
	    if ($main::cpp) {
	        print " " . $ftarget . "\n";
	    } else {
	        print " " . $source . "\n";
	    }
	}
    }
}

__END__

sfmakedepend - Fortran Dependency Checker

=head1 SYNOPSIS

sfmakedepend [--help] [--man] [--file=file]
             [-I dir] [--srcdir dir]  [--modsearch dir]
             [--modcase=upper|lower|asis]
             [--moddir dir]
             [--fext ext] [--objext ext] [--modext ext]
             [--addext ext] [--case=upper|lower|asis]
             [--compiler=crayold|cray|sgiold|sgi|nag|ibm|
                                 parasoft|hp|absoft|sun]
             [--depend=mod|obj] [--cpp] [--libdeps]
             [--drop]
             [file ...]

=head1 DESCRIPTION

This is a makedepend script for Fortran, including Fortran 90.
It searches for Fortran style includes, C preprocessor includes,
and module dependencies to the extent that I understand them.

Your files must have an extension listed in the @suffixes list
in the code. You might also want to modify $compile_string;
the compiler is called $(F90).

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<--help>

Print more details about the arguments.

=item I<--man>

Print a full man page.

=item I<--file file>

Change the name of the current Makefile (default is Makefile).
Use "--file -" to write to stdout.

=item I<-I dir>

Look in alternate directories for the include files. There can be
several "-I dir" options used at once. The current directory is
still searched first.

=item I<--srcdir dir>

Look in alternate directories for the source files, much like VPATH.
It can be used multiple times to search in more than one directory.

=item I<--modsearch dir>

Look in alternate directories for external module files, if
--depend=obj. It can be used multiple times to search in more than
one directory.

=item I<--modcase upper|lower|asis>

Make module file names of external (option --modsearch)
module file dependencies upper or lower case.

=item I<--moddir dir>

Tells sfmakedepend to prepend moddir to all module references. This
is required if you use a common module library directory for a
multi-directory project (+moddir= .. option on HP, eg.).

=item I<--fext>

This is used with the --cpp switch for compilers
which expect an extension other than .f on source files.  For
instance, one might choose to use "--fext f90".

=item I<--objext ext>

Tells sfmakedepend what extension to use for object files. The
default is "o", but "obj", for instance, is
appropriate on MS-DOG etc.

=item I<--modext ext>

Specifies the extension to use for Fortran 90 module files. The default
extension is "mod" since this seems to be an emerging standard. Let
me know if other compilers use a different filename for the module
information (keep that --compiler option up to date).

=item I<--addext ext>

Tells sfmakedepend to add targets with extension add_ext to the rules
for object files. For instance, to operate with (f77) ftnchek .prj
files, you could use

`--addext prj' to get rules like:
foo.prj foo.o: ... 

=item I<--case upper|lower|asis>

Controls case of module names when generating module file names. Only
relevant where the module file name is named after the module rather
than after the source file.

=item I<--compiler=crayold|cray|sgiold|sgi|nag|ibm|parasoft|hp|absoft|sun>

Controls the type of target compiler, setting the module name and other
options appropriately. The cray option assumes that FFLAGS includes
-e m for creating the .mod file, while crayold refers to the default of
including that information in the object file.

=item I<--depend=mod|obj>

Whether to use the module information file or the module object file
in dependencies.

=item I<--cpp>

There are times when one might choose to run .F files through cpp and
to keep the .f files (for the debugger or for a custom preprocessor).
In that case, make must be told about the cpp include dependencies of
the .f files. This option will provide those dependencies.

=item I<--libdeps>

Generate dependencies on modules for which source code is not
available. Otherwise a warning is issued, but the dependency is not
listed.

=item I<--drop>

Drop module dependencies (my_mod.mod: my_mod.o). This is also done when
--depend=obj.

=item I<[file ...]>

The list of source files to search for dependencies. If it is omitted,
all .f90 .f .F .f95 files in the --srcdir directories will be used.

=back

=head1 EXAMPLE

Search for include files in /usr/local/include:

   sfmakedepend --cpp --fext=f90 -I /usr/local/include *.F

Example usage in gnuMakefile:

   SRCDIRS= srcdir1 srcdir2 srcdir3 ...
   FSRCS0 := $(foreach DIR, . $(SRCDIRSH),$(wildcard $(DIR)/*.f))
   FSRCS  := $(sort $(notdir $(FSRCS0)))
   
   F_makedepend=sfmakedepend  --file - $(addprefix --srcdir ,$(SRCDIRSH)) \
                                      $(subst -I,-I ,$(Includes))
   depend $(MAKEFILE_INC):
             $(F_makedepend) $(FSRCS) >> $(MAKEFILE_INC)

   include $(MAKEFILE_INC)

=head1 AUTHOR

Kate Hedstrom (kate@arsc.edu)
   First Perl 5 Fortran 90 version November 1994.

=head1 CONTRIBUTORS

   Dave Love         (d.love@dl.ac.uk)
       Added the --objext and --addext options (1996).

   Patrick Jessee
       Added hp support (1997, now in --compiler option).

   Sergio Gelato     (gelato@tac.dk)
       Added the --compiler, --depend, --case options, and
       --(no)libdeps (1998).

   Tobias Buchal     (buchal@ifh.bau-verm.uni-karlsruhe.de)
       Added the --srcdir and --file - options (1999).

   Klaus Ramstock    (klaus@tdm234.el.utwente.nl)
       Added the --moddir option (1999).

   Sandra Schroedter (sandra@fsg-ship.de)
       Fix to preserve Makefile links (1999).

   Holger Bauer      (bauer@itsm.uni-stuttgart.de)
       Added the --drop option (2000).

   Patrick J?ckel    (joeckel@mpch-mainz.mpg.de)
       Made filelist optional (2004).
       Added the --modsearch option (2004).
       Added the --modcase option (2008).

Others I've doubtless forgotten.

=cut
#
# NOTES
#	This makedepend script is my first attempt at using perl 5
#	objects.  Therefore, it may not be the best example of how
#	to do this.  Also, it requires perl 5 and will die if you
#	to use it with an older perl.  The latest version is
#	available from:
#
#		http://www.arsc.edu/~kate/Perl/
#		ftp://ahab.rutgers.edu/pub/perl/sfmakedepend
#
#	Fortran 90 introduces some interesting dependencies.  Two
#	compilers I have access to (NAG f90 and IBM xlf) produce a
#	private "mod_name.mod" file if you define "module mod_name"
#	in your code.  This file is used by the compiler when you
#	use the module as a consistency check (type-safe).  On the
#	other hand, the Cray and Parasoft compilers store the module
#	information in the object file and then files which use the
#	modules need to be compiled with extra flags pointing to the
#	module object files.
#
#	This script assumes that all the files using and defining
#	modules are in the same directory and are all in the list of
#	files to be searched.  It seems that the industry has not
#	settled on a practical way to deal with a separate modules
#	directory, anyway.
#
#	I sometimes include non-existent files as a compile time
#	consistency check:
#
#	   #ifndef PLOTS
#	   #include "must_define_PLOTS"       /* bogus include */
#	   #endif
#
#	This program warns about include files it can't find, but
#	not if there is a "bogus" on the same line.
#
#     *	The f90 module dependencies can confuse some versions of
#	make, especially of the System V variety.  We use gnu
#	make because it has no problems with these dependencies.
#
# BUGS
#	It can sometimes produce duplicate dependencies.
#
#	It treats C preprocessor includes the same as Fortran
#	includes.  This can add unnecessary dependencies if you
#	use the -s flag and both kinds of includes.
#
#       Please let me know if you find any others.
#	Kate Hedstrom
#	kate@arsc.edu
