#
#  $Id: regres.pl,v 1.1 2009-04-21 01:21:42 eneil Exp $
#
#  This is the portion of a script this is included recursively
#

#
# Each test has a name, input file name, output file name,
# and the good output that is tested against.
#
@regres_list = 
    (
     {
	 exec_path   => "$top_builddir/mainprogs/main" , 
	 execute     => "chroma" , 
	 input       => "$test_dir/chroma/pbp/pbp.ini.xml" , 
	 output      => "pbp.candidate.xml",
	 metric      => "$test_dir/chroma/pbp/pbp.metric.xml" ,
	 controlfile => "$test_dir/chroma/pbp/pbp.out.xml" ,
     }
     );
