#
#  $Id: regres.pl,v 1.2 2007-11-10 21:33:20 edwards Exp $
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
	 input       => "$test_dir/chroma/hadron/gauge_transf/gauge_transf.ini.xml" , 
	 output      => "gauge_transf.candidate.xml",
	 metric      => "$test_dir/chroma/hadron/gauge_transf/gauge_transf.metric.xml" ,
	 controlfile => "$test_dir/chroma/hadron/gauge_transf/gauge_transf.out.xml" ,
     }
     );
