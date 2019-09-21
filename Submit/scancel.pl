use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 12856000; ($val <= 12856390); $val+=1)
{
	system "scancel $val";
}