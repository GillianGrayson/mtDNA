use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 14971249; ($val <= 14971365); $val+=1)
{
	system "scancel $val";
}