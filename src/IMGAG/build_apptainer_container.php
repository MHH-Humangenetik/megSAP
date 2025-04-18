<?php
/** 
	@page build_apptainer_container
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("build_apptainer_container", "Builds an apptainer container for a specified tool.");
$parser->addString("tool", "Tool to build a container for", false);
$parser->addString("tag", "Tag to build.", true, "master");
extract($parser->parse($argv));

//prevent undersore in tag
if (contains($tag, "_"))
{
	$tag = strtr($tag, "_", "-");
	print "Note: tool version/tag in container names must not contain underscore - using '{$tag}'\n";
}

//build
$sif = "{$tool}_{$tag}.sif";
$log = "{$tool}_{$tag}.log";
print "Building container {$sif} - in case of error see {$log}\n";

// write server name, user and date to logfile as header
$server = php_uname('n');
$user = getenv("USER") ?: getenv("LOGNAME");
$date = date("Y-m-d H:i:s");

$log_header = "Built on: {$server}\nBuilt by: {$user}\nBuild date: {$date}\n\n";
file_put_contents($log, $log_header);

exec2("apptainer build {$sif} data/tools/container_recipes/{$tool}_{$tag}.def >> $log 2>&1");
exec2("chmod 777 {$sif}");
print "Building container finished.\n";

if ($tool =="ngs-bits")
{
	//determine version
	list($stdout) = exec2("apptainer exec {$sif} MappingQC --version");
	$version = trim(strtr(implode("", $stdout), ["MappingQC"=>""]));
	print "ngs-bits version determined from container: {$version}\n";

	//determine final container name
	if($tag=="master")
	{
		$sif2 = "ngs-bits_master-".strtr($version, "_", "-").".sif";

		//move to container repo
		$container_repo = "/mnt/storage2/megSAP/tools/apptainer_container/";
		print "Deploying container {$sif2} to {$container_repo}\n";
		exec2("mv {$sif} {$container_repo}/{$sif2}");
		print "Deploying finished\n";
	}
}

?>